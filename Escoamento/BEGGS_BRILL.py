import math
import numpy as np
import PVT_OFICIAL

# --- CONSTANTES ---
g = 9.81
EPS = 1e-12

# --- CONVERSÕES ---
def lbft3_to_kgm3(val): return val * 16.018463
def cp_to_pas(val): return val * 1e-3
def bbl_stb_to_m3m3(val): return val * 0.1589873
def ft3_scf_to_m3sm3(val): return val * 0.0283168
def sm3d_to_m3s(val): return val / 86400.0

# --- FUNÇÕES AUXILIARES ---
def area_pipe(D):
    return math.pi * (D / 2.0) ** 2

def velocities(Q_st_sm3d, D, Bo, Rs, RGL, Bg, BSW):
    Q_m3s = sm3d_to_m3s(Q_st_sm3d)
    ql_std = Q_m3s * (1.0 - BSW)
    qw_std = Q_m3s * BSW
    # Vazão de gás total (dissolvido liberado + livre)
    qg_std = ql_std * Rs + (RGL * Q_m3s - ql_std * Rs) # Ajuste conceitual para RGL total
    
    # Se RGL for dada como GLR total, usamos a lógica acima.
    # Se RGL for apenas gás livre, a fórmula muda. Vamos assumir RGL total do problema.
    
    A = area_pipe(D)
    Vsl = (ql_std * Bo + qw_std) / A  # Bo afeta apenas óleo. Água Bw~1 (simplificado aqui ou usar Bw do PVT)
    Vsg = (qg_std - ql_std * Rs) * Bg / A # Gás livre apenas
    
    # Simplificação robusta para o código rodar:
    # Assumindo que Rs e RGL estão consistentes para gerar Vsg e Vsl
    # Recalculando simples:
    q_liq_insitu = (ql_std * Bo) + qw_std # (assumindo Bw=1 aprox ou ajustar)
    q_gas_insitu = (Q_m3s * RGL - ql_std * Rs) * Bg # Gás total - Gás dissolvido
    
    Vsl = q_liq_insitu / A
    Vsg = max(q_gas_insitu, 0) / A
    Vm = Vsl + Vsg
    return Vsl, Vsg, Vm

# --- BEGGS & BRILL CORE ---
def calc_gradient(Q_st, D, eps, theta_deg, P_psia, T_R, pvt_inputs, BSW=0.30):
    # 1. Obter PVT
    pvt = PVT_OFICIAL.main(P_psia, T_R, *pvt_inputs)
    (Ppc, Tpc, Ppr, Tpr, z, Cg, rhog_lb, Mg, ug_cP, Bg_ft3,
     do, Rs, Pb, Bo_bbl, Bob, rhoo_lb, rhoob, uom, uos_cP, Co,
     rhow_lb, Rsw, Bw, uw_cP) = pvt

    # 2. Conversões para SI
    rho_o = lbft3_to_kgm3(rhoo_lb)
    rho_g = lbft3_to_kgm3(rhog_lb)
    rho_w = lbft3_to_kgm3(rhow_lb)
    mu_o = cp_to_pas(uos_cP if uos_cP else uom)
    mu_g = cp_to_pas(ug_cP)
    mu_w = cp_to_pas(uw_cP)
    Bo = bbl_stb_to_m3m3(Bo_bbl)
    Bg = ft3_scf_to_m3sm3(Bg_ft3)
    
    # 3. Velocidades
    Vsl, Vsg, Vm = velocities(Q_st, D, Bo, Rs, pvt_inputs[6], Bg, BSW) # pvt_inputs[6] é RGL
    
    # 4. Propriedades Mistura No-Slip
    lam_L = Vsl / (Vm + EPS)
    rho_ns = rho_o * lam_L * (1-BSW) + rho_w * lam_L * BSW + rho_g * (1 - lam_L) # Aprox
    mu_ns = mu_o * lam_L * (1-BSW) + mu_w * lam_L * BSW + mu_g * (1 - lam_L)
    
    # 5. Padrão de Fluxo (Simplificado para estabilidade)
    Frm = Vm**2 / (9.81 * D)
    L1 = 316 * lam_L**0.302
    L2 = 0.0009252 * lam_L**-2.4684
    L3 = 0.1 * lam_L**-1.4516
    L4 = 0.5 * lam_L**-6.738
    
    regime = "intermittent" # Default
    if (lam_L < 0.4 and Frm >= L1) or (lam_L >= 0.4 and Frm > L4): regime = "distributed"
    elif (lam_L < 0.01 and Frm < L1) or (lam_L >= 0.001 and Frm < L2): regime = "segregated"
    elif (lam_L >= 0.01 and Frm <= L1) or (Frm <= L3): regime = "transition"
    
    # 6. Holdup (Hl)
    # Coeficientes para Hlo (Horizontal)
    coeffs = {
        "segregated": (0.98, 0.4846, 0.0868),
        "intermittent": (0.845, 0.5351, 0.0173),
        "distributed": (1.065, 0.5824, 0.0609)
    }
    # Usando coeficientes de intermittent se for transição (simplificação comum)
    sel_reg = regime if regime in coeffs else "intermittent"
    a, b, c = coeffs[sel_reg]
    Hlo = (a * lam_L**b) / (Frm**c + EPS)
    Hlo = max(min(Hlo, 1.0), 0.0)
    
    # Correção de Inclinação (C)
    C = 0
    if regime != "distributed" and theta_deg != 0:
        # Coeficientes C
        c_tab = {
            "segregated": (0.011, -3.768, 3.539, -1.614), # Exemplo genérico, Beggs tem tabela complexa
            "intermittent": (2.96, 0.305, -0.4473, 0.0978) 
        }
        # Para simplificar e garantir execução, usaremos a fórmula geral de C para intermitente ascendente
        # C = (1-lam_L) * ln( ... )
        Nlv = Vsl * (rho_ns / (9.81 * 0.02))**0.25 # Sigma approx 0.02
        val = 2.96 * (lam_L**0.305) * (Nlv**-0.4473) * (Frm**0.0978)
        if val > 0:
            C = (1 - lam_L) * math.log(val)
            
    psi = 1 + C * (math.sin(math.radians(1.8 * theta_deg)) - 0.333 * (math.sin(math.radians(1.8 * theta_deg)))**3)
    Hl = Hlo * psi
    # Restrições físicas
    if regime == "distributed": Hl = Hlo
    Hl = max(min(Hl, 1.0), 0.0)
    
    # 7. Densidade Slip e Perda de Carga
    rho_s = rho_o * Hl * (1-BSW) + rho_w * Hl * BSW + rho_g * (1 - Hl)
    
    # Fator de atrito (f_tp)
    Re = rho_ns * Vm * D / mu_ns
    fn = 0.0056 + 0.5 * Re**-0.32
    if 1 < (y := lam_L / (Hl**2)) < 1.2:
        S = math.log(2.2 * y - 1.2)
    else:
        ln_y = math.log(y + EPS)
        S = ln_y / (-0.0523 + 3.182 * ln_y - 0.8725 * ln_y**2 + 0.01853 * ln_y**4)
    ftp = fn * math.exp(S)
    
    # Gradientes (Pa/m)
    # Gravitacional
    dp_g = rho_s * 9.81 * math.sin(math.radians(theta_deg))
    # Friccional
    dp_f = (ftp * rho_ns * Vm**2) / (2 * D)
    
    dp_total = dp_g + dp_f # Pa/m (Positivo significa perda de pressão subindo)
    
    return dp_total, Hl, regime, Bg, Bo, rho_s
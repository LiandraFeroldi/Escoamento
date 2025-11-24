# BEGGS_BRILL.py  (versão final compatível com PVT_OFICIAL fornecido)
import math
import numpy as np
import PVT_OFICIAL
# constantes e utilitários
# -------------------------
g = 9.81
EPS = 1e-12

# conversões
def lbft3_to_kgm3(val): return val * 16.018463
def cp_to_pas(val): return val * 1e-3
def bbl_stb_to_m3m3(val): return val * 0.1589873
def ft3_scf_to_m3sm3(val): return val * 0.0283168
def sm3d_to_m3s(val): return val / 86400.0

def area_pipe(D):
    return math.pi * (D / 2.0) ** 2

# tabelas do método (ABC e CEFG)
ABC_TABLE = {
    "segregated":    (0.98,  0.4846, 0.0868),
    "intermittent":  (0.845, 0.5351, 0.0173),
    "distributed":   (1.065, 0.5824, 0.0609)
}

CEFG_TABLE_ASC = {
    "segregated":    (0.011, -3.768,  3.539,  -1.614),
    "intermittent":  (2.960, 0.305,  -0.4473, 0.0978),
    # distributed não é usado (C = 0)
}

# ----------------------------------------------------------------
# funções auxiliares Beggs & Brill (L_params, friction, regime)

def L_params(lambda_L):
    L1 = 316.0      * (lambda_L**0.302)
    L2 = 0.0009252  * (lambda_L**-2.4684)
    L3 = 0.1        * (lambda_L**-1.4516)
    L4 = 0.5        * (lambda_L**-6.738)
    return L1, L2, L3, L4

def friction_factor_hall(eps, D, Re):
    term = 2e4 * (eps / (D + EPS)) + 1e6 / (Re + EPS)
    return 0.0055 * (1+  (term ** (1/3)))

def froude_mixture_squared(Vm, D):
    return Vm**2 / (g * D + EPS)

def lambda_no_slip(Vsl, Vm):
    return max(0.0, min(Vsl / (Vm + EPS), 1.0))

def determine_regime(lambda_L, Frm2, L1, L2, L3, L4):
    if (lambda_L < 0.4 and Frm2 >= L1) or (lambda_L >= 0.4 and Frm2 > L4):
        return "distributed"
    if (lambda_L < 0.01 and Frm2 < L1) or (lambda_L >= 0.001 and Frm2 < L2):
        return "segregated"
    if (0.01 <= lambda_L < 0.4 and L3 <= Frm2 <= L1) or (lambda_L >= 0.4 and L3 <= Frm2 <= L4):
        return "transition"
    return "intermittent"

def compute_HLo_raw(lambda_L, Frm2, regime):
    a, b, c = ABC_TABLE[regime]
    return (a * (lambda_L**b)) / (Frm2**c + EPS)

def compute_C_parameter(lambda_L, Vsl, Frm2, rho_l, sigma_gl, regime, is_descending):
    # C = 0 para distributed ou trechos descendentes
    if regime == "distributed" or is_descending:
        return 0.0
    if regime not in CEFG_TABLE_ASC:
        return 0.0
    dprime, e, f, gexp = CEFG_TABLE_ASC[regime]
    denom = g * (sigma_gl + EPS)
    Nlv = Vsl * ( (rho_l / denom) ** 0.25 )
    inside = max(dprime * (lambda_L**e) * (Nlv**f) * (Frm2**gexp), EPS)
    C = (1 - lambda_L) * math.log(inside)
    return C

def psi_inclination(C, theta_deg):
    th = math.radians(theta_deg)
    return 1.0 + C * (math.sin(1.8 * th) - 0.333 * (math.sin(1.8 * th)**3))

# ----------------------------------------
# função de velocidades (robusta e clara)

def velocities_from_q_sm3d(Q_st_sm3d, D, Bo, Rs, RGL, Bg, BSW, Bw=1.0):
    Q_m3s = sm3d_to_m3s(Q_st_sm3d)
    q_oil_std = Q_m3s * (1.0 - BSW)
    q_water_std = Q_m3s * BSW

    # gás total no separador (std m3/s)
    gas_from_RGL = max(RGL * Q_m3s, 0.0)
    gas_dissolved = max(q_oil_std * Rs, 0.0)
    qg_free_std = max(gas_from_RGL - gas_dissolved, 0.0)

    q_liq_insitu = q_oil_std * Bo + q_water_std * Bw
    q_gas_insitu = qg_free_std * Bg

    A = area_pipe(D)
    Vsl = q_liq_insitu / (A + EPS)
    Vsg = q_gas_insitu / (A + EPS)
    Vm = Vsl + Vsg + EPS

    return Vsl, Vsg, Vm, q_liq_insitu, q_gas_insitu

# -------------------------------------------------------
# Função que o MAIN_SIMULATOR espera (calc_gradient)

def calc_gradient(Q_st_sm3d, D, eps, theta_deg, P_psia, T_R, pvt_inputs, BSW=0.30):
    """
    Interface compatível com MAIN_SIMULATOR:
      Q_st_sm3d : produção em sm3/d
      D         : diâmetro (m)
      eps       : rugosidade (m)
      theta_deg : ângulo trecho (graus)
      P_psia    : pressão absoluta em psia
      T_R       : temperatura em Rankine
      pvt_inputs: tupla com parâmetros que seu MAIN passa para PVT_OFICIAL
                 (o código tenta extrair RGL de pvt_inputs[5] mas tem fallback)
      BSW       : fração de água padrão
    Retorna:
      dp_total (Pa/m), Hl, regime, Bg (m3/std m3), Bo (m3/m3), rho_slip (kg/m3)
    """

    # -------------------------
    # 1) Chama PVT_OFICIAL com os inputs recebidos
    # -------------------------
    # Espera-se que MAIN forneça pvt_inputs apropriado para PVT_OFICIAL.main
    pvt = PVT_OFICIAL.main(P_psia, T_R, *pvt_inputs)

    # desempacotar conforme PVT_OFICIAL.main que você forneceu
    (Ppc, Tpc, Ppr, Tpr, z, Cg, rhog_lb, Mg, ug_cP, Bg_ft3,
     do, Rs, Pb, Bo_bbl, Bob, rhoo_lb, rhoob, uom, uos_cP, Co,
     rhow_lb, Rsw, Bw_pvt, uw_cP) = pvt

    # -------------------------
    # 2) conversões SI
    rho_o = lbft3_to_kgm3(rhoo_lb)
    rho_g = lbft3_to_kgm3(rhog_lb)
    rho_w = lbft3_to_kgm3(rhow_lb)
    mu_o = cp_to_pas(uos_cP if uos_cP else uom)
    mu_g = cp_to_pas(ug_cP)
    mu_w = cp_to_pas(uw_cP)
    Bo = bbl_stb_to_m3m3(Bo_bbl)
    Bg = ft3_scf_to_m3sm3(Bg_ft3)
    Bw = Bw_pvt if Bw_pvt is not None else 1.0

    # -------------------------
    # 3) Extrair RGL do pvt_inputs (robusto)
    RGL = None
    if isinstance(pvt_inputs, (list, tuple)) and len(pvt_inputs) > 5:
        # heurística: posição 5 costuma ser RGL nas versões anteriores do MAIN
        candidate = pvt_inputs[5]
        try:
            RGL = float(candidate)
        except:
            RGL = None

    # fallback: busca no pvt_inputs por valor plausível ( 0 <= RGL <= 1e6 )
    if RGL is None:
        for v in pvt_inputs:
            try:
                fv = float(v)
                if 0.0 <= fv <= 1e6:
                    # evita pegar dg ou TF (valores pequenos); preferir valores > 1
                    if fv > 1.0:
                        RGL = fv
                        break
            except:
                continue
    if RGL is None:
        # ultimo recurso, assume 0 (sem gás livre)
        RGL = 0.0

    # ------------------------------------------------
    # 4) Velocidades e vazões in-situ

    Vsl, Vsg, Vm, q_liq_insitu, q_gas_insitu = velocities_from_q_sm3d(
        Q_st_sm3d, D, Bo, Rs, RGL, Bg, BSW, Bw
    )

    # ------------------------------------------------
    # 5) No-slip mixture properties e Fr

    lambda_L = lambda_no_slip(Vsl, Vm)

    vol_oil = ( (Q_m3s := sm3d_to_m3s(Q_st_sm3d)) * (1 - BSW) * Bo ) / (q_liq_insitu + EPS)
    vol_water = 1.0 - vol_oil
    rho_liq_insitu = 1.0 / ( (vol_oil/(rho_o + EPS)) + (vol_water/(rho_w + EPS)) )

    # sanity fallback
    if not (50 < rho_liq_insitu < 2000):
        rho_liq_insitu = rho_o*(1-BSW) + rho_w*BSW

    rho_ns = rho_liq_insitu * lambda_L + rho_g * (1 - lambda_L)
    mu_ns = mu_o * lambda_L * (1-BSW) + mu_w * lambda_L * BSW + mu_g * (1 - lambda_L)

    Frm2 = froude_mixture_squared(Vm, D)
    L1, L2, L3, L4 = L_params(lambda_L)
    regime = determine_regime(lambda_L, Frm2, L1, L2, L3, L4)
    is_descending = (theta_deg < 0)

    # ------------------------------------------------
    # 6) HLo (raw) e aplicação da regra HLo > lambda_L

    if regime == "transition":
        HLo_seg_raw = compute_HLo_raw(lambda_L, Frm2, "segregated")
        HLo_int_raw = compute_HLo_raw(lambda_L, Frm2, "intermittent")
        HLo_seg = max(HLo_seg_raw, lambda_L + 1e-6)
        HLo_int = max(HLo_int_raw, lambda_L + 1e-6)
    else:
        HLo_raw = compute_HLo_raw(lambda_L, Frm2, regime)
        HLo = max(HLo_raw, lambda_L + 1e-6)

    # ------------------------------------------------
    # 7) Correção de inclinação C e psi; cálculo final HL
    # ------------------------------------------------
    # tensão superficial usada: tenta obter de PVT, fallback 0.02
    sigma_gl = 0.02

    if regime == "transition":
        C_seg = compute_C_parameter(lambda_L, Vsl, Frm2, rho_liq_insitu, sigma_gl, "segregated", False)
        psi_seg = psi_inclination(C_seg, theta_deg)
        HL_seg = HLo_seg * psi_seg

        C_int = compute_C_parameter(lambda_L, Vsl, Frm2, rho_liq_insitu, sigma_gl, "intermittent", False)
        psi_int = psi_inclination(C_int, theta_deg)
        HL_int = HLo_int * psi_int

        A = (L3 - Frm2) / (L3 - L2 + EPS)
        A = max(min(A, 1.0), 0.0)
        HL = max(min(A * HL_seg + (1 - A) * HL_int, 0.9999), 1e-6)
    else:
        if regime == "distributed":
            HL = HLo
        else:
            C = compute_C_parameter(lambda_L, Vsl, Frm2, rho_liq_insitu, sigma_gl, regime, is_descending)
            psi = psi_inclination(C, theta_deg)
            HL = max(min(HLo * psi, 0.9999), 1e-6)

    # ------------------------------------------------
    # 8) Densidade Slip e Re, fator de atrito

    rho_slip = rho_liq_insitu * HL + rho_g * (1 - HL)
    Re_ns = max(rho_ns * Vm * D / (mu_ns + EPS), 1e-9)
    f_n = friction_factor_hall(eps, D, Re_ns)

    # cálculo y e s
    y = lambda_L / (HL**2 + EPS)
    if 1.0 < y < 1.2:
        s = math.log(max(2.2 * y - 1.2, EPS))
    else:
        ly = math.log(max(y, EPS))
        den = -0.0523 + 3.182 * ly - 0.8725 * ly**2 + 0.01853 * ly**4
        s = (ly / den) if abs(den) > EPS else 0.0

    f_tp = f_n * math.exp(s)

    # perdas fricção e gravidade (Pa/m)
    dp_fric = (f_tp * rho_ns * Vm**2) / (2.0 * D)
    dp_grav = rho_slip * g * math.sin(math.radians(theta_deg))

    # ---------------
    # 9) termo de aceleração Ek e dp_total (Pa/m)
    P_abs_Pa = P_psia * 6894.76  # psia -> Pa
    Ek = (rho_slip * Vm * Vsg) / (P_abs_Pa + EPS)
    dp_total = (dp_fric + dp_grav) / max(1.0 - Ek, 1e-6)   # Pa/m, positivo = perda subindo

    # Retorna na ordem que o MAIN espera:
    # dp_total (Pa/m), HL, regime, Bg (m3/std per std), Bo (m3/m3), rho_slip (kg/m3)
    return dp_total, HL, regime, Bg, Bo, rho_slip

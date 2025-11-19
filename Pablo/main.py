# src/main.py
"""
Simulador inicial de perfil termo-hidráulico (versão estrutural).
Substitua as funções em calc_pvt(...) e beggs_brill(...) pelas equações completas do PDF.
"""
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from pathlib import Path

# -------------------------
# Config / entradas (A3)
# -------------------------
depth_total_m = 2100.0        # m (do enunciado)
Q_o_sc_m3day = 10000.0        # oil standard m3/day
Q_o_sc = Q_o_sc_m3day / 86400.0   # m3/s (std)
BSW = 0.30                    # fração de água no óleo (30%)
RGL = 150.0                   # sm3/sm3
API = 25.0
dg = 0.75                     # grav.espec. gás
P_res_bar = 550.0
T_res_C = 80.0
TEC_poco = 2.0  # W/mK (poço)
TEC_mar = 1.0   # W/mK (ambiente)
sigma_og = 0.00841
sigma_wg = 0.004
Dh = 0.10     # m (diâmetro hidráulico aproximado) -> ajuste conforme seu caso
DL = 10.0     # passo em m (escolha numérico)
nsteps = int(depth_total_m / DL)

# -------------------------
# Funções PVT (placeholder)
# -------------------------
def calc_pvt(P_bar, T_C, API, dg, RGL, BSW):
    """
    Entrada:
      P_bar, T_C
    Retorna dicionário com Bg, Bo, Bw, mu_g (cP), mu_o (cP), mu_w (cP), rho_g, rho_o, rho_w, Rs
    SUBSTITUIR por correlações do PDF (Papay, Standing, Lee, Beal, etc.)
    """
    # conversões (simples)
    P_psia = P_bar * 14.504
    T_R = (T_C * 9/5) + 491.67

    # --- simplificações temporárias (substituir) ---
    Z = max(0.1, 1.0 - 0.0001 * (P_psia/1000.0))   # aproximação de Z
    Bg = (14.7/ (P_psia)) * (T_R) * Z * 1e-3      # FVF gas (m3/sc) - escala arbitrária
    Bo = 1.0 + 1e-4*(100 - P_psia/100)            # FVF oil (placeholder)
    Bw = 1.0 + 1e-6 * P_psia                      # FVF water (placeholder)

    # viscosidades (cP)
    mu_g = 0.02  # cP (placeholder)
    mu_o = max(0.5, 10.0 * math.exp(-0.0005 * P_psia))  # cP (placeholder)
    mu_w = 0.5

    # densidades (kg/m3)
    rho_w = 1000.0
    rho_o = 1000.0 * (141.5/(API+131.5))  # mesma fórmula do pdf aproximada
    rho_g = dg * 1.225  # aproximado

    # solubilidade Rs placeholder
    Rs = max(0.0, 50.0 * (P_psia/1000.0))

    return {
        'Z': Z, 'Bg': Bg, 'Bo': Bo, 'Bw': Bw,
        'mu_g': mu_g, 'mu_o': mu_o, 'mu_w': mu_w,
        'rho_g': rho_g, 'rho_o': rho_o, 'rho_w': rho_w,
        'Rs': Rs
    }

# -------------------------
# Beggs & Brill simplificado (placeholder)
# -------------------------
def beggs_brill(props, Dh, Q_o_sc, Q_w_sc, Q_g_sc, theta_rad):
    """
    props: dicionario retornado por calc_pvt
    Retorna: dP_dL (bar/m), HL (liquid holdup)
    Substituir por implementação completa (L1..L4, padrão, tabela coeficientes).
    """
    # converte vazões padrão para in-situ (uso Bo, Bw, Bg)
    Q_o = Q_o_sc * props['Bo']   # m3/s
    Q_w = Q_w_sc * props['Bw']
    Q_g = Q_g_sc * props['Bg']   # nota: Bg placeholder

    A = math.pi * (Dh/2.0)**2
    V_L = (Q_o + Q_w) / A
    V_G = Q_g / A
    V_m = V_L + V_G + 1e-12

    lambda_L = V_L / V_m
    Frm = V_m**2 / (9.81 * Dh)

    # Holdup simplificado (empírico)
    HL0 = max(0.0, min(1.0, lambda_L * (1.0 + 0.1 * Frm)))
    # inclinação correction psi ~ 1 for now
    psi = 1.0
    HL = HL0 * psi

    # mixture density
    rho_mix = HL * props['rho_o'] + (1-HL) * props['rho_g']
    # mixture viscosity (linear rule-of-thumb)
    mu_mix = HL * props['mu_o'] + (1-HL) * props['mu_g']
    # Reynolds number
    Re = rho_mix * V_m * Dh / ( (mu_mix/1000.0) + 1e-12 )  # mu_mix cP -> Pa.s via /1000
    
    # friction factor (Blasius / laminar)
    if Re > 4000:
        f = 0.079 * Re**(-0.25)
    else:
        f = 64.0 / max(Re,1.0)

    # Darcy-Weisbach friction gradient (Pa/m)
    dP_dL_fric_Pa = f * rho_mix * V_m**2 / (2.0 * Dh)
    # gravity
    dP_dL_grav_Pa = rho_mix * 9.81 * math.sin(theta_rad)

    dP_dL_total_Pa = dP_dL_fric_Pa + dP_dL_grav_Pa
    dP_dL_bar = dP_dL_total_Pa / 1e5  # Pa -> bar

    return dP_dL_bar, HL

# -------------------------
# Simple thermal update
# -------------------------
def update_temperature(T_C, T_env_C, DL, Q_L, Dh, TEC=1.0):
    """
    modelo simples de aproximação: T_next = T_env + (T - T_env)*exp(-k*DL)
    k dependente de área e fluxo (simplificado).
    """
    A = math.pi*(Dh/2.0)**2
    v = Q_L / max(A, 1e-12)
    if v <= 0:
        return T_C
    tau = DL / v
    k = 0.1  # coef de resfriamento arbitrário -> ajustar via TEC
    T_next = T_env_C + (T_C - T_env_C) * math.exp(-k * DL / (tau + 1e-12))
    return T_next

# -------------------------
# Main loop
# -------------------------
def run_simulation():
    x = np.linspace(0, depth_total_m, nsteps+1)
    P_bar = np.zeros_like(x)
    T_C = np.zeros_like(x)
    HL = np.zeros_like(x)
    Bg = np.zeros_like(x)
    Bo = np.zeros_like(x)

    P_bar[0] = P_res_bar
    T_C[0] = T_res_C

    # vazões std:
    Q_o_sc_flow = Q_o_sc * (1.0 - BSW)
    Q_w_sc_flow = Q_o_sc * BSW
    Q_g_sc_flow = RGL * Q_o_sc  # sm3/s approximado -> usar conversões corretas

    # ambiente: assume leito marinho 4°C até manifold, depois varia (simplificado)
    T_env_C = 4.0

    theta_rad = math.radians(90.0 - 37.0) # exemplo (inclinação na seção) -> ajuste por trecho

    for i in range(len(x)-1):
        props = calc_pvt(P_bar[i], T_C[i], API, dg, RGL, BSW)
        dP_dL_bar, hl = beggs_brill(props, Dh, Q_o_sc_flow, Q_w_sc_flow, Q_g_sc_flow, theta_rad)

        # atualizar pressão
        P_bar[i+1] = max(0.001, P_bar[i] - dP_dL_bar * DL)

        # update temp
        Q_L_in_situ = (Q_o_sc_flow * props['Bo'] + Q_w_sc_flow * props['Bw'])
        T_C[i+1] = update_temperature(T_C[i], T_env_C, DL, Q_L_in_situ, Dh, TEC=TEC_poco)

        # salvar
        HL[i] = hl
        Bg[i] = props['Bg']
        Bo[i] = props['Bo']

    # último nó
    HL[-1] = HL[-2] if len(HL)>1 else HL[0]
    Bg[-1] = Bg[-2]
    Bo[-1] = Bo[-2]

    # montar dataframe e salvar
    df = pd.DataFrame({
        'x_m': x,
        'P_bar': P_bar,
        'T_C': T_C,
        'HL': HL,
        'Bg': Bg,
        'Bo': Bo
    })
    outdir = Path(__file__).resolve().parents[1] / 'data'
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / 'perfil.xlsx'
    df.to_excel(outpath, index=False)
    print("Salvo em:", outpath)

    # plots simples
    plt.figure()
    plt.plot(x, P_bar); plt.xlabel('x (m)'); plt.ylabel('P (bar)'); plt.title('Perfil de Pressão')
    plt.grid(True); plt.savefig(outdir / 'P_profile.png')

    plt.figure()
    plt.plot(x, T_C); plt.xlabel('x (m)'); plt.ylabel('T (°C)'); plt.title('Perfil de Temperatura')
    plt.grid(True); plt.savefig(outdir / 'T_profile.png')

    plt.figure()
    plt.plot(x, HL); plt.xlabel('x (m)'); plt.ylabel('HL'); plt.title('Holdup de Líquido')
    plt.grid(True); plt.savefig(outdir / 'HL_profile.png')

    plt.figure()
    plt.plot(x, Bg, label='Bg'); plt.plot(x, Bo, label='Bo'); plt.xlabel('x (m)'); plt.ylabel('Fatores de Formação'); plt.legend()
    plt.grid(True); plt.savefig(outdir / 'BF_profile.png')

    print("Gráficos salvos em:", outdir)

if __name__ == "__main__":
    run_simulation()

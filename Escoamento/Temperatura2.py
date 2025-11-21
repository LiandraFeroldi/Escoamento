# Temperatura2.py — versão reorganizada com perfil de temperatura por L
# Temperatura2.py — versão reorganizada com perfil de temperatura por L
import math
import numpy as np

import PVT_OFICIAL
from PVT_OFICIAL import main

# ---------------------------------------------------------------------------
# Funções auxiliares de unidade
# ---------------------------------------------------------------------------
def lbft3_to_kgm3(val_lbft3):
    return val_lbft3 * 16.018463

def pa_to_psia(p_pa):
    return p_pa / 6894.76

def psia_to_pa(p_psia):
    return p_psia * 6894.76

# ---------------------------------------------------------------------------
# Funções hidrodinâmicas / térmicas
# ---------------------------------------------------------------------------
def velocities_from_q(Q_st, D, Bo, Rs, RGL, Bg, BSW):
    A = np.pi * (D / 2) ** 2
    ql_std = Q_st * (1 - BSW)
    qw_std = Q_st * BSW
    qg_std = ql_std * Rs + RGL * qw_std

    Vsl = ql_std / (A * Bo)
    Vsg = qg_std / (A * Bg)
    Vm = Vsl + Vsg
    return Vsl, Vsg, Vm, ql_std, qw_std, qg_std


def massica(do_lbft3, ql_std):
    return ql_std * do_lbft3


def calor_especifico(do_lbft3, T_sup_C):
    Cp = ((2e-3) * T_sup_C - 1.429) * do_lbft3 + (2.67e-3) * T_sup_C + 3.049
    return Cp

# ---------------------------------------------------------------------------
# PERFIL DE TEMPERATURA POR L (Linear + interação térmica real)
# ---------------------------------------------------------------------------

def calcular_perfil_trecho(
    L_total,
    n_steps,
    T_start_R,
    T_end_ext_R,
    theta_rad,
    TEC,
    qm,
    Cp,
    P_start_psia,
    Q_st,
    D,
    RGL,
    BSW,
    other_pvt_inputs,
):
    dg, R, Mar, Api, TF_initial_F, S = other_pvt_inputs

    x = np.linspace(0, L_total, n_steps)
    dL = L_total / (n_steps - 1)

    # Perfil linear inicial T(L)
    T_linear = np.linspace(T_start_R, T_end_ext_R, n_steps)

    # Arrays resultado
    T_field = np.zeros(n_steps)
    P_field = np.zeros(n_steps)

    # Condições iniciais
    T_local = T_start_R
    P_local = P_start_psia

    for i in range(n_steps):
        # PVT
        TF_local = T_local - 459.67
        pvt = PVT_OFICIAL.main(P_local, T_local, dg, T_local, R, Mar, Api, RGL, TF_local, S)

        (
            Ppc,
            Tpc,
            Ppr,
            Tpr,
            z,
            Cg,
            rhog,
            Mg,
            ug,
            Bg,
            do,
            Rs,
            Pb,
            Bo,
            Bob,
            rhoo,
            rhoob,
            uom,
            uos,
            Co,
            rhow,
            Rsw,
            Bw,
            uw,
        ) = pvt

        # Densidades SI
        rhoo_si = lbft3_to_kgm3(rhoo)
        rhog_si = lbft3_to_kgm3(rhog)
        rhow_si = lbft3_to_kgm3(rhow)

        # Frações (estimativa)
        _, _, _, ql_std, qw_std, qg_std = velocities_from_q(Q_st, D, Bo, Rs, RGL, Bg, BSW)
        vol_tot = max(ql_std + qw_std + qg_std, 1e-12)
        fL = ql_std / vol_tot
        fW = qw_std / vol_tot
        fG = qg_std / vol_tot
        rho_mix = fL * rhoo_si + fW * rhow_si + fG * rhog_si

        # Modelo térmico real
        expo = math.exp(-(TEC * dL) / (qm * Cp)) if qm * Cp != 0 else 0
        T_env = T_linear[i]
        T_new = T_env - (T_env - T_local) * expo

        T_field[i] = T_new

        # Pressão hidrostática
        dz = dL * math.sin(theta_rad)
        P_new_pa = psia_to_pa(P_local) - rho_mix * 9.81 * dz
        P_local = pa_to_psia(P_new_pa)
        P_field[i] = P_local

        T_local = T_new

    return {"x": x, "T": T_field, "P": P_field}

# ---------------------------------------------------------------------------
# Função Geral — 3 Trechos
# ---------------------------------------------------------------------------

def calcular_perfil_completo(params):
    dg, R, Mar, Api, TF_init_F, S = params["other_pvt_inputs"]

    # Trecho 1
    t1 = calcular_perfil_trecho(
        params["TDVpoco"],
        params["n_steps"],
        params["T1_R"],
        params["T2_R"],
        params["theta_vertical"],
        params["TEC_poco"],
        params["qm"],
        params["Cp"],
        params["P1_psia"],
        params["Q_st"],
        params["D"],
        params["RGL"],
        params["BSW"],
        params["other_pvt_inputs"],
    )

    # Trecho 2
    t2 = calcular_perfil_trecho(
        params["L_bomba"],
        params["n_steps"],
        t1["T"][-1],
        params["T3_R"],
        params["theta_inclin"],
        params["TEC_poco"],
        params["qm"],
        params["Cp"],
        t1["P"][-1],
        params["Q_st"],
        params["D"],
        params["RGL"],
        params["BSW"],
        params["other_pvt_inputs"],
    )

    # Trecho 3
    t3 = calcular_perfil_trecho(
        params["L_manifold"],
        params["n_steps"],
        t2["T"][-1],
        params["T_surface_R"],
        0.0,
        params["TEC_marinho"],
        params["qm"],
        params["Cp"],
        t2["P"][-1],
        params["Q_st"],
        params["D"],
        params["RGL"],
        params["BSW"],
        params["other_pvt_inputs"],
    )

    return {"trecho1": t1, "trecho2": t2, "trecho3": t3}

# -------------------------------------------------------------
# Função de Plotagem dos Perfis T(L) e P(L) — convertendo unidades
# -------------------------------------------------------------
import matplotlib.pyplot as plt

def plotar_perfil(resultado):
    # Trecho 1
    x1 = resultado["trecho1"]["x"]
    T1_R = resultado["trecho1"]["T"]
    P1_psia = resultado["trecho1"]["P"]

    # Trecho 2
    x2 = resultado["trecho2"]["x"] + x1[-1]
    T2_R = resultado["trecho2"]["T"]
    P2_psia = resultado["trecho2"]["P"]

    # Trecho 3
    x3 = resultado["trecho3"]["x"] + x2[-1]
    T3_R = resultado["trecho3"]["T"]
    P3_psia = resultado["trecho3"]["P"]

    # Converte temperaturas (Rankine -> Celsius): Tc = (Tr - 491.67) * 5/9
    T1_C = (T1_R - 491.67) * 5.0/9.0
    T2_C = (T2_R - 491.67) * 5.0/9.0
    T3_C = (T3_R - 491.67) * 5.0/9.0

    # Converte pressão (psia -> bar): 1 bar = 14.503773773748295 psi
    psi_to_bar = 1.0 / 14.503773773748295
    P1_bar = P1_psia * psi_to_bar
    P2_bar = P2_psia * psi_to_bar
    P3_bar = P3_psia * psi_to_bar

    # --- Gráfico de Temperatura (°C) ---
    plt.figure(figsize=(12, 6))
    plt.plot(x1, T1_C, label="Trecho 1 — Vertical")
    plt.plot(x2, T2_C, label="Trecho 2 — Inclinado")
    plt.plot(x3, T3_C, label="Trecho 3 — Horizontal")
    plt.xlabel("Comprimento (m)")
    plt.ylabel("Temperatura (°C)")
    plt.title("Perfil de Temperatura ao longo do escoamento")
    plt.grid(True)
    plt.legend()
    plt.show()

    # --- Gráfico de Pressão (bar) ---
    plt.figure(figsize=(12, 6))
    plt.plot(x1, P1_bar, label="Trecho 1 — Vertical")
    plt.plot(x2, P2_bar, label="Trecho 2 — Inclinado")
    plt.plot(x3, P3_bar, label="Trecho 3 — Horizontal")
    plt.xlabel("Comprimento (m)")
    plt.ylabel("Pressão (bar)")
    plt.title("Perfil de Pressão ao longo do escoamento")
    plt.grid(True)
    plt.legend()
    plt.show()

# -------------------------------------------------------------
# Bloco MAIN — definir parâmetros antes de rodar
# -------------------------------------------------------------
if __name__ == "__main__":

    # -------- ENTRADAS DO SEU PROBLEMA --------
    dg = 0.75
    Api = 25
    R = 10.73
    Mar = 0.029
    RGL = 150
    S = 0

    T_sup_C = 17
    T_surface_R = T_sup_C * 9/5 + 491.67

    # Comprimentos
    TDVpoco = 450
    L_bomba = 1050
    L_manifold = 700

    # Temperaturas
    T1_R = 80 * (9/5) + 491.67     # fundo
    T2_R = 4 * (9/5) + 491.67      # manif.
    T3_R = 15 * (9/5) + 491.67     # subida

    # Pressões
    P1_psia = 550 * 14.504

    # Geometria e fluido
    D = 6 * 0.0254
    BSW = 0.25
    Q_st = 10000

    # Coeficientes térmicos
    TEC_poco = 2
    TEC_marinho = 1

    # Ângulos
    theta_vertical = math.radians(90)
    theta_inclin = math.radians(37)

    # Para o cálculo de qm e Cp precisamos rodar uma chamada previa de PVT
    TF_local = T1_R - 459.67
    PVT0 = PVT_OFICIAL.main(P1_psia, T1_R, dg, T1_R, R, Mar, Api, RGL, TF_local, S)
    (_,_,_,_,_,_,_,_,_,_, Bg, do, Rs, Pb, Bo, *_ ) = PVT0

    # Vazões padrão -> cálculo de qm e Cp
    _,_,_, ql_std, qw_std, qg_std = velocities_from_q(Q_st, D, Bo, Rs, RGL, Bg, BSW)

    qm = massica(do, ql_std)
    Cp = calor_especifico(do, T_sup_C)

    # -------- DEFINIR TODOS OS PARAMS --------
    params = {
        "TDVpoco": TDVpoco,
        "L_bomba": L_bomba,
        "L_manifold": L_manifold,

        "T1_R": T1_R,
        "T2_R": T2_R,
        "T3_R": T3_R,
        "T_surface_R": T_surface_R,

        "P1_psia": P1_psia,

        "TEC_poco": TEC_poco,
        "TEC_marinho": TEC_marinho,

        "theta_vertical": theta_vertical,
        "theta_inclin": theta_inclin,

        "Q_st": Q_st,
        "D": D,
        "RGL": RGL,
        "BSW": BSW,

        "qm": qm,
        "Cp": Cp,
        "n_steps": 500,

        "other_pvt_inputs": (dg, R, Mar, Api, TF_local, S),
    }

    # -------- EXECUTAR --------
    resultado = calcular_perfil_completo(params)
    plotar_perfil(resultado)

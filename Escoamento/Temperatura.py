import math
import numpy as np
from Escoamento import PVT_OFICIAL


def velocities_from_q(Q_st, D, Bo, Rs, RGL, Bg, BSW):
    A = np.pi * (D/2)**2
    ql_std = Q_st * (1 - BSW)
    qw_std = Q_st * BSW
    qg_std = ql_std * Rs + RGL * qw_std

    Vsl = ql_std / (A * Bo)
    Vsg = qg_std / (A * Bg)
    Vm  = Vsl + Vsg
    return Vsl, Vsg, Vm,ql_std,qw_std,qg_std  

def vazões_situ(ql_std,qw_std,qg_std,Rs,Rsw,Bg,Bo,Bw): 
    qg=(qg_std-(ql_std*Rs)-(qw_std*Rsw))*Bg
    ql=ql_std*Bo + qw_std*Bw
    return qg,ql

def massica(do,ql_std):
    qm = ql_std * do
    return qm

def calor_especifico(do, T_sup):
    Cp = ((2e-3) * T_sup - 1.429) * do + (2.67e-3) * T_sup + 3.049
    # Imprime valor de Cp (unidade: verificar se é J/(kg·K) ou outra)
    print(f'Calor Específico Óleo (Cp): {round(Cp, 5)}  (verifique unidade)')
    return Cp


def comprimento_poço(TDVpoco,theta1):
    Lpoco = TDVpoco /np.sin(theta1)
    print(f'Poço - Manifold: {round(Lpoco, 5)}m')
    return Lpoco

def malha(Lpoco):
    n=500
    deltaL=Lpoco/n
    print(f'Distância entre pontos Poço: {round(deltaL, 5)}m')
    return

def deltaT(T1,T2,n):
    dT =(T1 - T2) / n
    print(f'Distância entre pontos Temperatura: {round(dT, 5)}°R')


def calcular_temperaturas(
        TVDpoço, L_manifold, L_bomba,
        T1, T2, T3, Tpc,
        Cp, qm, g,
        theta1,
        TEC_poco, TEC_marinho,
        n=500
    ):

    # --- TRECHO INCLINADO ---
    x = np.linspace(0, TVDpoço, n)
    dT1 = (T1 - T2) / n
    T_novo = np.full_like(x, T1)
    T_old = [T1]
    T_L = []

    for i in range(n):
        T_novo[i] -= dT1 * i
        expo = np.exp(-(TEC_poco / (qm * Cp)) * TVDpoço)
        T_L_value = T_novo[i] - ((qm * g * np.sin(theta1)) / TEC_poco) * expo * \
                    (T_novo[i] - (qm * g * np.sin(theta1)) / TEC_poco - T_old[i])
        T_L.append(T_L_value)
        T_old.append(T_L_value)

    T_pr_inc = [T / Tpc for T in T_L]

    # --- TRECHO HORIZONTAL ---
    x_h = np.linspace(0, L_manifold, n)
    T_L2 = [T2 for _ in range(n)]
    T_pr_hor = [T2 / Tpc for _ in range(n)]

    # --- TRECHO VERTICAL ---
    y = np.linspace(TVDpoço, L_bomba, n)
    dT3 = (T3 - T2) / n
    T_novo = np.full_like(y, T2)
    T_old = [T2]
    T_L3 = []

    for i in range(n):
        T_novo[i] += dT3 * i
        expo = np.exp(-(TEC_marinho / (qm * Cp)) * L_bomba)
        T_L_value = T_novo[i] - ((qm * g) / TEC_marinho) * expo * \
                    (T_novo[i] - (qm * g) / TEC_marinho - T_old[i])
        T_L3.append(T_L_value)
        T_old.append(T_L_value)

    T_pr_vert = [T / Tpc for T in T_L3]

    return {
        "T_inclinado": T_L,
        "T_pr_inclinado": T_pr_inc,
        "T_horizontal": T_L2,
        "T_pr_horizontal": T_pr_hor,
        "T_vertical": T_L3,
        "T_pr_vertical": T_pr_vert
    }

if __name__ == "__main__":
    dg = 0.75 # Densidade relativa do gás
    Api = 25 # Grau API
    RGL = 150 # Razão gás óleo
    R = 10.73 # Constante universal dos gases [ft^3.psi/°R/lb.mol]
    Mar = 0.029 # Massa molecular do ar [kg/mol]
    P_sc = 14.7 # Pressão na condição padrão [psia]
    T_sc = 60 # Temperatura na condição padrão [°F]
    TEC_marinho = 1 #[w/mk]
    TEC_poco = 2
    T_sup = 17 #C

    L_bomba = 1050 #m
    L_manifold = 700 #m
    TDVpoco = 450

    d = 6 * 0.0254 # diâmetro máxima para obter surgencia
    e = 3 * 10**(-6)
    rho_w = 1000
    rho_ar = 1.225 #kg/m**3
    P1 = 550 * 14.504 #psi
    P1_bar = 550 #bar
    P3 = 5 * 14.504 #psi
    P3_bar = 5 #bar
    bsw = 0.25
    g = 9.81 #m/s**2
    T1 = 80 * (9/5) + 491.67
    T2_C = 4 #°C
    T1_F = T1 - 459.67
    T2 = 4 * (9/5) + 491.67
    T3 = 15 * (9/5) + 491.67
    sigma_wg = 0.004 #N/m
    sigma_og = 0.00841 #N/m
    theta1 = math.radians(90)
    theta_marinho = math.radians(37)
    theta3 = math.radians(90)
    qlsc_d=10000 #sm³/d
    S=0

    # --- CHAMAR PVT ---
TR = T_sup * 9/5 + 491.67     # Rankine
TF = T_sup * 9/5 + 32         # Fahrenheit

(Ppc, Tpc, Ppr, Tpr, z, Cg, rhog, Mg, ug, Bg,
 do, Rs, Pb, Bo, Bob, rhoo, rhoob, uom, uos, Co,
 rhow, Rsw, Bw, uw) = PVT_OFICIAL.main(
        P1, TR, dg, TR, R, Mar, Api, RGL, TF, S
)

# --- VELOCIDADES ---
Vsl, Vsg, Vm, ql_std, qw_std, qg_std = velocities_from_q(
    Q_st=qlsc_d, D=d,
    Bo=Bo, Rs=Rs, RGL=RGL,
    Bg=Bg, BSW=bsw
)

# --- VAZÕES IN-SITU ---
qg, ql = vazões_situ(ql_std, qw_std, qg_std, Rs, Rsw, Bg, Bo, Bw)

# --- MASSA E CP ---
qm = massica(do, ql_std)
Cp = calor_especifico(do, T_sup)

# --- TEMPERATURAS ---
resultado = calcular_temperaturas(
        TDVpoco, L_manifold, L_bomba,
        T1, T2, T3, Tpc,
        Cp, qm, g,
        theta1, theta_marinho,
        TEC_poco, TEC_marinho
)

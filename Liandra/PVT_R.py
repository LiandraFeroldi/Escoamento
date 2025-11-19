import math
import numpy as np

dg = 0.7 # Densidade relativa do gás
API = 30 # Grau API
RGL = 150 # Razão gás óleo
R = 10.73 # Constante universal dos gases [ft^3.psi/°R/lb.mol]
M_ar = 0.029 # Massa molecular do ar [kg/mol]
P_sc = 14.7 # Pressão na condição padrão [psia]
T_sc = 60 # Temperatura na condição padrão [°F]
TEC_marinho = 1 #[w/mk]
TEC_poco = 2
T_sup = 17
L_bomba = 1050 #m
L_manifold = 850 #m
Z_poco = 1
d = 6 * 0.0254 # diâmetro m
e = 3 * 10**(-6)
rho_w = 1000
rho_ar = 1.225 #kg/m**3
P1 = 550 * 14.504 #psi
P1_bar = 550 #bar
P3 = 15 * 14.504 #psi
P3_bar = 15 #bar
bsw = 0.30
g = 9.81 #m/s**2
T1 = 80 * (9/5) + 491.67
T2_C = 4 #°C
T1_F = T1 - 459.67
T2 = 4 * (9/5) + 491.67
T3 = 15 * (9/5) + 491.67
sigma_wg = 0.004 #N/m
sigma_og = 0.00841 #N/m
theta1= math.radians(90)
theta2 = math.radians(37)
theta3 = math.radians(90)

print(T1, T2, T3)
print(P1)
print(T1_F)

def calculate_pvt_res():
    # Área de passagem
    Ap = np.pi * (d / 2)**2
    print(f'Área de Passagem: {round(Ap, 5)} m^2')
    
    # Densidade do óleo in situ
    do = 141.5 / (API + 131.5)
    print(f'Densidade do óleo: {round(do, 5)} Kg/m^3')
    rho_o = rho_w * do
    print(f'Massa específica do óleo: {round(rho_o, 5)}')

    # Densidade do gás
    rho_g = dg * rho_ar
    print(rho_g, 'massa específica do gás')
    Q_I = 10000 / 86400 # [m^3/s]
    print(f'Vazão Volumétrica: {round(Q_I, 5)} kg/m^3')
    q_m = Q_I * do
    
    print(f'Vazão Mássica: {round(q_m, 5)} kg/m^3')
    Cp = ((2 * 10**-3) * T_sup - 1.429) * do + (2.67 * 10**-3) * T_sup + 3.049
    print(f'Calor Específico Óleo: {round(q_m, 5)} kg/m^3') # O print usa q_m, pode ser um typo no original
    return Ap, do,rho_g,q_m,Cp
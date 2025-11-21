import numpy as np
import math
import matplotlib.pyplot as plt
from PVT import (dados_oleo,dados_gas,vazao_liquido_std,vazao_std,vazao_insitu,propriedades_agua,propriedades_oleo,propriedades_gas,Razao_de_Solubilidade_agua,Volume_formação_agua,viscosidade_agua,massa_especifica_agua,viscosidade_oleo_saturado,obter_pressao_bolha,massa_especifica_oleo,fator_volume_formacao_oleo,razao_solubilidade_gas_oleo,propriedades_pseudocriticas,propriedades_pseudoreduzidas,fator_compressibilidade_papay,dados_gas,viscosidade_gas_lee,fator_formacao_gas)



def area(d):
     Ap = np.pi * (d / 2)**2
     return Ap








dg = 0.7 # Densidade relativa do gás
Api = 30 # Grau API
RGL = 250 # Razão gás óleo
R = 10.73 # Constante universal dos gases [ft^3.psi/°R/lb.mol]
Mar = 0.029 # Massa molecular do ar [kg/mol]
P_sc = 14.7 # Pressão na condição padrão [psia]
T_sc = 60 # Temperatura na condição padrão [°F]
TEC_marinho = 1 #[w/mk]
TEC_poco = 2
T_sup = 15
L_bomba = 1050 #m
L_manifold = 700 #m
Z_poco = 450
d = 3.5 * 0.0254 # diâmetro máxima para obter surgencia
e = 3 * 10**(-6)
rho_w = 1000
rho_ar = 1.225 #kg/m**3
P1 = 350 * 14.504 #psi
P1_bar = 350 #bar
P3 = 5 * 14.504 #psi
P3_bar = 5 #bar
bsw = 25
g = 9.81 #m/s**2
T1 = 80 * (9/5) + 491.67
T2_C = 4 #°C
T1_F = T1 - 459.67
T2 = 4 * (9/5) + 491.67
T3 = 15 * (9/5) + 491.67
sigma_wg = 0.004 #N/m
sigma_og = 0.00841 #N/m
theta1 = math.radians(90)
theta2 = math.radians(37)
theta3 = math.radians(90)
qlsc_d=10000 #sm³/d
S=0


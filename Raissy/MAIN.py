# MAIN
import pandas as pd
import math
import numpy as np
from temp import pvt_reser
from temp_2 import PVT

# Supondo que pvt_reser e PVT sejam módulos com as funções necessárias já definidas
T_L, T_L2, T_L3, T_pr_hor, T_pr_vert, T_pr_inc, e, d, Q_I, Bo, Bw, RGL, Rs, RSW, Bg, \
rho_w, rho_o, rho_g1, σwg, σog, μw_si, μo_si, g, teta_rad1, rho_g, µg, P1, L_poco, theta_3, \
P3, L_bomba, L_manifold, P1_bar, n, P_pc, P3_bar, do, API, bsw, Z_poco = \
    pvt_reser.calculate_pvt_reser()

P_L1, P_L2, P_L3, P_pr_hor, P_pr_inc, P_pr_vert, Z2_hor, Z3_vert, Z1_inc = \
    pvt_reser.calculate_bb(T_L, T_L2, T_L3, T_pr_hor, T_pr_vert, T_pr_inc, e, d, Q_I, Bo, Bw,
                         RGL, RS, RSW, Bg, rho_w, rho_o, σwg, σog, μw_si, μo_si, g, teta_rad1, rho_g, µg, P1,
                         L_poco, theta_3, P3, L_bomba, L_manifold, P1_bar, n, P_pc, P3_bar, Z_poco)

tabela = []

for i in range(len(T_L)):
    T1 = T_L[i]
    P1 = P_L1[i] * 14.504
    Z = Z1_inc[i]
    Ppr = P_pr_inc[i]
    Tpr = T_pr_inc[i]
    Ppr, Tpr, Z, Pb, µg, Bg, Eg, dg = PVT.fase_gas(T1, P1, Z, Ppr, Tpr, rho_g1)
    Rs, Βο, μο, Bt = PVT.fase_oleo(do, dg, API, T1, P1, RGL, Pb, Bg)
    rho_w, RSW, Bw, μw = PVT.fase_agua(T1, P1, bsw)
    tabela.append((P1, T1, Z, Pb, µg, Bg, Eg, Rs, Bo, Bt, μο, rho_w, RSW, Bw, μw))

for i in range(len(T_L)):
    T1 = T_L2[i]
    P1 = P_L2[i] * 14.504
    Z = Z2_hor[i]
    Ppr = P_pr_hor[i]
    Tpr = T_pr_hor[i]
    Ppr, Tpr, Z, Pb, µg, Bg, Eg, dg = PVT.fase_gas(T1, P1, Z, Ppr, Tpr, rho_g1)
    Rs, Βο, μο, Bt = PVT.fase_oleo(do, dg, API, T1, P1, RGL, Pb, Bg)
    rho_w, RSW, Bw, μw = PVT.fase_agua(T1, P1, bsw)
    tabela.append((P1, T1, Z, Pb, µg, Bg, Eg, Rs, Bo, Bt, μο, rho_w, RSW, Bw, μw))

for i in range(len(T_L)):
    T1 = T_L3[i]
    P1 = P_L3[i] * 14.504
    Z = Z3_vert[i]
    Ppr = P_pr_vert[i]
    Tpr = T_pr_vert[i]
    Ppr, Tpr, Z, Pb, µg, Bg, Eg, dg = PVT.fase_gas(T1, P1, Z, Ppr, Tpr, rho_g1)
    Rs, Βο, μο, Bt = PVT.fase_oleo(do, dg, API, T1, P1, RGL, Pb, Bg)
    rho_w, RSW, Bw, μw = PVT.fase_agua(T1, P1, bsw)
    tabela.append((P1, T1, Z, Pb, µg, Bg, Eg, Rs, Bo, Bt, μο, rho_w, RSW, Bw, μw))

# Create column names based on variable names (remove or modify as needed)
column_names = [
    "P", "T", "Z", "Pb", "µg", "Bg", "Eg", "Rs", "Bo", "Bt",
    "μο", "rho_w", "RSW", "Bw", "μω"
]

# Create a Pandas DataFrame from the list of lists
df = pd.DataFrame(tabela, columns=column_names)

# Export the DataFrame to an Excel file (replace 'tabela.xlsx' with your desired filename)
df.to_excel("tabela.xlsx", index=False)
print("Tabela successfully exported to tabela.xlsx!")
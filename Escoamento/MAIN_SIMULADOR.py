# Nome do Arquivo: MAIN_SIMULATOR.py
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Importação dos módulos auxiliares (certifique-se que estão na mesma pasta)
import PVT_OFICIAL
import BEGGS_BRILL as BB

# ============================================================================
# 1. DADOS DE ENTRADA (BASEADO NO PDF - GRUPO 2)
# ============================================================================
print("--- INICIANDO SIMULAÇÃO GRUPO 2 ---")

# --- Fluido ---
Q_std_d = 10000.0       # sm3/d (Meta de produção)
BSW = 0.30              # 30%
RGL = 150.0             # sm3/sm3
API = 25.0              # Grau API
dg = 0.75               # Gravidade específica gás

# Cálculo da densidade relativa do óleo
do = 141.5 / (API + 131.5) 

# --- Reservatório ---
P_res_bar = 380.0       # bar
T_res_C = 80.0          # °C

# --- Geometria do Duto ---
D_pol = 3.5             # Diâmetro
D_m = D_pol * 0.0254    # m
Rugosidade = 4.5e-5     # m
Area_Duto = math.pi * (D_m / 2)**2

# --- Propriedades Térmicas ---
TEC_poco = 2            # W/m.K
TEC_marinho = 1         # W/m.K

# --- Temperaturas Ambientais ---
T_fundo_mar_C = 4.0     # °C
T_superficie_C = 17.0   # °C

# ============================================================================
# 2. DEFINIÇÃO DA TOPOLOGIA
# ============================================================================
dL_step = 10.0 
sections = []

# Poço (2100->950)
for i in range(int((2100 - 950) / dL_step)):
    frac = i / int((2100 - 950) / dL_step)
    T_amb = 80.0 + (4.0 - 80.0) * frac
    sections.append({"theta": 90.0, "T_amb_C": T_amb, "TEC": TEC_poco, "dL": dL_step})

# Flowline (100m vert / sin(37))
L_flow = 100.0 / math.sin(math.radians(37.0))
for i in range(int(L_flow / dL_step)):
    sections.append({"theta": 37.0, "T_amb_C": T_fundo_mar_C, "TEC": TEC_marinho, "dL": dL_step})

# Riser (850m)
for i in range(int(850.0 / dL_step)):
    frac = i / int(850.0 / dL_step)
    T_amb = 4.0 + (17.0 - 4.0) * frac
    sections.append({"theta": 90.0, "T_amb_C": T_amb, "TEC": TEC_marinho, "dL": dL_step})

L_total_sistema = sum([s['dL'] for s in sections])

# ============================================================================
# 3. LOOP DE SIMULAÇÃO
# ============================================================================

P_atual_psia = P_res_bar * 14.5038
T_atual_R = (T_res_C * 9.0/5.0) + 491.67
L_acumulado = 0.0

TF_start = T_atual_R - 459.67
Pb_inicial_psia = PVT_OFICIAL.obter_pressao_bolha(RGL, dg, API, TF_start)

# Pegar densidades iniciais
pvt_init = PVT_OFICIAL.main(P_atual_psia, T_atual_R, dg, T_atual_R, 10.73, 28.96, API, RGL, TF_start, 0)
rho_g_init = pvt_init[6] * 16.018463
rho_o_init = pvt_init[15] * 16.018463

# Dicionário COMPLETO (com dp_fric e dp_grav inicializados)
dados_simulacao = {
    "L_m": [0.0], "P_bar": [P_res_bar], "T_C": [T_res_C], "Holdup": [0.0],
    "Regime": ["N/A"], "Bo": [0.0], "Bg": [0.0], "Pb_bar": [Pb_inicial_psia / 14.5038],
    "Vsg": [0.0], "dp_dL": [0.0], 
    "rho_o": [rho_o_init], "rho_g": [rho_g_init],
    "dp_fric": [0.0], "dp_grav": [0.0]  # <--- ADICIONADO
}

q_m3s = Q_std_d / 86400.0
qm_kg_s = q_m3s * 850.0 

print(f"Calculando {len(sections)} passos de simulação...")

for step in sections:
    dL = step['dL']
    theta = step['theta']
    T_amb_C = step['T_amb_C']
    TEC_linear = step['TEC']
    
    if P_atual_psia < 15: break # Trava de segurança

    pvt_inputs = (dg, T_atual_R, 10.73, 28.96, API, RGL, (T_atual_R - 459.67), 0)
    
    # Inicializa variáveis
    dp_total_Pa_m = 0.0
    dp_fric_Pa_m = 0.0
    dp_grav_Pa_m = 0.0
    rho_g_step = 0
    rho_o_step = 0

    try:
        # CAPTURA OS 8 VALORES DO BEGGS (incluindo fricção e gravidade)
        dp_total_Pa_m, Hl, regime, Bg, Bo, rho_mix_si, dp_fric_Pa_m, dp_grav_Pa_m = BB.calc_gradient(
            Q_std_d, D_m, Rugosidade, theta, P_atual_psia, T_atual_R, pvt_inputs, BSW
        )
        
        pvt_res = PVT_OFICIAL.main(P_atual_psia, T_atual_R, *pvt_inputs)
        rho_g_step = pvt_res[6] * 16.018463
        rho_o_step = pvt_res[15] * 16.018463
        
    except Exception as e:
        print(f"ERRO CRÍTICO em L={L_acumulado:.1f}m: {e}")
        break

    # Vsg Manual
    TF_atual = T_atual_R - 459.67
    Pb_local = PVT_OFICIAL.obter_pressao_bolha(RGL, dg, API, TF_atual)
    Rs_local = PVT_OFICIAL.razao_solubilidade_gas_oleo(P_atual_psia, dg, API, TF_atual, Pb_local)
    
    Q_oleo_std_s = (Q_std_d / 86400.0) * (1 - BSW)
    Gas_livre_std = max(0, (RGL - Rs_local) * Q_oleo_std_s)
    Q_gas_insitu = Gas_livre_std * Bg 
    Vsg_local = Q_gas_insitu / Area_Duto

    # Atualiza Pressão
    dp_psi_m = dp_total_Pa_m * 1.45038e-4
    P_new_psia = P_atual_psia - (dp_psi_m * dL)
    
    # Térmico
    T_amb_K = T_amb_C + 273.15
    T_old_K = (T_atual_R - 491.67) * 5.0/9.0 + 273.15
    theta_rad = math.radians(theta)
    Cp_kJ = ((2e-3) * T_amb_C - 1.429) * do + (2.67e-3) * T_amb_C + 3.049
    Cp_oleo = Cp_kJ * 1000.0 
    term_A = (qm_kg_s * 9.81 * math.sin(theta_rad)) / TEC_linear
    fator_exp = math.exp(-(TEC_linear * dL) / (qm_kg_s * Cp_oleo))
    T_new_K = (T_amb_K - term_A) - fator_exp * (T_amb_K - term_A - T_old_K)
    T_new_C = T_new_K - 273.15
    T_new_R = (T_new_C * 9.0/5.0) + 491.67
    
    # C. Pressão de Bolha
    TF_new = T_new_R - 459.67
    Pb_atual_psia = PVT_OFICIAL.obter_pressao_bolha(RGL, dg, API, TF_new)
    
    # D. ATUALIZAÇÃO E ARMAZENAMENTO
    L_acumulado += dL
    P_atual_psia = P_new_psia
    T_atual_R = T_new_R
    
    dados_simulacao["L_m"].append(L_acumulado)
    dados_simulacao["P_bar"].append(P_atual_psia / 14.5038)
    dados_simulacao["T_C"].append(T_new_C)
    dados_simulacao["Holdup"].append(Hl)
    dados_simulacao["Regime"].append(regime)
    dados_simulacao["Bo"].append(Bo)
    dados_simulacao["Bg"].append(Bg)
    dados_simulacao["Pb_bar"].append(Pb_atual_psia / 14.5038)
    dados_simulacao["Vsg"].append(Vsg_local) 
    dados_simulacao["dp_dL"].append(dp_total_Pa_m * 1e-5) 
    dados_simulacao["rho_o"].append(rho_o_step)
    dados_simulacao["rho_g"].append(rho_g_step)
    # SALVANDO COMPONENTES (CONVERTENDO PARA bar/m)
    dados_simulacao["dp_fric"].append(dp_fric_Pa_m * 1e-5) 
    dados_simulacao["dp_grav"].append(dp_grav_Pa_m * 1e-5)

# ============================================================================
# 4. RESULTADOS E PLOTAGEM
# ============================================================================
df = pd.DataFrame(dados_simulacao)

# --- PLOTAGEM (5 LINHAS PARA CABER O 9º GRÁFICO) ---
fig, axs = plt.subplots(5, 2, figsize=(14, 25)) # Aumentei altura para 25
fig.suptitle(f'Simulação de Escoamento - Grupo 2 (Q={Q_std_d} sm3/d)', fontsize=16)

# 1. Pressão
axs[0, 0].plot(df["L_m"], df["P_bar"], 'b-', lw=2); axs[0, 0].set_title('Pressão'); axs[0, 0].grid(True)
axs[0, 0].set_xlim(0, L_total_sistema)

# 2. Temperatura
axs[0, 1].plot(df["L_m"], df["T_C"], 'r-', lw=2); axs[0, 1].set_title('Temperatura'); axs[0, 1].grid(True)
axs[0, 1].set_xlim(0, L_total_sistema)

# 3. Holdup
axs[1, 0].plot(df["L_m"], df["Holdup"], 'g-', lw=2); axs[1, 0].set_title('Holdup'); axs[1, 0].grid(True)
axs[1, 0].set_xlim(0, L_total_sistema)

# 4. Fatores Volume
axs[1, 1].plot(df["L_m"], df["Bo"], 'b-', label='Bo'); axs[1, 1].set_title('Bo e Bg'); axs[1, 1].grid(True)
ax4_twin = axs[1, 1].twinx(); ax4_twin.plot(df["L_m"], df["Bg"], 'r--', label='Bg')
axs[1, 1].set_xlim(0, L_total_sistema)

# 5. Pb
axs[2, 0].plot(df["L_m"], df["Pb_bar"], 'purple', lw=2); axs[2, 0].set_title('Pb'); axs[2, 0].grid(True)
axs[2, 0].set_xlim(0, L_total_sistema)

# 6. Vsg
axs[2, 1].plot(df["L_m"], df["Vsg"], 'orange', lw=2); axs[2, 1].set_title('Vsg'); axs[2, 1].grid(True)
axs[2, 1].set_xlim(0, L_total_sistema)

# 7. Gradiente Total
axs[3, 0].plot(df["L_m"], df["dp_dL"], 'brown', lw=2); axs[3, 0].set_title('Gradiente Total (dP/dL)'); axs[3, 0].grid(True)
axs[3, 0].set_xlim(0, L_total_sistema)

# 8. Densidades
axs[3, 1].plot(df["L_m"], df["rho_o"], 'b-', label='Oleo'); axs[3, 1].set_title('Densidades'); axs[3, 1].grid(True)
ax8_twin = axs[3, 1].twinx(); ax8_twin.plot(df["L_m"], df["rho_g"], 'r--', label='Gas')
axs[3, 1].set_xlim(0, L_total_sistema)

# 9. NOVO GRÁFICO: Componentes de Perda
# Plotamos Fricção no eixo esquerdo e Gravidade no direito (pois gravidade pode ser bem maior no inicio)
axs[4, 0].plot(df["L_m"], df["dp_fric"], 'r-', label='Fricção')
axs[4, 0].set_xlabel('Comprimento (m)')
axs[4, 0].set_ylabel('Fricção (bar/m)', color='r')
axs[4, 0].tick_params(axis='y', labelcolor='r')
axs[4, 0].set_title('Componentes de Perda (Fricção vs Gravidade)')
axs[4, 0].grid(True)
axs[4, 0].set_xlim(0, L_total_sistema)

ax9_twin = axs[4, 0].twinx()
ax9_twin.plot(df["L_m"], df["dp_grav"], 'g--', label='Gravitacional')
ax9_twin.set_ylabel('Gravidade (bar/m)', color='g')
ax9_twin.tick_params(axis='y', labelcolor='g')

lines_9a, labels_9a = axs[4, 0].get_legend_handles_labels()
lines_9b, labels_9b = ax9_twin.get_legend_handles_labels()
axs[4, 0].legend(lines_9a + lines_9b, labels_9a + labels_9b, loc='best')

# Remove o gráfico vazio (10)
fig.delaxes(axs[4, 1])

plt.tight_layout()
plt.show()
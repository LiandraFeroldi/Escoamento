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
RGL = 300.0             # sm3/sm3
API = 25.0              # Grau API
dg = 0.75               # Gravidade específica gás

# Cálculo da densidade relativa do óleo
do = 141.5 / (API + 131.5) 

# --- Reservatório ---
P_res_bar = 300.0       # bar
T_res_C = 80.0          # °C

# --- Geometria do Duto ---
D_pol = 4             # Diâmetro
D_m = D_pol * 0.0254    # m
Rugosidade = 4.5e-5     # m
Area_Duto = math.pi * (D_m / 2)**2

# --- Propriedades Térmicas ---
TEC_poco = 10            # W/m.K
TEC_marinho =   25      # W/m.K

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
rho_mix_init = rho_o_init * (1-BSW) + 1000 * BSW

# Viscosidades iniciais (Estimativa para plot)
mu_o_init = (pvt_init[18] if pvt_init[18] else pvt_init[17]) * 0.001 # cP -> Pa.s
mu_g_init = pvt_init[8] * 0.001
mu_mix_init = mu_o_init

# Dicionário COMPLETO
dados_simulacao = {
    "L_m": [0.0], "P_bar": [P_res_bar], "T_C": [T_res_C], "Holdup": [0.0],
    "Regime": ["N/A"], "Bo": [0.0], "Bg": [0.0], "Pb_bar": [Pb_inicial_psia / 14.5038],
    "Vsg": [0.0], "Vsl": [0.0], "Vm": [0.0], 
    "dp_dL": [0.0], 
    "rho_o": [rho_o_init], "rho_g": [rho_g_init],
    "rho_mix": [rho_mix_init], 
    "rho_ns": [rho_mix_init], 
    "dp_fric": [0.0], "dp_grav": [0.0],
    "mu_o": [mu_o_init], "mu_g": [mu_g_init], "mu_mix": [mu_mix_init] # NOVAS COLUNAS
}

q_m3s = Q_std_d / 86400.0
qm_kg_s = q_m3s * 850.0 

print(f"Calculando {len(sections)} passos de simulação...")

for step in sections:
    dL = step['dL']
    theta = step['theta']
    T_amb_C = step['T_amb_C']
    TEC_linear = step['TEC']
    
    if P_atual_psia < 15: break 

    pvt_inputs = (dg, T_atual_R, 10.73, 28.96, API, RGL, (T_atual_R - 459.67), 0)
    
    # Inicializa variáveis
    dp_total_Pa_m = 0.0
    dp_fric_Pa_m = 0.0
    dp_grav_Pa_m = 0.0
    rho_mix_si = 0.0
    rho_ns_step = 0.0
    rho_g_step = 0.0
    rho_o_step = 0.0
    mu_o_val = 0.0; mu_g_val = 0.0; mu_mix_val = 0.0

    try:
        # CAPTURA OS 9 VALORES DO BEGGS
        dp_total_Pa_m, Hl, regime, Bg, Bo, rho_mix_si, dp_fric_Pa_m, dp_grav_Pa_m, rho_ns_step = BB.calc_gradient(
            Q_std_d, D_m, Rugosidade, theta, P_atual_psia, T_atual_R, pvt_inputs, BSW
        )
        
        pvt_res = PVT_OFICIAL.main(P_atual_psia, T_atual_R, *pvt_inputs)
        rho_g_step = pvt_res[6] * 16.018463
        rho_o_step = pvt_res[15] * 16.018463
        
        # --- CÁLCULO MANUAL DAS VISCOSIDADES PARA O GRÁFICO ---
        # (Já que não estamos pegando do BB, pegamos do PVT aqui)
        mu_o_cP = pvt_res[18] if pvt_res[18] else pvt_res[17] # uos ou uom
        mu_g_cP = pvt_res[8] # ug
        mu_w_cP = pvt_res[23] # uw
        
        mu_o_val = mu_o_cP * 0.001 # Pa.s
        mu_g_val = mu_g_cP * 0.001 # Pa.s
        mu_w_val = mu_w_cP * 0.001 # Pa.s
        
        # Viscosidade da Mistura No-Slip (Aprox)
        # mu_ns = mu_o * lam + mu_g * (1-lam)
        # Precisamos calcular lambda
        Pb_loc = PVT_OFICIAL.obter_pressao_bolha(RGL, dg, API, T_atual_R - 459.67)
        Rs_loc = PVT_OFICIAL.razao_solubilidade_gas_oleo(P_atual_psia, dg, API, T_atual_R - 459.67, Pb_loc)
        q_oil = (Q_std_d/86400)*(1-BSW)
        q_gas = max(0, (RGL - Rs_loc)*q_oil) * Bg
        vsg_tmp = q_gas / Area_Duto
        vsl_tmp = ((q_oil * Bo) + (Q_std_d/86400)*BSW) / Area_Duto
        lam_tmp = vsl_tmp / (vsl_tmp + vsg_tmp + 1e-9)
        
        mu_liq_mix = mu_o_val*(1-BSW) + mu_w_val*BSW
        mu_mix_val = mu_liq_mix * lam_tmp + mu_g_val * (1 - lam_tmp)
        
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

    # Líquido (Vsl)
    Q_liq_insitu = (Q_oleo_std_s * Bo) + (Q_std_d / 86400.0) * BSW 
    Vsl_local = Q_liq_insitu / Area_Duto

    # Mistura (Vm)
    Vm_local = Vsl_local + Vsg_local

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
    
    # Pb
    TF_new = T_new_R - 459.67
    Pb_atual_psia = PVT_OFICIAL.obter_pressao_bolha(RGL, dg, API, TF_new)
    
    # Armazenamento
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
    dados_simulacao["Vsl"].append(Vsl_local)
    dados_simulacao["Vm"].append(Vm_local) 
    dados_simulacao["dp_dL"].append(dp_total_Pa_m * 1e-5) 
    dados_simulacao["rho_o"].append(rho_o_step)
    dados_simulacao["rho_g"].append(rho_g_step)
    dados_simulacao["rho_mix"].append(rho_mix_si)
    dados_simulacao["rho_ns"].append(rho_ns_step)
    dados_simulacao["dp_fric"].append(dp_fric_Pa_m * 1e-5) 
    dados_simulacao["dp_grav"].append(dp_grav_Pa_m * 1e-5)
    dados_simulacao["mu_o"].append(mu_o_val) 
    dados_simulacao["mu_g"].append(mu_g_val) 
    dados_simulacao["mu_mix"].append(mu_mix_val)
# ============================================================================
# 4. RESULTADOS E PLOTAGEM - GRÁFICOS SEPARADOS
# ============================================================================

df = pd.DataFrame(dados_simulacao)

print("\n" + "="*60)
print(f" RELATÓRIO DE INTERPRETAÇÃO (L = {L_acumulado:.1f} m)")
if len(df) > 1:
    P_chegada = df["P_bar"].iloc[-1]
    if P_chegada >= 15.0:
        print(f"   >>> STATUS: SUCESSO! O poço produz.")
    else:
        print(f"   >>> STATUS: FALHA. Pressão insuficiente.")
print("="*60 + "\n")

try:
    df.to_excel("Relatorio_Grupo2.xlsx", index=False)
    print("Arquivo 'Relatorio_Grupo2.xlsx' gerado com sucesso.")
except Exception as e:
    print(f"Erro ao salvar Excel: {e}")

# ============================================================================
# GRÁFICO 1: PRESSÃO
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["P_bar"], 'b-', lw=2, label='Pressão')
plt.plot(df["L_m"], df["Pb_bar"], 'k--', lw=1.5, label='Pressão de Bolha')
plt.axhline(15, color='r', linestyle='--', label='Pressão Separador')
plt.title('Pressão ao Longo da Linha (bar)')
plt.xlabel('Comprimento (m)'); plt.ylabel('Pressão (bar)')
plt.grid(True); plt.legend()
plt.show()

# ============================================================================
# GRÁFICO 2: TEMPERATURA
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["T_C"], 'r-', lw=2)
plt.title('Temperatura (°C)')
plt.xlabel('Comprimento (m)')
plt.ylabel('T (°C)')
plt.grid(True)
plt.show()

# ============================================================================
# GRÁFICO 3: HOLDUP
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["Holdup"], 'g-', lw=2)
plt.title('Holdup Líquido')
plt.xlabel('Comprimento (m)'); plt.ylabel('Holdup (-)')
plt.grid(True)
plt.show()

# ============================================================================
# GRÁFICO 4: FATORES DE VOLUME (Bo e Bg)
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["Bo"], 'b-', label='Bo')
plt.plot(df["L_m"], df["Bg"], 'r--', label='Bg')
plt.title('Fatores de Volume')
plt.xlabel('Comprimento (m)')
plt.grid(True)
plt.legend()
plt.show()

# ============================================================================
# GRÁFICO 5: PRESSÃO DE BOLHA
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["Pb_bar"], 'purple', lw=2)
plt.title('Pressão de Bolha (bar)')
plt.xlabel('Comprimento (m)'); plt.ylabel('Pb (bar)')
plt.grid(True)
plt.show()

# ============================================================================
# GRÁFICO 6: VELOCIDADES SUPERFICIAIS
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["Vsg"], 'orange', lw=2, label='Vsg')
plt.plot(df["L_m"], df["Vsl"], 'c-', lw=2, label='Vsl')
plt.plot(df["L_m"], df["Vm"], 'k--', lw=2, label='Vm')
plt.title('Velocidades Superficiais (m/s)')
plt.xlabel('Comprimento (m)'); plt.ylabel('Velocidade (m/s)')
plt.grid(True); plt.legend()
plt.show()

# ============================================================================
# GRÁFICO 7: GRADIENTE TOTAL
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["dp_dL"], 'brown', lw=2)
plt.title('Gradiente Total (bar/m)')
plt.xlabel('Comprimento (m)'); plt.ylabel('dp/dL (bar/m)')
plt.grid(True)
plt.show()

# ============================================================================
# GRÁFICO 8: DENSIDADES
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["rho_o"], 'b-', label='Óleo')
plt.plot(df["L_m"], df["rho_mix"], 'k-', linewidth=2, label='Mistura')
plt.plot(df["L_m"], df["rho_g"], 'r--', label='Gás')
plt.title('Densidades (kg/m³)')
plt.xlabel('Comprimento (m)')
plt.ylabel('Densidade (kg/m³)')
plt.grid(True); plt.legend()
plt.show()

# ============================================================================
# GRÁFICO 9: PERDAS (FRICÇÃO E GRAVIDADE)
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["dp_fric"], 'r-', label='Perda por Fricção')
plt.plot(df["L_m"], df["dp_grav"], 'g--', label='Perda por Gravidade')
plt.title('Perdas por Fricção e Gravidade (bar/m)')
plt.xlabel('Comprimento (m)')
plt.ylabel('Perda (bar/m)')
plt.grid(True); plt.legend()
plt.show()

# ============================================================================
# GRÁFICO 10: VISCOSIDADES
# ============================================================================
plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["mu_o"], 'b-', label='Óleo')
plt.plot(df["L_m"], df["mu_mix"], 'k-', linewidth=2, label='Mistura')

plt.title('Viscosidades (Pa.s)')
plt.xlabel('Comprimento (m)')
plt.ylabel('Viscosidade (Pa.s)')
plt.grid(True)
plt.legend()

plt.figure(figsize=(10,5))
plt.plot(df["L_m"], df["mu_g"], 'r--', label='Gás')
plt.title('Viscosidade do Gás (Pa.s)')
plt.xlabel('Comprimento (m)')
plt.ylabel('Viscosidade (Pa.s)')
plt.grid(True)
plt.legend()
plt.show()

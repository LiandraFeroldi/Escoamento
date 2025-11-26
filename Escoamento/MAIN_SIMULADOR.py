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
RGL = 200.0             # sm3/sm3
API = 25.0              # Grau API
dg = 0.75               # Gravidade específica gás

# Cálculo da densidade relativa do óleo (necessário para o Cp das meninas)
do = 141.5 / (API + 131.5) 

# --- Reservatório (Condição de Contorno Inicial - Fundo) ---
P_res_bar = 380.0       # bar
T_res_C = 80.0          # °C

# --- Geometria do Duto ---
D_pol = 3.5             # Diâmetro (suposição padrão, ajuste se necessário)
D_m = D_pol * 0.0254    # m
Rugosidade = 4.5e-5     # m (Aço carbono típico)
Area_Duto = math.pi * (D_m / 2)**2

# --- Propriedades Térmicas ---
# O Cp será calculado dinamicamente no loop
TEC_poco = 2            # W/m.K (Linear)
TEC_marinho = 1         # W/m.K (Linear)

# --- Temperaturas Ambientais (Condições de Contorno Externas) ---
T_fundo_mar_C = 4.0     # °C
T_superficie_C = 17.0   # °C

# ============================================================================
# 2. DEFINIÇÃO DA TOPOLOGIA (GEOMETRIA DO SISTEMA)
# ============================================================================
# Discretização (tamanho do passo de cálculo)
dL_step = 10.0 # metros (quanto menor, mais preciso)

sections = []

# --- TRECHO 1: POÇO (Reservatório -> ANM) ---
L_poco = 2100 - 950
n_steps1 = int(L_poco / dL_step)

for i in range(n_steps1):
    z_atual = 2100 - (i * dL_step)
    frac = i / n_steps1 
    T_amb = 80.0 + (4.0 - 80.0) * frac
    sections.append({"theta": 90.0, "T_amb_C": T_amb, "TEC": TEC_poco, "dL": dL_step})

# --- TRECHO 2: FLOWLINE (ANM -> Manifold) ---
L_flowline = 100.0 / math.sin(math.radians(37.0))
n_steps2 = int(L_flowline / dL_step)

for i in range(n_steps2):
    sections.append({"theta": 37.0, "T_amb_C": T_fundo_mar_C, "TEC": TEC_marinho, "dL": dL_step})

# --- TRECHO 3: RISER (Manifold -> Plataforma) ---
L_riser = 850.0
n_steps3 = int(L_riser / dL_step)

for i in range(n_steps3):
    z_atual = 850 - (i * dL_step)
    frac = i / n_steps3
    T_amb = 4.0 + (17.0 - 4.0) * frac
    sections.append({"theta": 90.0, "T_amb_C": T_amb, "TEC": TEC_marinho, "dL": dL_step})

# CALCULA O COMPRIMENTO TOTAL PREVISTO (Para fixar o gráfico)
L_total_sistema = sum([s['dL'] for s in sections])

# ============================================================================
# 3. LOOP DE SIMULAÇÃO (MARCHA PxT)
# ============================================================================

# Conversões Iniciais
P_atual_psia = P_res_bar * 14.5038
T_atual_R = (T_res_C * 9.0/5.0) + 491.67
L_acumulado = 0.0

# Cálculo inicial do Pb para o ponto de partida
TF_start = T_atual_R - 459.67
Pb_inicial_psia = PVT_OFICIAL.obter_pressao_bolha(RGL, dg, API, TF_start)

# Chamada inicial ao PVT para pegar densidades de partida
pvt_init = PVT_OFICIAL.main(P_atual_psia, T_atual_R, dg, T_atual_R, 10.73, 28.96, API, RGL, TF_start, 0)
# Indices baseados no retorno do PVT_OFICIAL: rhog (6), rhoo (15)
rho_g_init = pvt_init[6] * 16.018463 # lb/ft3 -> kg/m3
rho_o_init = pvt_init[15] * 16.018463

# Listas para armazenar resultados (para gráficos e Excel)
dados_simulacao = {
    "L_m": [0.0],
    "P_bar": [P_res_bar],
    "T_C": [T_res_C],
    "Holdup": [0.0],
    "Regime": ["N/A"],
    "Bo": [0.0],
    "Bg": [0.0],
    "Pb_bar": [Pb_inicial_psia / 14.5038],
    "Vsg": [0.0],
    "dp_dL": [0.0],
    "rho_o": [rho_o_init], # Densidade Oleo (kg/m3)
    "rho_g": [rho_g_init]  # Densidade Gas (kg/m3)
}

# Vazão Mássica Aproximada (para equação de temperatura)
q_m3s = Q_std_d / 86400.0
qm_kg_s = q_m3s * 850.0 # Estimativa inicial

print(f"Calculando {len(sections)} passos de simulação...")

for step in sections:
    dL = step['dL']
    theta = step['theta']
    T_amb_C = step['T_amb_C']
    TEC_linear = step['TEC']
    
    # --- A. CÁLCULO HIDRÁULICO (Pressão) ---
    pvt_inputs = (dg, T_atual_R, 10.73, 28.96, API, RGL, (T_atual_R - 459.67), 0)
    
    # Inicializa variável de perda caso dê erro
    dp_total_Pa_m = 0.0

    try:
        # Chama Beggs & Brill
        dp_total_Pa_m, Hl, regime, Bg, Bo, rho_mix_si = BB.calc_gradient(
            Q_std_d, D_m, Rugosidade, theta, P_atual_psia, T_atual_R, pvt_inputs, BSW
        )
        
        # RECALCULA PVT AQUI PARA PEGAR AS DENSIDADES SEPARADAS (Já que o BB retorna misturado)
        # Isso garante que temos rho_oleo e rho_gas exatos para o plot
        pvt_res = PVT_OFICIAL.main(P_atual_psia, T_atual_R, *pvt_inputs)
        rho_g_step = pvt_res[6] * 16.018463
        rho_o_step = pvt_res[15] * 16.018463
        
    except Exception as e:
        print(f"ERRO CRÍTICO em L={L_acumulado:.1f}m: {e}")
        break

    # --- CÁLCULO MANUAL DA VSG (Para o gráfico) ---
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
    
    # --- B. CÁLCULO TÉRMICO (Fórmula da Imagem + Cp das Meninas) ---
    T_amb_K = T_amb_C + 273.15
    T_old_K = (T_atual_R - 491.67) * 5.0/9.0 + 273.15
    theta_rad = math.radians(theta)

    # CÁLCULO DO Cp (IGUAL AO DELAS)
    Cp_kJ = ((2e-3) * T_amb_C - 1.429) * do + (2.67e-3) * T_amb_C + 3.049
    Cp_oleo = Cp_kJ * 1000.0 # Converte kJ para J

    # Termo A e Exponencial
    term_A = (qm_kg_s * 9.81 * math.sin(theta_rad)) / TEC_linear
    arg_exp = (TEC_linear * dL) / (qm_kg_s * Cp_oleo)
    fator_exp = math.exp(-arg_exp)

    # Fórmula Térmica
    T_new_K = (T_amb_K - term_A) - fator_exp * (T_amb_K - term_A - T_old_K)

    # Converte de volta para Rankine
    T_new_C = T_new_K - 273.15
    T_new_R = (T_new_C * 9.0/5.0) + 491.67
    
    # --- C. CÁLCULO DA PRESSÃO DE BOLHA (NOVO) ---
    TF_new = T_new_R - 459.67
    Pb_atual_psia = PVT_OFICIAL.obter_pressao_bolha(RGL, dg, API, TF_new)
    
    # --- D. ATUALIZAÇÃO E ARMAZENAMENTO ---
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
    dados_simulacao["dp_dL"].append(dp_total_Pa_m * 1e-5) # Armazena Perda em bar/m
    dados_simulacao["rho_o"].append(rho_o_step)
    dados_simulacao["rho_g"].append(rho_g_step)

# ============================================================================
# 4. RESULTADOS, RELATÓRIO E GRÁFICOS
# ============================================================================
df = pd.DataFrame(dados_simulacao)

print("\n" + "="*60)
print(f" RELATÓRIO DE INTERPRETAÇÃO (L = {L_acumulado:.1f} m)")
print("="*60)

if len(df) > 1:
    P_chegada = df["P_bar"].iloc[-1]
    T_chegada = df["T_C"].iloc[-1]
    Pb_chegada = df["Pb_bar"].iloc[-1]
    Regimes_todos = df["Regime"].unique()
    Hl_final = df["Holdup"].iloc[-1]

    # 1. ANÁLISE DE SURGÊNCIA
    print(f"\n1. ANÁLISE DE SURGÊNCIA (Capacidade de Produção):")
    print(f"   - Pressão Disponível na Plataforma: {P_chegada:.2f} bar")
    print(f"   - Pressão Mínima Requerida:         15.00 bar")
    
    if P_chegada >= 15.0:
        margem = P_chegada - 15.0
        print(f"   >>> STATUS: SUCESSO! O poço produz naturalmente.")
        print(f"   >>> Margem de segurança: {margem:.2f} bar disponíveis.")
    else:
        print(f"   >>> STATUS: FALHA. Pressão insuficiente.")
        print(f"   >>> Déficit: faltam {15.0 - P_chegada:.2f} bar para produzir.")

    # 2. ANÁLISE DE FASES (Líquido vs Gás)
    print(f"\n2. COMPORTAMENTO DE FASES (Ponto de Bolha):")
    print(f"   - Pressão de Bolha na Chegada: {Pb_chegada:.2f} bar")
    
    if P_chegada < Pb_chegada:
        try:
            idx_cruzamento = df[df["P_bar"] < df["Pb_bar"]].index[0]
            L_cruzamento = df["L_m"].iloc[idx_cruzamento]
            print(f"   >>> ALERTA: O fluido se tornou BIFÁSICO dentro da tubulação.")
            print(f"   >>> O gás começou a liberar em L = {L_cruzamento:.1f} m.")
            print(f"   >>> Holdup Final: {Hl_final:.2f} (Significa que {(1-Hl_final)*100:.1f}% do volume do tubo é gás).")
        except:
            pass
    else:
        print(f"   >>> SITUAÇÃO: O fluido permaneceu MONOFÁSICO (Líquido) em todo o trajeto.")
        print(f"   >>> A Pressão do Fluxo nunca caiu abaixo da Pressão de Bolha.")
        print(f"   >>> Holdup constante em 1.00 (Tubo cheio de óleo).")

    # 3. ANÁLISE DE PERDAS
    print(f"\n3. BALANÇO DE PERDAS (Energia e Calor):")
    delta_P = df["P_bar"].iloc[0] - P_chegada
    delta_T = df["T_C"].iloc[0] - T_chegada
    print(f"   - Perda Total de Pressão (Atrito + Gravidade): {delta_P:.2f} bar")
    print(f"   - Perda Total de Temperatura (Resfriamento):   {delta_T:.2f} °C")

    # 4. REGIMES DE ESCOAMENTO
    print(f"\n4. REGIMES DE ESCOAMENTO IDENTIFICADOS:")
    for r in Regimes_todos:
        print(f"   - {r}")

else:
    print("Simulação falhou no primeiro passo.")
print("="*60 + "\n")

try:
    df.to_excel("Relatorio_Grupo2.xlsx", index=False)
    print("Arquivo 'Relatorio_Grupo2.xlsx' gerado com sucesso.")
except Exception as e:
    print(f"Erro ao salvar Excel: {e}")

# --- PLOTAGEM ---
# 8 GRÁFICOS (4 LINHAS X 2 COLUNAS)
fig, axs = plt.subplots(4, 2, figsize=(14, 20))
fig.suptitle(f'Simulação de Escoamento - Grupo 2 (Q={Q_std_d} sm3/d)', fontsize=16)

# 1. Pressão
axs[0, 0].plot(df["L_m"], df["P_bar"], 'b-', lw=2, label="Pressão Fluxo")
axs[0, 0].plot(df["L_m"], df["Pb_bar"], 'k--', lw=1.5, label="Pressão Bolha")
axs[0, 0].axhline(15, color='r', linestyle='--', label='P Separador')
axs[0, 0].set_xlabel('Comprimento (m)')
axs[0, 0].set_ylabel('Pressão (bar)')
axs[0, 0].set_title('Perfil de Pressão')
axs[0, 0].grid(True)
axs[0, 0].legend()
axs[0, 0].set_xlim(0, L_total_sistema)

# 2. Temperatura
axs[0, 1].plot(df["L_m"], df["T_C"], 'r-', lw=2)
axs[0, 1].set_xlabel('Comprimento (m)')
axs[0, 1].set_ylabel('Temperatura (°C)')
axs[0, 1].set_title('Perfil de Temperatura')
axs[0, 1].grid(True)
axs[0, 1].set_xlim(0, L_total_sistema)

# 3. Holdup
axs[1, 0].plot(df["L_m"], df["Holdup"], 'g-', lw=2)
axs[1, 0].set_xlabel('Comprimento (m)')
axs[1, 0].set_ylabel('Holdup Líquido (-)')
axs[1, 0].set_title('Fração de Líquido')
axs[1, 0].grid(True)
axs[1, 0].set_ylim(0, 1.1)
axs[1, 0].set_xlim(0, L_total_sistema)

# 4. Fatores Volume (Bo e Bg)
axs[1, 1].plot(df["L_m"], df["Bo"] / 0.1589873, 'b-', label='Bo')
axs[1, 1].set_xlabel('Comprimento (m)')
axs[1, 1].set_ylabel('Bo', color='b')
axs[1, 1].tick_params(axis='y', labelcolor='b')
axs[1, 1].set_title('Fatores Volume (Bo e Bg)')
axs[1, 1].grid(True)
axs[1, 1].set_xlim(0, L_total_sistema)

ax4_twin = axs[1, 1].twinx()
ax4_twin.plot(df["L_m"], df["Bg"]/ 0.1589873, 'r--', label='Bg')
ax4_twin.set_ylabel('Bg', color='r')
ax4_twin.tick_params(axis='y', labelcolor='r')
lines_1, labels_1 = axs[1, 1].get_legend_handles_labels()
lines_2, labels_2 = ax4_twin.get_legend_handles_labels()
axs[1, 1].legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')

# 5. Perfil de Pressão de Bolha Isolado
axs[2, 0].plot(df["L_m"], df["Pb_bar"], 'purple', lw=2)
axs[2, 0].set_xlabel('Comprimento (m)')
axs[2, 0].set_ylabel('Pressão de Bolha (bar)')
axs[2, 0].set_title('Perfil de Pressão de Bolha')
axs[2, 0].grid(True)
axs[2, 0].set_xlim(0, L_total_sistema)

# 6. Perfil de Velocidade Superficial Gás (Vsg)
axs[2, 1].plot(df["L_m"], df["Vsg"], 'orange', lw=2)
axs[2, 1].set_xlabel('Comprimento (m)')
axs[2, 1].set_ylabel('Vsg (m/s)')
axs[2, 1].set_title('Velocidade Superficial do Gás')
axs[2, 1].grid(True)
axs[2, 1].set_xlim(0, L_total_sistema)

# 7. Perda de Carga por Comprimento
axs[3, 0].plot(df["L_m"], df["dp_dL"], 'brown', lw=2)
axs[3, 0].set_xlabel('Comprimento (m)')
axs[3, 0].set_ylabel('Perda de Carga (bar/m)')
axs[3, 0].set_title('Gradiente de Pressão (dP/dL)')
axs[3, 0].grid(True)
axs[3, 0].set_xlim(0, L_total_sistema)

# 8. NOVO GRÁFICO: Densidades
axs[3, 1].plot(df["L_m"], df["rho_o"], 'b-', label='Óleo')
axs[3, 1].set_xlabel('Comprimento (m)')
axs[3, 1].set_ylabel('Densidade Óleo (kg/m³)', color='b')
axs[3, 1].tick_params(axis='y', labelcolor='b')
axs[3, 1].set_title('Densidades dos Fluidos')
axs[3, 1].grid(True)
axs[3, 1].set_xlim(0, L_total_sistema)

ax8_twin = axs[3, 1].twinx()
ax8_twin.plot(df["L_m"], df["rho_g"], 'r--', label='Gás')
ax8_twin.set_ylabel('Densidade Gás (kg/m³)', color='r')
ax8_twin.tick_params(axis='y', labelcolor='r')

lines_8a, labels_8a = axs[3, 1].get_legend_handles_labels()
lines_8b, labels_8b = ax8_twin.get_legend_handles_labels()
axs[3, 1].legend(lines_8a + lines_8b, labels_8a + labels_8b, loc='best')

plt.tight_layout()
plt.show()
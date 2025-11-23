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
dg = 0.75               # Gravidade específica gás (CORRIGIDO DE 'g' PARA 'dg')

# --- Reservatório (Condição de Contorno Inicial - Fundo) ---
P_res_bar = 550.0       # bar
T_res_C = 80.0          # °C

# --- Geometria do Duto ---
D_pol = 3           # Diâmetro (suposição padrão, ajuste se necessário)
D_m = D_pol * 0.0254    # m
Rugosidade = 4.5e-5     # m (Aço carbono típico)

# --- Propriedades Térmicas ---
Cp_oleo = 2200.0        # J/kg.K (Valor estimado para óleo cru)
TEC_poco = 2.0          # W/m.K (Linear)
TEC_marinho = 1.0       # W/m.K (Linear)

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
# Vertical (90 graus). O fluido SUBINDO.
# Z vai de 2100m para 950m. Delta Z = 1150m.
L_poco = 2100 - 950
n_steps1 = int(L_poco / dL_step)

for i in range(n_steps1):
    # Cota atual (começa em 2100 e diminui subindo)
    z_atual = 2100 - (i * dL_step)
    
    # [cite_start]Temperatura Externa: Varia linearmente de 80°C (fundo) a 4°C (ANM) [cite: 12]
    # Fração percorrida (0 no fundo, 1 na ANM)
    frac = i / n_steps1 
    T_amb = 80.0 + (4.0 - 80.0) * frac
    
    sections.append({
        "nome": "Poço",
        "theta": 90.0,       # Vertical subindo
        "T_amb_C": T_amb,
        "TEC": TEC_poco,
        "dL": dL_step
    })

# --- TRECHO 2: FLOWLINE (ANM -> Manifold) ---
# Inclinado (Leito Marinho). ANM (950m) -> Manifold (850m).
# Delta Z = 100m. [cite_start]Ângulo = 37 graus[cite: 12].
# Comprimento da tubulação = Delta Z / sin(37)
L_flowline = 100.0 / math.sin(math.radians(37.0))
n_steps2 = int(L_flowline / dL_step)

for i in range(n_steps2):
    # [cite_start]Temperatura Externa: Constante 4°C no leito marinho [cite: 13]
    sections.append({
        "nome": "Flowline",
        "theta": 37.0,       # Inclinado subindo
        "T_amb_C": T_fundo_mar_C,
        "TEC": TEC_marinho,
        "dL": dL_step
    })

# --- TRECHO 3: RISER (Manifold -> Plataforma) ---
# Vertical. Manifold (850m) -> Superfície (0m).
L_riser = 850.0
n_steps3 = int(L_riser / dL_step)

for i in range(n_steps3):
    z_atual = 850 - (i * dL_step)
    
    # [cite_start]Temperatura Externa: Varia linearmente de 4°C (Manifold) a 17°C (Superfície) [cite: 13]
    frac = i / n_steps3
    T_amb = 4.0 + (17.0 - 4.0) * frac
    
    sections.append({
        "nome": "Riser",
        "theta": 90.0,       # Vertical subindo
        "T_amb_C": T_amb,
        "TEC": TEC_marinho,
        "dL": dL_step
    })
# ============================================================================
# 3. LOOP DE SIMULAÇÃO (MARCHA PxT)
# ============================================================================

# Conversões Iniciais
P_atual_psia = P_res_bar * 14.5038
T_atual_R = (T_res_C * 9.0/5.0) + 491.67
L_acumulado = 0.0

# Listas para armazenar resultados (para gráficos e Excel)
dados_simulacao = {
    "L_m": [0.0],
    "P_bar": [P_res_bar],
    "T_C": [T_res_C],
    "Holdup": [0.0],
    "Regime": ["N/A"],
    "Bo": [0.0],
    "Bg": [0.0]
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
    # CORREÇÃO AQUI: Adicionado T_atual_R na lista de argumentos (antes do R=10.73)
    # A ordem correta do PVT é: dg, T, R, Mar, API, RGL, TF, S
    pvt_inputs = (dg, T_atual_R, 10.73, 28.96, API, RGL, (T_atual_R - 459.67), 0)
    
    try:
        # Chama a função do arquivo BEGGS_BRILL
        dp_total_Pa_m, Hl, regime, Bg, Bo, rho_mix_si = BB.calc_gradient(
            Q_std_d, D_m, Rugosidade, theta, P_atual_psia, T_atual_R, pvt_inputs, BSW
        )
    except Exception as e:
        print(f"ERRO CRÍTICO em L={L_acumulado:.1f}m: {e}")
        # Se der erro, para o loop para não repetir 200 vezes
        break

    # Atualiza Pressão: P_new = P_old - (Gradiente * dL)
    dp_psi_m = dp_total_Pa_m * 1.45038e-4
    P_new_psia = P_atual_psia - (dp_psi_m * dL)
    
    # --- B. CÁLCULO TÉRMICO (Fórmula da Imagem) ---
    T_amb_K = T_amb_C + 273.15
    T_old_K = (T_atual_R - 491.67) * 5.0/9.0 + 273.15
    theta_rad = math.radians(theta)

    # Termo A (Fonte de Calor Potencial - qm * g * sin(theta) / TEC)
    term_A = (qm_kg_s * 9.81 * math.sin(theta_rad)) / TEC_linear

    # Termo Exponencial (Relaxamento Térmico - exp(-TEC * dL / (qm * Cp)))
    arg_exp = (TEC_linear * dL) / (qm_kg_s * Cp_oleo)
    fator_exp = math.exp(-arg_exp)

    # Aplicação da Fórmula: T_novo = (T_amb - TermoA) - Exp * (T_amb - TermoA - T_old)
    # Esta linha implementa a fórmula da sua imagem:
    T_new_K = (T_amb_K - term_A) - fator_exp * (T_amb_K - term_A - T_old_K)

    # Converte de volta para Rankine
    T_new_C = T_new_K - 273.15
    T_new_R = (T_new_C * 9.0/5.0) + 491.67
    
    # --- C. ATUALIZAÇÃO E ARMAZENAMENTO ---
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

# ============================================================================
# 4. RESULTADOS E GRÁFICOS
# ============================================================================
df = pd.DataFrame(dados_simulacao)

print("-" * 40)
print(f"RESULTADO FINAL (L = {L_acumulado:.1f} m)")
if len(df) > 1:
    P_chegada = df["P_bar"].iloc[-1]
    T_chegada = df["T_C"].iloc[-1]
    print(f"Pressão de Chegada: {P_chegada:.2f} bar")
    print(f"Temperatura de Chegada: {T_chegada:.2f} °C")

    if P_chegada >= 15.0:
        print("STATUS: ESCOAMENTO GARANTIDO (P > 15 bar)")
    else:
        print("STATUS: FALHA DE ESCOAMENTO (P < 15 bar)")
else:
    print("Simulação falhou no primeiro passo.")
print("-" * 40)

# Tenta Exportar Excel (Protegido contra erro de instalação)
try:
    df.to_excel("Relatorio_Grupo2.xlsx", index=False)
    print("Arquivo 'Relatorio_Grupo2.xlsx' gerado com sucesso.")
except ImportError:
    print("AVISO: Biblioteca 'openpyxl' não instalada. O Excel não foi gerado, mas os gráficos aparecerão.")
except Exception as e:
    print(f"Erro ao salvar Excel: {e}")

# --- PLOTAGEM ---
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(f'Simulação de Escoamento - Grupo 2 (Q={Q_std_d} sm3/d)', fontsize=16)

# 1. Pressão
ax1.plot(df["L_m"], df["P_bar"], 'b-', lw=2)
ax1.axhline(15, color='r', linestyle='--', label='P Separador (15 bar)')
ax1.set_xlabel('Comprimento (m)')
ax1.set_ylabel('Pressão (bar)')
ax1.set_title('Perfil de Pressão')
ax1.grid(True)
ax1.legend()

# 2. Temperatura
ax2.plot(df["L_m"], df["T_C"], 'r-', lw=2)
ax2.set_xlabel('Comprimento (m)')
ax2.set_ylabel('Temperatura (°C)')
ax2.set_title('Perfil de Temperatura')
ax2.grid(True)

# 3. Holdup
ax3.plot(df["L_m"], df["Holdup"], 'g-', lw=2)
ax3.set_xlabel('Comprimento (m)')
ax3.set_ylabel('Holdup Líquido (-)')
ax3.set_title('Fração de Líquido')
ax3.grid(True)
ax3.set_ylim(0, 1.1)

# 4. Fatores Volume (Bo e Bg)
ax4.plot(df["L_m"], df["Bo"] / 0.1589873, 'b-', label='Bo')
ax4.set_xlabel('Comprimento (m)')
ax4.set_ylabel('Bo', color='b')
ax4.tick_params(axis='y', labelcolor='b')
ax4.set_title('Fatores Volume (Bo e Bg)')
ax4.grid(True)

ax4_twin = ax4.twinx()
ax4_twin.plot(df["L_m"], df["Bg"]/ 0.1589873, 'r--', label='Bg')
ax4_twin.set_ylabel('Bg', color='r')
ax4_twin.tick_params(axis='y', labelcolor='r')

lines_1, labels_1 = ax4.get_legend_handles_labels()
lines_2, labels_2 = ax4_twin.get_legend_handles_labels()
ax4.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')

plt.tight_layout()
plt.show()

# Adicione isso no final do seu código MAIN_SIMULATOR.py para tirar a prova real
print(f"DEBUG: Pressão de Chegada: {P_chegada:.2f} bar")
# Estimar Pb rapidinho com os dados de entrada
Pb_estimado = 18.2 * ((RGL / dg) ** 0.83) * (10 ** (0.00091 * T_chegada - 0.0125 * API)) - 1.4
print(f"DEBUG: Ponto de Bolha Estimado na Chegada: {Pb_estimado * 0.0689:.2f} bar")

if P_chegada > (Pb_estimado * 0.0689):
    print("CONCLUSÃO: O fluido chega pressurizado acima do Ponto de Bolha.")
    print("          Por isso o Holdup é constante (quase tudo líquido) e a reta é linear.")
import numpy as np
import math

# =======================
# 1. PARÂMETROS GERAIS
# =======================
g = 9.81

# Comprimentos
L_poco = 1500        # m (trecho inclinado)
L_horiz = 6000       # m
L_marinho = 30       # m

# Temperaturas do ambiente (°C → Rankine)
def to_R(C):
    return (C * 9/5) + 491.67

T1 = to_R(80)   # entrada do poço
T2 = to_R(4)    # início do manifold
T3 = to_R(15)   # superfície

# Propriedades do fluido / térmicas
TEC_poco = 2       # coeficiente de troca térmica
TEC_marinho = 1
do = 0.8           # densidade relativa (exemplo)
Q_l = 10000/86400  # vazão em m³/s
q_m = Q_l * do     # vazão mássica aproximada
Cp = 2500          # calor específico (J/kg.K) — coloque o seu valor

# Ângulo do trecho inclinado
theta = math.radians(60)

# =======================
# 2. MALHA DE CÁLCULO
# =======================
n = 500
x_inc = np.linspace(0, L_poco, n)
x_hor = np.linspace(0, L_horiz, n)
x_mar = np.linspace(0, L_marinho, n)

# =======================
# 3. FUNÇÃO (modelo das meninas)
# =======================

def calcular_trecho_inclinado(T_in, T_out, L, theta, q_m, Cp, TEC):
    """Replica a fórmula exatamente como no material."""
    n = 500
    x = np.linspace(0, L, n)

    dT = (T_in - T_out) / n     # gradiente linear ambiente
    T_amb = T_in - dT * np.arange(n)

    T_calc = []
    T_old = T_in

    for i in range(n):

        termo_grav = (q_m * g * math.sin(theta)) / TEC

        expoente = math.exp(-(TEC / (q_m * Cp)) * L)

        T_i = (
            T_amb[i]
            - termo_grav
            - expoente * (T_amb[i] - termo_grav - T_old)
        )

        T_calc.append(T_i)
        T_old = T_i

    return np.array(T_calc)


def calcular_trecho_vertical(T_in, T_out, L, q_m, Cp, TEC):
    """Mesma lógica do inclinado, mas sem seno(theta)."""
    n = 500
    x = np.linspace(0, L, n)

    dT = (T_in - T_out) / n
    T_amb = T_in - dT * np.arange(n)

    T_calc = []
    T_old = T_in

    for i in range(n):

        termo_grav = (q_m * g) / TEC   # sem sin(theta)

        expoente = math.exp(-(TEC / (q_m * Cp)) * L)

        T_i = (
            T_amb[i]
            - termo_grav
            - expoente * (T_amb[i] - termo_grav - T_old)
        )

        T_calc.append(T_i)
        T_old = T_i

    return np.array(T_calc)


# =======================
# 4. CALCULANDO TODOS OS TRECHOS
# =======================

# Trecho inclinado
T_inc = calcular_trecho_inclinado(
    T1, T2, L_poco, theta,
    q_m, Cp, TEC_poco
)

# Trecho horizontal (temperatura ambiente constante)
T_horiz = np.full_like(x_hor, T2)

# Trecho vertical marinho
T_marinho = calcular_trecho_vertical(
    T2, T3, L_marinho,
    q_m, Cp, TEC_marinho
)

# Temperatura total
T_total = np.concatenate([T_inc, T_horiz, T_marinho])
x_total = np.concatenate([x_inc, x_hor + L_poco, x_mar + L_poco + L_horiz])


# =======================
# 5. EXEMPLO DE PLOT
# =======================
import matplotlib.pyplot as plt

plt.plot(x_total, T_total)
plt.xlabel("Comprimento (m)")
plt.ylabel("Temperatura (°R)")
plt.title("Perfil Térmico (Modelo Igual ao das Meninas)")
plt.grid()
plt.show()

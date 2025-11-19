import numpy as np

# --- DADOS INICIAIS ---
dg = 0.84  # Gravidade relativa do gás
do = 0.86  # Gravidade relativa do óleo
Pb = 5000  # psia - Pressão de bolha (este valor será redefinido pelo cálculo da função)
P = 3626   # psia - Pressão atual
Tf = 122   # ºF - Temperatura em Fº
Tsep = 80  # ºF - Temperatura do separador
Psep = 100 # psia - Pressão do separador
Mar = 28.96  # lb/lbmol - Massa molar do ar
R = 10.73  # psi*ft³/lbmol.ºR - Constante universal do gás
rhoob =float('nan')   # no caso da nossa aplicação nao possui
Bob = float('nan')    # no caso da nossa aplicação nao possui

# --- CÁLCULOS INICIAIS ---

# Mudar de ºF para ºR
def calculate_TR(Tf):
    T = (Tf + 460)  # temperatura em Rº
    return T

T = calculate_TR(Tf)

# Grau API do óleo
def calculate_API(do):
    API = ((141.5 / do) - 131.5)
    return API

API = calculate_API(do)

# Cálculo da gravidade específica (recalculado com base no API)
def calculate_do(API):
    do_calc = 141.5 / (API + 131.5)
    return do_calc

do = calculate_do(API)

# Massa molecular do gás
def calculate_Mg(Mar, dg):
    Mg = Mar * dg  # lb/lbmol
    return Mg

Mg = calculate_Mg(Mar, dg)
#razão de solubilidade gás-oleo
def calc_rs(P, API, Tf, dg):
    Rs_calc = ((((P / 18.2) + 1.4) * ((10**(0.0125 * API)) / (10**(0.00091 * Tf))))**(1 / 0.83)) * dg
    return Rs_calc



# Definir pressão de bolha pela correlação de Standing
def pressao_bolha(Rs, API, Tf, dg):
    apb = (0.00091 * Tf) - (0.0125 * API)
    Pb_calc = 18.2 * (((Rs / dg)**0.83) * (10**apb) - 1.4)
    return Pb_calc


# --- CÁLCULO DAS PROPRIEDADES PSEUDOCRÍTICAS (Ppc, Tpc) ---
def P_pseudocritica1(dg):
    Ppc = 677 + 15 * dg - 37.5 * dg**2  # Ppseudocritica em psia
    return Ppc

def T_pseudocritica1(dg):
    Tpc = 168 + 325 * dg - 12.5 * dg**2  # Tpseudocritica em ºR
    return Tpc

def P_pseudocritica2(dg):
    Ppc = 706 - 51.7 * dg - 11.1 * dg**2  # Ppseudocritica em psia
    return Ppc

def T_pseudocritica2(dg):
    Tpc = 187 + 330 * dg - 71.5 * dg**2  # Tpseudocritica em ºR
    return Tpc

if dg < 0.75:
    Ppc = P_pseudocritica1(dg)
    Tpc = T_pseudocritica1(dg)
else:
    Ppc = P_pseudocritica2(dg)
    Tpc = T_pseudocritica2(dg)

# --- Cálculo da Pressão e Temperatura Pseudoreduzidas ---
def P_pseudoreduzida(P, Ppc):
    Ppr = P / Ppc  # Ppseudoreduzida em psi
    return Ppr

Ppr = P_pseudoreduzida(P, Ppc)

def T_pseudoreduzida(T, Tpc):
    Tpr = T / Tpc  # Tpseudoreduzida em ºR
    return Tpr

Tpr = T_pseudoreduzida(T, Tpc)

# --- Usando correlação de Papay para encontrar fator Z ---
def calculate_Z_papay(Ppr, Tpr):
    Z = 1 - ((3.53 * Ppr) / (10**(0.9813 * Tpr))) + ((0.274 * (Ppr**2)) / (10**(0.8157 * Tpr)))
    return Z

Z = calculate_Z_papay(Ppr, Tpr)

# --- Fator de Volume de Formação do Gás ---)
def fator_Vformacao_gas(Tf, P, Z):
    Bg = (14.7 / 60) * (Z * (Tf  / P))
    return Bg

Bg = fator_Vformacao_gas(Tf, P, Z)

# --- CÁLCULO DA MASSA ESPECÍFICA DO GÁS ---
def calculate_rho_g(P, Mg, Z, R, T):
    rho_g = (P * Mg) / (Z * R * T)
    return rho_g

rho_g = calculate_rho_g(P, Mg, Z, R, T)

# --- VISCOSIDADE DA FASE GÁS (Correlação de Lee) ---
def calculate_u_g_lee(T, Mg, rho_g):
    xv = 3.448 + (986.4 / T) + (0.01009 * Mg)  # parâmetro
    yv = 2.4 - (0.2 * xv)  # parâmetro
    kv = ((9.379 + (0.0160 * Mg)) * (T**1.5)) / (209.2 + (19.26 * Mg) + T)  # parâmetro
    ug = (10**(-4)) * kv * np.exp(xv * ((rho_g / 62.4)**yv))
    return ug

ug = calculate_u_g_lee(T, Mg, rho_g)  # viscosidade do gás em cP

# --- COMPRESSIBILIDADE DO GÁS ---
def compressibiliade_gas(Ppr, Tpr, Z, Ppc):
    dzdPpr = -(3.53 / (10**(0.9813 * Tpr))) + 2 * (0.274 / (10**(0.8157 * Tpr))) * Ppr
    Cpr = (1 / Ppr) - (1 / Z) * dzdPpr
    Cg = Cpr / Ppc
    return Cg

Cg = compressibiliade_gas(Ppr, Tpr, Z, Ppc)  # compress. do gás

# Calculo de Rs e Pb antes das condições
Rs = calc_rs(P, API, Tf, dg)
Pb = pressao_bolha(Rs, API, Tf, dg)

# Inicialização de variáveis para garantir que existam fora dos blocos if/else
Co = 0
Bo = 0
rho_o = 0
uobb = 0
uob = 0

if P > Pb:
    # --- Compressibilidade do Óleo (P > Pb) ---
    def compressibilidade_oleo_subsaturado(rhoob, P, Pb):
        termo1 = rhoob + 0.004347 * (P - Pb) - 79.1
        termo2 = 0.0007141 * (P - Pb) - 12.938
        Co_calc = (10**(-6)) * np.exp(termo1 / termo2)
        return Co_calc
    Co = compressibilidade_oleo_subsaturado(rhoob, P, Pb)

    # --- Cálculo do Fator Volume-Formação do Óleo (Correlação de Standing para P > Pb) ---
    def calculate_Bo_subsaturado(Bob, Co, Pb, P):
        Bo_calc = Bob * np.exp(-Co * (Pb - P))
        return Bo_calc
    Bo = calculate_Bo_subsaturado(Bob, Co, Pb, P)

    # --- Cálculo da Massa Específica do Óleo (P > Pb) ---
    def calculate_rho_o_subsaturado(rhoob, Co, P, Pb):
        rho_o_calc = rhoob * np.exp(Co * (P - Pb))
        return rho_o_calc
    rho_o = calculate_rho_o_subsaturado(rhoob, Co, P, Pb)

    # --- Viscosidade do Óleo Subsaturado por Beggs ---
    def calculate_u_o_beggs_subsaturado(API, Tf, Rs, P, Pb):
        Aa = 10**(3.0324 - 0.02023 * API)  # parâmetro
        uodd = 10**(Aa * Tf**(-1.163)) - 1  # viscosidade do óleo morto
        aa = 10.715 * (Rs + 100)**(-0.515)  # parâmetro
        bb = 5.44 * (Rs + 150)**(-0.338)  # parâmetro
        uobb_calc = aa * (uodd**bb)
        m = 2.6 * (P**(1.187)) * np.exp(-11.513 - (8.98 * 10**(-5)) * P)
        uo = uobb_calc * ((P / Pb)**m)
        return uo
    uo = calculate_u_o_beggs_subsaturado(API, Tf, Rs, P, Pb)
    uobb = uo # Usando uobb para representar a viscosidade do óleo para P > Pb neste caso
    uob = uo # Usando uob para representar a viscosidade do óleo para P > Pb neste caso

else: # P <= Pb
    # --- Cálculo do Fator Volume-Formação do Óleo (Correlação de Standing para P <= Pb) ---
    def calc_bo(Rs, Tf, do, dg):
        Bo_calc = 0.9759 + 0.00012 * (Rs * ((dg / do)**0.5) + 1.25 * Tf)**1.2
        return Bo_calc
    Bo = calc_bo(Rs, Tf, do, dg)

    # --- Cálculo da Massa Específica do Óleo (P <= Pb) ---
    def calculate_rho_o_saturado(Rs, do, dg, Bo):
        rho_o_calc = (62.4 * do + 0.0136 * Rs * dg) / Bo
        return rho_o_calc
    rho_o = calculate_rho_o_saturado(Rs, do, dg, Bo)

    # --- Compressibilidade do Óleo por Standing (P <= Pb) ---
    def compressibilidade_oleo_saturado(P, API, dg, do, Tf):
        delta_P = 50  # Incremento pequeno de pressão (psia)
        # Derivadas numéricas centradas
        Rs_up = calc_rs(P + delta_P, API, Tf, dg)
        Rs_down = calc_rs(P - delta_P, API, Tf, dg)
        dRs_dP = (Rs_up - Rs_down) / (2 * delta_P)

        Bo_up = calc_bo(Rs_up, Tf, do, dg)
        Bo_down = calc_bo(Rs_down, Tf, do, dg)
        dBo_dP = (Bo_up - Bo_down) / (2 * delta_P)
        Co_calc = ((-1 / Bo) * dBo_dP) + ((Bg / Bo) * dRs_dP)
        return Co_calc
    Co = compressibilidade_oleo_saturado(P, API, dg, do, Tf)  # compress. do oleo

    # --- Viscosidade do Óleo pela Correlação de Beggs (Óleo Saturado) ---
    def calculate_u_o_beggs_saturado(API, Tf, Rs):
        Aa = 10**(3.0324 - 0.02023 * API)  # parâmetro
        uodd = 10**(Aa * Tf**(-1.163)) - 1  # viscosidade do óleo morto
        aa = 10.715 * (Rs + 100)**(-0.515)  # parâmetro
        bb = 5.44 * (Rs + 150)**(-0.338)  # parâmetro
        uobb_calc = aa * (uodd**bb)
        return uobb_calc
    uobb = calculate_u_o_beggs_saturado(API, Tf, Rs)  # viscosidade oleo saturado em cP

    # --- Viscosidade do Óleo pela Correlação de Standing (Óleo Saturado) ---
    def calculate_u_o_standing_saturado(API, T, Rs):
        A = 10**(0.43 + (8.33 / API))  # parametro
        uod = 0.32 + ((1.8 * (10**7)) / (API**4.53)) * ((360 / (T - 260))**A)  # viscosidade oleo morto
        a = 10**(((-7.4 * 10**(-4)) * Rs) + (2.2 * (10**(-7)) * (Rs**2)))
        b1 = (0.68 / (10**(Rs * (8.62 * 10**(-5)))))
        b2 = (0.25 / (10**(Rs * (1.1 * 10**(-3)))))
        b3 = (0.062 / (10**(Rs * (3.74 * 10**(-3)))))
        b = b1 + b2 + b3
        uob_calc = a * (uod**b)
        return uob_calc
    uob = calculate_u_o_standing_saturado(API, T, Rs)  # viscosidade oleo saturado em cP

    # --- RESULTADOS ---
print(f"Rs = {Rs:.2f} scf/STB")
print(f"Bo = {Bo:.4f} bbl/STB")
print(f"Viscosidade do óleo (Standing) = {uob:.2f} cP")
print(f"Viscosidade do óleo (Beggs & Robinson) = {uobb:.2f} cP")
print(f"Massa específica do óleo = {rho_o:.2f} lb/ft³")
print(f"Massa específica do gás = {rho_g:.2f} lb/ft³")
print(f"Viscosidade do gás (Lee) = {ug:.5f} cP")
print(f"Compressibilidade do gás (Cg) = {Cg:.6e} psi⁻¹")
print(f"Compressibilidade do óleo com gás em solução (Co) = {Co:.6e} psi⁻¹")
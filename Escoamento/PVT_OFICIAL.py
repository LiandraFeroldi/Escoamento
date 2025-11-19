
import math
import numpy as np
   

def propriedades_pseudocriticas(dg):
    if dg < 0.75:
        Ppc = 677 - 15.0 * dg - 37.5 * dg**2
        Tpc = 168 + 325 * dg - 12.5 * dg**2
    else:
        Ppc = 706 - 51.7 * dg - 11.1 * dg**2
        Tpc = 187 + 330 * dg - 71.5 * dg**2
    return Ppc, Tpc

def propriedades_pseudoreduzidas(P, T, Ppc, Tpc):
    Ppr = P / Ppc
    Tpr = T / Tpc
    return Ppr, Tpr

def fator_compressibilidade_papay(Ppr, Tpr):
    # Fator de compressibilidade (Papay, 1985)
    z = 1 - ((3.53 * Ppr) / (10 ** (0.9813 * Tpr))) + ((0.274 * Ppr**2) / (10 ** (0.8157 * Tpr)))
    return z

def compressibilidade_gas_INSITU(P, Ppc, z, Tk, dg, Ppr, Tpr):
    Mar = 28.96
    Mg = Mar * dg
    R = 10.7316
    rhog = (P * Mg) / (z * R * Tk)

    dZ_dPpr = -3.53 / (10 ** (0.9813 * Tpr)) + (2 * 0.274 * Ppr) / (10 ** (0.8157 * Tpr))
    dZ_dP = dZ_dPpr / Ppc

    Cg = (1 / P) - (1 / z) * dZ_dP
    return Cg, rhog, Mg
#--------------------------------------------------------------------

# Dados do gás 
def dados_gas(Mar, dg, P, z, R, TR):
    Mg = Mar * dg
    # T deve estar em Rankine se R=10.7316 (ft·lb/(lb·°R)) e P em psi
    rhog = (P * Mg) / (z * R * TR)    # lb/ft3 (aprox), verifique unidades!
    return Mg, rhog

def viscosidade_gas_lee(Mg, T, rhog):
    kv = ((9.379 + 0.0160 * Mg) * T ** 1.5) / (209.2 + 19.26 * Mg + T)
    xv = 3.448 + (986.4 / T) + 0.01009 * Mg
    yv = 2.4 - 0.2 * xv
    ug = 1e-4 * kv * math.exp(xv * (rhog / 62.4) ** yv)
    return ug

def fator_formacao_gas(z, TR, P):
    # T em Rankine, P em psi — verifique formula original
    Bg = (14.7 / 60) * z * (TR / P)
    return Bg

#-------------------------------------------------------------------------------------------------
# PROPRIEDADES DO ÓLEO
def dados_oleo(Api):
    do = (141.5 / (131.5 + Api))
    rhoo=do*1000 # N SEI SE A GENTE PODE USAR ASSIM PRA ACHAR O RHO DO ÓLEO
    return do, rhoo

def obter_pressao_bolha(RGL, dg, Api, TF):
    #Correlação de Standing
    a = 0.00091 * TF - 0.0125 * Api
    Pb = 18.8 * (((RGL / dg) ** 0.83) * 10 ** a - 1.4)
    return Pb

def razao_solubilidade_gas_oleo(P, dg, Api, TF, Pb):
    if P >= Pb:
        Rs = dg * (((Pb / 18.2) + 1.4) * 10 ** (0.0125 * Api - 0.00091 * TF)) ** (1 / 0.83)
    else:
        Rs = dg * (((P / 18.2) + 1.4) * 10 ** (0.0125 * Api - 0.00091 * TF)) ** (1 / 0.83)
    return Rs

def fator_volume_formacao_oleo(Rs, dg, do, TF, P, Pb, RGL, Co):
    Bob = 0.9759 + 0.00012 * ((RGL * ((dg / do) ** 0.5) + 1.25 * TF) ** 1.2)
    if P > Pb:
        Bo = Bob * math.exp(-Co * (P - Pb))
    else:
        Bo = 0.9759 + 0.00012 * ((Rs * ((dg / do) ** 0.5) + 1.25 * TF) ** 1.2)
    return Bo, Bob

def massa_especifica_oleo_INSITU(Rs, dg, do, Bo, P, Pb, RGL, Co):
    if P > Pb:
        rhoob = (62.4 * do + 0.0136 * RGL * dg) / Bo
        rhoo = rhoob * math.exp(Co * (P - Pb))
        return rhoo, rhoob
    else:
        rhoo = (62.4 * do + 0.0136 * Rs * dg) / Bo
        return rhoo, None

def viscosidade_oleo_morto(Api, TF):
    #Viscosidade do Óleo Morto Correlação de Beggs e Robinson (1975)
    A = 10 ** (3.0324 - 0.02023 * Api)
    uom = 10 ** (A * TF ** -1.163) - 1
    return uom

def viscosidade_oleo_saturado(P, Pb, Rs, uom):
    #Viscosidade do Óleo saturado Correlação de Beggs e Robinson (1975)
    if P <= Pb:
        c = 10.715 * (Rs + 100) ** (-0.515)
        b = 5.44 * (Rs + 150) ** (-0.338)
        uos = c * uom ** b
    else:
        # para P > Pb use aproximação (original tinha m dependente de P)
        m = 2.68 * (P ** 1.187) * math.exp(-11.513 - (8.98e-5) * P)
        uos = uom * (P / Pb) ** m
    return uos

def compressibilidade_oleo(P, TF, dg, do, Rs, Bo, Api, Bg, Pb, rhoob):
    if P >= Pb:
        Co = 10**6 * np.exp((rhoob + 0.004347 * (P - Pb) - 79.1) / (0.0007141 * (P - Pb) - 12.938))
    else:
        dRs_dP = dg * (1 / 0.83) * (((P / 18.2) + 1.4) * (10 ** (0.0125 * Api - 0.00091 * TF))) ** (1 / 0.83 - 1) * (1 / 18.2) * (10 ** (0.0125 * Api - 0.00091 * TF))
        dBo_dP = 0.00012 * 1.2 * (Rs * ((dg / do) ** 0.5) + 1.25 * TF) ** 0.2 * ((dg / do) ** 0.5) * dRs_dP
        Co = -(1 / Bo) * dBo_dP + (Bg / Bo) * dRs_dP
    return Co


# --- Água ---
def massa_especifica_agua(S):
    rhow = 62.368 + 0.438603 * S + 1.60074e-3 * (S**2)
    return rhow

def Razao_de_Solubilidade_agua(P, TF):
    A0 = 8.15839
    A1 = -6.12265e-2
    A2 = 1.91663e-4
    A3 = -2.1654e-7
    B0 = 1.01021e-2
    B1 = -7.44241e-5
    B2 = 3.05553e-7
    B3 = -2.94883e-10
    C0 = -9.02505
    C1 = 0.130237
    C2 = -8.53425e-4
    C3 = 2.34122e-6
    C4 = -2.37049e-9

    A = A0 + A1*TF + A2*TF**2 + A3*TF**3
    B = B0 + B1*TF + B2*TF**2 + B3*TF**3
    C_interno = C0 + C1*TF + C2*TF**2 + C3*TF**3 + C4*TF**4
    C = C_interno * 1e-7
    Rsw = A + B*P + C*P**2
    return Rsw

def Volume_formacao_agua(P, TF):
    delta_VwT = -1.0001e-2 + 1.33391e-4 * TF+ 5.50654e-7 * (TF ** 2)
    termo1 = -1.95301e-9 * P * TF
    termo2 = -1.72834e-13 * (P**2) * TF
    termo3 = -3.58922e-7 * P
    termo4 = -2.25341e-10 * (P**2)
    delta_Vwp = termo1 + termo2 + termo3 + termo4
    Bw = (1 + delta_VwT) * (1 + delta_Vwp)
    return Bw

def viscosidade_agua(P, TF, S):
    A0 = 109.527
    A1 = -8.40564
    A2 = 0.313314
    A3 = 8.72213e-3
    B0 = -1.12166
    B1 = 2.63951e-2
    B2 = -6.7946e-4
    B3 = -5.47119e-5
    B4 = -1.55586e-6

    A = A0 + A1*S + A2*(S**2) + A3*(S**3)
    B = B0 + B1*S + B2*(S**2) + B3*(S**3) + B4*(S**4)
    mu_w1 = A * (TF**B)
    pressao_factor = 0.9994 + 4.0295e-5*P + 3.1062e-9*(P**2)
    uw = mu_w1 * pressao_factor
    return uw


#====================================================================================

def propriedades_gas(P, TR, dg,T,R):
    Ppc, Tpc = propriedades_pseudocriticas(dg)
    Ppr, Tpr = propriedades_pseudoreduzidas(P, T, Ppc, Tpc)
    z = fator_compressibilidade_papay(Ppr, Tpr)
    Mg=dados_gas(Mar,dg,z,R,T)
    ug = viscosidade_gas_lee(Mg, T, rhog)
    Bg = fator_formacao_gas(z, T, P)
    return Ppc, Tpc, Ppr, Tpr, z, rhog, Mg, ug, Bg

def propriedades_oleo(P, TF, dg, do, Pb, Rs, z, Bob, rhoob):
    Bg = fator_formacao_gas(z, TF, P)
    Bo = fator_volume_formacao_oleo(Rs, dg, do, TF, P, Pb, Bob)
    rhoo = massa_especifica_oleo_INSITU(Rs, dg, do, Bo, P, Pb, rhoob)
    uos = viscosidade_oleo_saturado(P, Pb, Rs)

    return Api, Rs, Pb, Bo, rhoo, uos, Bg

def propriedades_agua(S,P,TF):
    rhow=massa_especifica_agua(S)
    Rsw=Razao_de_Solubilidade_agua(P,TF)
    Bw=Volume_formacao_agua(P,TF)
    uw=viscosidade_agua(P,TF,S)

    return rhow, Rsw, Bw, uw
#---------------------------------------------------------------


#---------Dados de entrada---------------------------------------------------------
def propriedades_gas(P, T, dg):
    Ppc, Tpc = propriedades_pseudocriticas(dg)
    Ppr, Tpr = propriedades_pseudoreduzidas(P, T, Ppc, Tpc)
    z = fator_compressibilidade_papay(Ppr, Tpr)
    Cg, rhog, Mg = compressibilidade_gas_INSITU(P, Ppc, z, T, dg, Ppr, Tpr)
    ug = viscosidade_gas_lee(Mg, T, rhog)
    Bg = fator_formacao_gas(z, T, P)
    return Ppc, Tpc, Ppr, Tpr, z, Cg, rhog, Mg, ug, Bg

def propriedades_oleo(P, TF, dg, do, Pb, Rs, z, Bob, rhoob):

    Bg = fator_formacao_gas(z, TF, P)
    Co = compressibilidade_oleo(P, TF, dg, do, Rs, Bob, z, Api, Bg, Pb, rhoob)
    Bo = fator_volume_formacao_oleo(Rs, dg, do, TF, P, Pb, Bob, Co)
    rhoo = massa_especifica_oleo_INSITU(Rs, dg, do, Bo, P, Pb, rhoob, Co)
    uom = viscosidade_oleo_morto(Api, TF)
    uos = viscosidade_oleo_saturado(P, Pb, Rs, uom)

    return Api, Rs, Pb, Bo, rhoo, uom, uos, Co, Bg, Rs,

def propriedades_gas_unificado(P, TF_degR, dg, Mar, R):
    Ppc, Tpc = propriedades_pseudocriticas(dg)
    Ppr, Tpr = propriedades_pseudoreduzidas(P, TF_degR, Ppc, Tpc)
    z = fator_compressibilidade_papay(Ppr, Tpr)
    Mg, rhog = dados_gas(Mar, dg, P, z, R, TF_degR)
    ug = viscosidade_gas_lee(Mg, TF_degR, rhog)
    Bg = fator_formacao_gas(z, TF_degR, P)
    return Ppc, Tpc, Ppr, Tpr, z, rhog, Mg, ug, Bg

# ------------------ Bloco principal (exemplo de uso coerente) ------------------
if __name__ == "__main__":
    dg = 0.75
    qlsc_d = 10000.0  # sm3/d 
    bsw = 0.3
    RGL = 150.0
    Api = 25.0
    S = 0.0
    T=50
    TF =50
    TR=50
    Tk=50
    P = 7977.08   # psi (verifique se é psi mesmo)
    Mar = 28.96
    R = 10.7316
    d= 6 * 0.0254 # diâmetro m
    e = 3 * 10**(-6)
    T_sup=17 #C
    T1=80*(9/5) + 491.67
    T2= 4 *(9/5) + 491.67
    T3=17 *(9/5) + 491.67
    P_sc = 14.7 # Pressão na condição padrão [psia]
    T_sc = 60 # Temperatura na condição padrão [°F]
    TEC_marinho = 1 #[w/mk]
    TEC_poco = 2
    L_bomba = 1050 #m
    L_manifold = 700 #m
    TVDpoco = 450
    rho_w = 1000
    rho_ar = 1.225 #kg/m**3
    P1 = 350 * 14.504 #psi
    P1_bar = 350 #bar
    P3 = 5 * 14.504 #psi
    P3_bar = 5 #bar
    g = 9.81 #m/s**2
    T1 = 80 *(9/5) + 491.67
    T2_C = 4 #°C
    T1_F = T1 - 459.67
    T2 = 4 *(9/5) + 491.67
    T3 = 15 *(9/5) + 491.67
    sigma_wg = 0.004 #N/m
    sigma_og = 0.00841 #N/
    thata1=math.radians(60)
    theta2=math.radians(30)
    theta_3 =math.radians(90)

    

    # Gás
    Ppc, Tpc, Ppr, Tpr, z, rhog, Mg, ug, Bg = propriedades_gas_unificado(P, TR, dg, Mar, R)
    print("gás: z, rhog, Mg, ug, Bg =", z, rhog, Mg, ug, Bg)

    # Óleo - calcular do (densidade relativa) e rhoo base
    do, rhoo_kg_m3 = dados_oleo(Api)
    # Escolher: usar Pb calculado por Standing com RGL ou calcular Rs primeiro.
    Pb = obter_pressao_bolha(RGL, dg, Api, TF)
    Rs = razao_solubilidade_gas_oleo(P, dg, Api, TF, Pb)
    Bo, Bob = fator_volume_formacao_oleo(Rs, dg, do, TF, P, Pb, RGL)
    rhoo, rhoob = massa_especifica_oleo_INSITU(Rs, dg, do, Bo, P, Pb, RGL)
    uom = viscosidade_oleo_morto(Api, TF)
    uos = viscosidade_oleo_saturado(P, Pb, Rs, uom)
    print("óleo: Pb, Rs, Bo, rhoo, uom, uos =", Pb, Rs, Bo, rhoo, uom, uos)

    # Água
    rhow = massa_especifica_agua(S)
    Rsw = Razao_de_Solubilidade_agua(P, TF)
    Bw = Volume_formacao_agua(P, TF)
    uw = viscosidade_agua(P, TF, S)
    print("água: rhow, Rsw, Bw, uw =", rhow, Rsw, Bw, uw)
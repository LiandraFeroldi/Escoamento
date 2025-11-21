import math
import numpy as np
    
# Funções de Propriedades Pseudocríticas e Reduzidas
def propriedades_pseudocriticas(dg):
    if dg < 0.75:
        Ppc = 677 + 15.0 * dg - 37.5 * dg**2 #psia
        Tpc = 168 + 325 * dg - 12.5 * dg**2 #ºR
    else:
        Ppc = 706 - 51.7 * dg - 11.1 * dg**2  #psia
        Tpc = 187 + 330 * dg - 71.5 * dg**2 #ºR
    return Ppc, Tpc

def propriedades_pseudoreduzidas(P, TR, Ppc, Tpc):
    Ppr = P / Ppc
    Tpr = TR/ Tpc
    return Ppr, Tpr

def calcular_z_brill(Ppr, Tpr):
    # Fator de compressibilidade (brill e e begs)
    aa=1.39*((Tpr-0.92)**0.5)- 0.36*Tpr- 0.101
    bb=(0.62 - 0.23 * Tpr)*Ppr + ((0.066 / (Tpr - 0.86) )- 0.037) * Ppr**2 + (0.32 * Ppr**6) / (10 ** (9 * (Tpr - 1)))
    cc=0.132 - 0.32 * math.log10(Tpr)
    D_exponente = 0.3106 - 0.49 * Tpr + 0.1824 * Tpr**2
    dd= 10 ** D_exponente
    z = aa + math.exp(-bb) + cc * (Ppr ** dd)
    return z


def compressibilidade_gas_INSITU(P, Ppc, z, TR, dg, Ppr, Tpr):
    Mar = 28.96
    Mg = Mar * dg
    R = 10.7316
    rhog = (P * Mg) / (z * R * TR) 
    
    h = 0.0001 
    
    z_plus = calcular_z_brill(Ppr + h, Tpr) 
    
    # Derivada numérica: (f(x+h) - f(x)) / h
    dZ_dPpr = (z_plus - z) / h 
    # ---------------------------------------------

    dZ_dP = dZ_dPpr / Ppc
    Cg = (1 / P) - (1 / z) * dZ_dP 
    
    return Cg, rhog, Mg
#--------------------------------------------------------------------

def viscosidade_gas_lee(Mg, TR, rhog):
    kv = ((9.379 + 0.0160 * Mg) * TR ** 1.5) / (209.2 + 19.26 * Mg + TR) 
    xv = 3.448 + (986.4 / TR) + 0.01009 * Mg
    yv = 2.4 - 0.2 * xv
    ug = 1e-4 * kv * math.exp(xv * (rhog / 62.4) ** yv)
    return ug 

def fator_formacao_gas(z, TR, P):
    Bg = (14.7 / 519.67) * z * (TR / P) 
    return Bg

#-------------------------------------------------------------------------------------------------
# PROPRIEDADES DO ÓLEO
def dados_oleo(Api):
    do = (141.5 / (131.5 + Api))
    return do

def obter_pressao_bolha(RGL, dg, Api, TF):
    #Correlação de Standing
    a = 0.00091 * TF - 0.0125 * Api
    Pb = 18.2 * (((RGL / dg) ** 0.83) * 10 ** a - 1.4) 
    return Pb 

def razao_solubilidade_gas_oleo(P, dg, Api, TF, Pb):
    if P > Pb: 
        Rs = dg * ((((Pb / 18.2) + 1.4) * 10 ** (0.0125 * Api - 0.00091 * TF)) ** (1 / 0.83))
    else:
        Rs = dg * ((((P / 18.2) + 1.4) * 10 ** (0.0125 * Api - 0.00091 * TF)) ** (1 / 0.83))
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
        # Calcular Bob para garantir que rhoob seja calculado corretamente,
        # pois Bob deve usar RGL, não Rs.
        Bob_for_rhoob = 0.9759 + 0.00012 * ((RGL * ((dg / do) ** 0.5) + 1.25 * TF) ** 1.2)
        rhoob = (62.4 * do + 0.0136 * RGL * dg) / Bob_for_rhoob
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
        # P > Pb: usa viscosidade na pressão de bolha (uob)
        Rs_pb = Rs # Rs em Pb é igual a RGL
        c_pb = 10.715 * (Rs_pb + 100) ** (-0.515) 
        b_pb = 5.44 * (Rs_pb + 150) ** (-0.338)
        uob = c_pb * uom ** b_pb
        
        m = 2.6 * (P ** 1.187) * math.exp(-11.513 - (8.98e-5) * P) 
        uos = uob * (P / Pb) ** m 
    return uos

def compressibilidade_oleo(P, TF, dg, do, Rs, Bob, Api, Bg, Pb, rhoob):
    if P >= Pb:
        # Assumindo que rhoob é a densidade na pressão de bolha e Bob é o FVF na bolha
        Co = (10**-6 )* np.exp((rhoob + 0.004347 * (P - Pb) - 79.1) / (0.0007141 * (P - Pb) - 12.938))
    else:
        # Correlação de Standing para dRs/dP
        # Usamos Bob, pois Bo=Bob na região P <= Pb 
        dRs_dP = dg * (1 / 0.83) * (((P / 18.2) + 1.4) * (10 ** (0.0125 * Api - 0.00091 * TF))) ** (1 / 0.83 - 1) * (1 / 18.2) * (10 ** (0.0125 * Api - 0.00091 * TF))
        dBo_dP = 0.00012 * 1.2 * (Rs * ((dg / do) ** 0.5) + 1.25 * TF) ** 0.2 * ((dg / do) ** 0.5) * dRs_dP
        Co = -(1 / Bob) * dBo_dP + (Bg / Bob) * dRs_dP
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

# ====================================================================================

# FUNÇÃO MAIN CORRIGIDA
# Recebe todos os argumentos necessários
def main(P, TR, dg, T, R, Mar, Api, RGL, TF, S):
    # GÁS
    Ppc, Tpc = propriedades_pseudocriticas(dg)
    Ppr, Tpr = propriedades_pseudoreduzidas(P, T, Ppc, Tpc)
    z = calcular_z_brill(Ppr, Tpr)
    Cg, rhog, Mg = compressibilidade_gas_INSITU(P, Ppc, z, TR, dg, Ppr, Tpr)
    ug = viscosidade_gas_lee(Mg, TR, rhog)
    Bg = fator_formacao_gas(z, TR, P)

    # ÓLEO - FLUXO DE CÁLCULO AJUSTADO
    do = dados_oleo(Api)
    Pb = obter_pressao_bolha(RGL, dg, Api, TF)
    Rs = razao_solubilidade_gas_oleo(P, dg, Api, TF, Pb)
    uom = viscosidade_oleo_morto(Api, TF)

    # 1. Calcule Bob e rhoob (necessários para Co na região P > Pb)
    Bob_calc = 0.9759 + 0.00012 * ((RGL * ((dg / do) ** 0.5) + 1.25 * TF) ** 1.2)
    rhoob_calc = (62.4 * do + 0.0136 * RGL * dg) / Bob_calc
    
    # 2. Calcule Co (Compressibilidade do Óleo)
    Co = compressibilidade_oleo(P, TF, dg, do, Rs, Bob_calc, Api, Bg, Pb, rhoob_calc) 
    
    # 3. Use Co para calcular Bo e rhoo
    Bo, Bob = fator_volume_formacao_oleo(Rs, dg, do, TF, P, Pb, RGL, Co)
    rhoo, rhoob = massa_especifica_oleo_INSITU(Rs, dg, do, Bo, P, Pb, RGL, Co)
    uos = viscosidade_oleo_saturado(P, Pb, Rs, uom)

    # ÁGUA
    rhow = massa_especifica_agua(S)
    Rsw = Razao_de_Solubilidade_agua(P, TF)
    Bw = Volume_formacao_agua(P, TF)
    uw = viscosidade_agua(P, TF, S)

    # Retorna todos os valores calculados, incluindo Co
    return Ppc, Tpc, Ppr, Tpr, z, Cg, rhog, Mg, ug, Bg, do, Rs, Pb, Bo,Bob, rhoo,rhoob, uom, uos, Co, rhow, Rsw, Bw, uw

# ---------------------------------------------------------------


#---------Dados de entrada e Chamada da Função---------------------------------------------------------
if __name__ == "__main__":
    # DADOS DE ENTRADA DO USUÁRIO
    dg = 0.75
    RGL = 150.0
    Api = 25.0
    S = 0.0 # Salinidade
    TF = 176.0 # Temperatura em °F (Reservatório)
    TR = 176.0 + 459.67 # Temperatura em °R (para as correlações de gás)
    T = TR # T usado em propriedades_pseudoreduzidas deve ser em Rankine
    P = 7977.08 # Pressão [psia]
    Mar = 28.96 #lb/lb.mol
    R = 10.7316
    
    # Outras variáveis definidas no seu bloco original (mantidas, mas não passadas para main)
    # qlsc_d, bsw, Tk, d, e, T_sup, T1, T2, T3, P_sc, T_sc, TEC_marinho, TEC_poco, L_bomba, L_manifold, TVDpoco, rho_w, rho_ar, P1, P1_bar, P3, P3_bar, g, T2_C, T1_F, sigma_wg, sigma_og, thata1, theta2, theta_3

    # CHAMADA DA FUNÇÃO com todos os argumentos
    resultados = main(P, TR, dg, T, R, Mar, Api, RGL, TF, S)
    
    # Desempacota o resultado para o print
    (Ppc, Tpc, Ppr, Tpr, z, Cg, rhog, Mg, ug, Bg, 
     do, Rs, Pb, Bo, Bob, rhoo, rhoob, uom, uos, Co, 
     rhow, Rsw, Bw, uw) = resultados
    
    # PRINTS FINAIS
    print("=" * 60)
    print("           RESULTADOS DAS PROPRIEDADES PVT           ")
    print("=" * 60)
    
    print("\n## GÁS")
    print(f"Ppc={Ppc:.2f} psia, Tpc={Tpc:.2f} °R, Ppr={Ppr:.3f}, Tpr={Tpr:.3f}, z={z:.4f}")
    print(f"Cg={Cg:.2e} psi⁻¹, rhog={rhog:.3f} lb/ft³, Mg={Mg:.2f}lb/lb.mol, ug={ug:.3e} cP, Bg={Bg:.4e} ft³/SCF")
    
    print("\n## ÓLEO")
    print(f"do={do:.3f}, Rs={Rs:.2f} SCF/STB, Pb={Pb:.2f} psia")
    print(f"Bo={Bo:.4f} bbl/STB, Bob={Bob:.4f} bbl/STB")
    print(f"rhoo={rhoo:.3f} lb/ft³, uom={uom:.3f} cP, uos={uos:.3f} cP, Co={Co:.2e} psi⁻¹")
    
    print("\n## ÁGUA")
    print(f"rhow={rhow:.3f} lb/ft³, Rsw={Rsw:.2f} SCF/STB, Bw={Bw:.4f} bbl/STB, uw={uw:.3f} cP")
    print("-" * 60)
from barril.units import Scalar
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
#--------------------------------------------------------------------
#GÁS
def dados_gas(Mar,dg,z,R,T):
    Mg = Mar * dg
    rhog = (P * Mg) / (z * R * T) #Acho que a P e T não estão na unidade certa
    #OLHAR DEPOIS ISSO
    return Mg, rhog

def viscosidade_gas_lee(Mg, T, rhog):
    #Viscosidade do gás (Lee et al.)
    kv = ((9.379 + 0.0160 * Mg) * T ** 1.5) / (209.2 + 19.26 * Mg + T)
    xv = 3.448 + (986.4 / T) + 0.01009 * Mg
    yv = 2.4 - 0.2 * xv
    ug = 1e-4 * kv * math.exp(xv * (rhog / 62.4) ** yv)
    return ug

def fator_formacao_gas(z, TF, P):
    Bg = (14.7 / 60) * z * (TF / P)
    return Bg

#-------------------------------------------------------------------------------------------------
# PROPRIEDADES DO ÓLEO
Pb=0
def dados_oleo(Api):
    do=(141.5/(131.5+ Api))
    rhoo=do*1000 # N SEI SE A GENTE PODE USAR ASSIM PRA ACHAR O RHO DO ÓLEO
    return do, rhoo

def obter_pressao_bolha(Rs, dg, Api, TF):
    #Correlação de Standing
    a = 0.00091 * TF - 0.0125 * Api
    Pb = 18.8 * (((Rs / dg) ** 0.83) * 10 ** a - 1.4)
    return Pb

def razao_solubilidade_gas_oleo(P, dg, Api, TF, Pb):
    # Rs - Razão de solubilidade do gás no óleo
    if P >= Pb:
        Rs = dg * (((Pb / 18.2) + 1.4) * 10 ** (0.0125 * Api - 0.00091 * TF)) ** (1 / 0.83)
    else:
        Rs = dg * (((P / 18.2) + 1.4) * 10 ** (0.0125 * Api - 0.00091 * TF)) ** (1 / 0.83)
    return Rs

def fator_volume_formacao_oleo(Rs, dg, do, TF, P, Pb, Bob, Co):
    if P > Pb:
        Bo = Bob * np.exp(-Co * (P - Pb))
    else:
        Bo = 0.9759 + 0.00012 * ((Rs * ((dg / do) ** 0.5) + 1.25 * TF) ** 1.2)
    return Bo

def massa_especifica_oleo(Rs, dg, do, Bo, P, Pb, rhoob, Co):
    if P > Pb:
        rhoo = rhoob * np.exp(Co * (P - Pb))
    else:
        rhoo = (62.4 * do + 0.0136 * Rs * dg) / Bo
    return rhoo

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
        m = 2.68 * (P ** 1.187) * np.exp(-11.513 - (8.98e-5) * P)
        uos = uom * (P / Pb) ** m
    return uos

def compressibilidade_oleo(P, TF, dg, do, Rs, Bo, z, Api, Bg, Pb, rhoob):
    if P >= Pb:
        Co = 10**6 * np.exp((rhoob + 0.004347 * (P - Pb) - 79.1) / (0.0007141 * (P - Pb) - 12.938))
    else:
        dRs_dP = dg * (1 / 0.83) * (((P / 18.2) + 1.4) * (10 ** (0.0125 * Api - 0.00091 * TF))) ** (1 / 0.83 - 1) * (1 / 18.2) * (10 ** (0.0125 * Api - 0.00091 * TF))
        dBo_dP = 0.00012 * 1.2 * (Rs * ((dg / do) ** 0.5) + 1.25 * TF) ** 0.2 * ((dg / do) ** 0.5) * dRs_dP
        Co = -(1 / Bo) * dBo_dP + (Bg / Bo) * dRs_dP
    return Co

#-----------------------------------------------------------------------------------
#água
#COMPONENTES DA ÁGUA
 
def massa_especifica_agua(S):
    rhow=62.368+(0.438603)*S +(1.60074*10**-3)*(S**2) #lb/SCF
    return rhow

def Razao_de_Solubilidade_agua(P,TF):
    # Coeficientes para o parâmetro A
    A0 = 8.15839
    A1 = -6.12265e-2  # 10^-2
    A2 = 1.91663e-4   # 10^-4
    A3 = -2.1654e-7   # 10^-7
    
    # Coeficientes para o parâmetro B
    B0 = 1.01021e-2   # 10^-2
    B1 = -7.44241e-5   # 10^-5
    B2 = 3.05553e-7   # 10^-7
    B3 = -2.94883e-10  # 10^-10
    
    # Coeficientes para o parâmetro C
    C0 = -9.02505
    C1 = 0.130237
    C2 = -8.53425e-4   # 10^-4
    C3 = 2.34122e-6   # 10^-6
    C4 = -2.37049e-9   # 10^-9
    A = A0 + A1*TF + A2*math.pow(TF, 2) + A3*math.pow(TF, 3)
    
    B = B0 + B1*TF + B2*math.pow(TF, 2) + B3*math.pow(TF, 3)
    C_interno = C0 + C1*TF+ C2*math.pow(TF, 2) + C3*math.pow(TF, 3) + C4*math.pow(TF, 4)
    C = C_interno * 1e-7  # 10^-7
    Rsw = A + B*P + C*math.pow(P, 2)
    return Rsw

def Volume_formação_agua(P,TF):
    delta_VwT = (
        -1.0001e-2 +
         1.33391e-4 * TF +
         5.50654e-7 * math.pow(T, 2)
    )
    termo1 = -1.95301e-9 * P * TF
    termo2 = -1.72834e-13 * math.pow(P, 2) * TF
    termo3 = -3.58922e-7 * P
    termo4 = -2.25341e-10 * math.pow(P, 2)
    
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

def propriedades_gas(P, T, dg):
    Ppc, Tpc = propriedades_pseudocriticas(dg)
    Ppr, Tpr = propriedades_pseudoreduzidas(P, T, Ppc, Tpc)
    z = fator_compressibilidade_papay(Ppr, Tpr)
    Mg=dados_gas(Mar,dg,z,R,T)
    ug = viscosidade_gas_lee(Mg, T, rhog)
    Bg = fator_formacao_gas(z, T, P)
    return Ppc, Tpc, Ppr, Tpr, z, rhog, Mg, ug, Bg

def propriedades_oleo(P, TF, dg, do, Pb, Rs, z, Bob, rhoob):

    if Pb is None and Rs is not None:
        Pb = obter_pressao_bolha(Rs, dg, Api, TF)
    elif Rs is None and Pb is not None:
        Rs = razao_solubilidade_gas_oleo(P, dg, Api, TF, Pb)
    elif Rs is None and Pb is None:
        raise ValueError("Forneça Rs ou Pb para calcular propriedades do óleo.")

    Bg = fator_formacao_gas(z, TF, P)

    Bo = fator_volume_formacao_oleo(Rs, dg, do, TF, P, Pb, Bob)
    rhoo = massa_especifica_oleo(Rs, dg, do, Bo, P, Pb, rhoob)
    uos = viscosidade_oleo_saturado(P, Pb, Rs)

    return Api, Rs, Pb, Bo, rhoo, uos, Bg

def propriedades_agua(S,P,TF):
    rhow=massa_especifica_agua(S)
    Rsw=Razao_de_Solubilidade_agua(P,TF)
    Bw=Volume_formação_agua(P,TF)
    uw=viscosidade_agua(P,TF,S)

    return rhow, Rsw, Bw, uw
#---------------------------------------------------------------


# Vazões -------------------------------------------------------
def vazao_liquido_std(qlsc_d):
    qlsc=(10000/(24*3600)) #m³/s
    return qlsc

def vazao_std(qlsc,bsw,RGL): #superficie
    qosc=qlsc*(1- bsw)
    qwsc=qlsc*bsw
    qgsc=RGL*qlsc
    return qosc,qwsc,qgsc
    
def vazao_insitu(qosc,Bo,qwsc,Rs,Rsw,Bg,qgsc):
    ql=qosc*Bo +qwsc
    qg=(qgsc-qosc*Rs-qwsc*Rsw)*Bg
    return ql,qg
#---------Dados de entrada---------------------------------------------------------

dg = 0.75
qlsc_d=10000 #sm³/d
bsw=0.3
RGL=150 #(sm³/sm³)
Api=25
S=bsw #Salinidade em % do peso dos sólido, por aproximação
Pb = 0 # P está abaixo da pressão de bolha
TF = Scalar(122, 'degF') # Temperatura em Fahrenheit
Tc=TF.GetValue('degC')  # Converte para Celsius
T = TF.GetValue('degR')  # Converte para Rankine
P = 7977.08  # psi ou 550 bar
Mar = 28.96
R = 10.7316


# Cálculo GÁS
Ppc, Tpc, Ppr, Tpr, z, Cg, rhog, Mg, ug, Bg = propriedades_gas(P, T, dg)

# Cálculo ÓLEO
Api, Rs, Pb, Bo, rhoo, uom, uos, Co, Bg_oleo = propriedades_oleo(P, TF.GetValue('degF'), dg, do, Pb, Rs, z, Bob, rhoob)

# Cálculo ÁGUA
rhow, Rsw, Bw, uw=propriedades_agua(S,P,TF)


qlsc=vazao_liquido_std(qlsc_d)
qosc,qwsc,qgsc=vazao_std(qlsc,bsw,RGL)
ql, qg = vazao_insitu(qosc,Bo,qwsc,Rs,Rsw,Bg,qgsc)

# --- Resultados ---
print("\n--- Vazões ---")
print(f"ql= {ql:.4f} m³ | Tpc = {qg:.4f} m³")


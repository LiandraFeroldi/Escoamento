from barril.units import Scalar
import math
import numpy as np

# --- Constantes ---
rhoob = 1  
Bob = 1    

# --- Funções de propriedades do GÁS ---

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

def compressibilidade_gas(P, Ppc, z, T, dg, Ppr, Tpr):
    Mar = 28.96
    Mg = Mar * dg
    R = 10.7316
    rhog = (P * Mg) / (z * R * T)

    dZ_dPpr = -3.53 / (10 ** (0.9813 * Tpr)) + (2 * 0.274 * Ppr) / (10 ** (0.8157 * Tpr))
    dZ_dP = dZ_dPpr / Ppc

    Cg = (1 / P) - (1 / z) * dZ_dP
    return Cg, rhog, Mg

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

# --- Funções de propriedades do ÓLEO ---

def grau_api(do):
    Api = (141.5 / do) - 131.5
    return Api

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



def propriedades_gas(P, T, dg):
    Ppc, Tpc = propriedades_pseudocriticas(dg)
    Ppr, Tpr = propriedades_pseudoreduzidas(P, T, Ppc, Tpc)
    z = fator_compressibilidade_papay(Ppr, Tpr)
    Cg, rhog, Mg = compressibilidade_gas(P, Ppc, z, T, dg, Ppr, Tpr)
    ug = viscosidade_gas_lee(Mg, T, rhog)
    Bg = fator_formacao_gas(z, T, P)
    return Ppc, Tpc, Ppr, Tpr, z, Cg, rhog, Mg, ug, Bg

def propriedades_oleo(P, TF, dg, do, Pb, Rs, z, Bob, rhoob):
    Api = grau_api(do)

    
    if Pb is None and Rs is not None:
        Pb = obter_pressao_bolha(Rs, dg, Api, TF)
    elif Rs is None and Pb is not None:
        Rs = razao_solubilidade_gas_oleo(P, dg, Api, TF, Pb)
    elif Rs is None and Pb is None:
        raise ValueError("Forneça Rs ou Pb para calcular propriedades do óleo.")

    Bg = fator_formacao_gas(z, TF, P)
    Co = compressibilidade_oleo(P, TF, dg, do, Rs, Bob, z, Api, Bg, Pb, rhoob)
    Bo = fator_volume_formacao_oleo(Rs, dg, do, TF, P, Pb, Bob, Co)
    rhoo = massa_especifica_oleo(Rs, dg, do, Bo, P, Pb, rhoob, Co)
    uom = viscosidade_oleo_morto(Api, TF)
    uos = viscosidade_oleo_saturado(P, Pb, Rs, uom)

    return Api, Rs, Pb, Bo, rhoo, uom, uos, Co, Bg


# Dados de entrada
dg = 0.84 
do = 0.86
Pb = 5000
Rs = None  
Tf = Scalar(122, 'degF') # Temperatura em Fahrenheit
T = Tf.GetValue('degR')  # Converte para Rankine
P = 3626                 # Pressão em psia

# Cálculo GÁS
Ppc, Tpc, Ppr, Tpr, z, Cg, rhog, Mg, ug, Bg = propriedades_gas(P, T, dg)

# Cálculo ÓLEO
Api, Rs, Pb, Bo, rhoo, uom, uos, Co, Bg_oleo = propriedades_oleo(P, Tf.GetValue('degF'), dg, do, Pb, Rs, z, Bob, rhoob)

# --- Resultados ---
print("\n--- GÁS ---")
print(f"Ppc = {Ppc:.4f} psia | Tpc = {Tpc:.4f} °R")
print(f"Ppr = {Ppr:.4f} | Tpr = {Tpr:.4f}")
print(f"Fator z = {z:.4f}")
print(f"Compressibilidade do gás = {Cg:.4e} psia⁻¹")
print(f"Massa específica do gás = {rhog:.4f} lbm/ft³")
print(f"Viscosidade do gás = {ug:.4f} cP")
print(f"Fator de volume de formação do gás (Bg) = {Bg:.4f} ft³/scf")

print("\n--- ÓLEO ---")
print(f"API = {Api:.4f}°")
print(f"Rs = {Rs:.4f} scf/stb")
print(f"Pb = {Pb:.2f} psia")
print(f"Bo = {Bo:.4f} bbl/stb")
print(f"Massa específica do óleo = {rhoo:.4f} lbm/ft³")
print(f"Viscosidade óleo morto = {uom:.4f} cP")
print(f"Viscosidade óleo saturado = {uos:.4f} cP")
print(f"Compressibilidade do óleo = {Co:.4e} psia⁻¹")
print(f"Fator de volume de formação do gás (Bg) = {Bg_oleo:.4f} ft³/scf")
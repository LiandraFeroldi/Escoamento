from barril.units import Scalar
import math
import numpy as np
from PVT import(vazao_liquido_std,vazao_std,vazao_insitu,propriedades_agua,propriedades_oleo,propriedades_gas,Razao_de_Solubilidade_agua,Volume_formação_agua,viscosidade_agua,massa_especifica_agua,viscosidade_oleo_saturado,obter_pressao_bolha,massa_especifica_oleo,fator_volume_formacao_oleo,razao_solubilidade_gas_oleo,propriedades_pseudocriticas,propriedades_pseudoreduzidas,fator_compressibilidade_papay,dados_gas,viscosidade_gas_lee,fator_formacao_gas)
#---------Dados de entrada---------------------------------------------------------

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

#DADOS DE ENTRADA   
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
Api, Rs, Pb, Bo, rhoo, uom, uos, Co, Bg_oleo = propriedades_oleo(P, TF.GetValue('degF'), dg, do, Pb, Rs, z, Bob, rhoob)
rhow, Rsw, Bw, uw=propriedades_agua(S,P,TF)

qlsc=vazao_liquido_std(qlsc_d)
qosc,qwsc,qgsc=vazao_std(qlsc,bsw,RGL)
ql, qg = vazao_insitu(qosc,Bo,qwsc,Rs,Rsw,Bg,qgsc)

# --- Resultados ---
print("\n--- Vazões ---")
print(f"ql= {ql:.4f} m³ | Tpc = {qg:.4f} m³")
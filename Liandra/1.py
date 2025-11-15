from barril.units import Scalar
import math
import numpy as np
from l import (propriedades_pseudocriticas,propriedades_pseudoreduzidas,fator_compressibilidade_papay)

# Dados de entrada
dg = 0.75
bsw=0.3
RGL=150 #(sm³/sm³)
api=25
Pb = 5000
Rs = None  
Tf = Scalar(122, 'degF') # Temperatura em Fahrenheit
T = Tf.GetValue('degR')  # Converte para Rankine
P = 7977.08  # psi ou 550 bar
Mar = 28.96
R = 10.7316

#--------------------------------------------------------------------
do=(141.5/(131.5+ api))
Mg = Mar * dg
rhog = (P * Mg) / (z * R * T) #Acho que a P e T não estão na unidade certa
#OLHAR DEPOIS ISSO
#--------------------------------------------------------------------

Bg=0
Rs=0
Rsw=0
Bo=0
# Vazões -------------------------------------------------------
qlsc_d=10000 #sm³/d
qlsc=(10000/(60*24*60)) #m³/s
#superficie
qosc=qlsc*(1- bsw)
qwsc=qlsc*bsw
qgsc=RGL*qlsc
#in situ
ql=qosc*Bo +qwsc
qg=(qgsc-qosc*Rs-qwsc*Rsw)*Bg
#------------------------------------------------------------------

Ppc, Tpc =propriedades_pseudocriticas(dg)
Ppr, Tpr =propriedades_pseudoreduzidas(P, T, Ppc, Tpc)
z =fator_compressibilidade_papay(Ppr, Tpr) # Pela correlação de Papay
 
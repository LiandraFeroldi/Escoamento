import math
import numpy as np
T1=80  #ºC
T_inf_1=80  #ºC
T_inf_2=4   #ºC
T_inf_3=4   #ºC
T_inf_4=15  #ºC
TEC_poco=2  #w/mK
TEC_anm=1   #w/mK

Told=T1
#trecho poço
while L>=950:
    L_poco=L
    deltal=10
    v_m=rhoom*Q_m
    cp=10
    L=2100 - deltaL
    Tnew=T_inf_1 + (T_inf_1-Told) * math.exp(-TEC_poco*L_poco/ (v_m * cp)) 
    Told=Tnew
    deltaL=deltaL+10

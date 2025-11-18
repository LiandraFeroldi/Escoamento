# CÓDIGO PVT GERAL
import numpy as np
import math
import matplotlib.pyplot as plt

class PVT():
    def fase_gas(T1, P1, Z, Ppr, Tpr, rho_g1):
        dg = 0.7 # Densidade relativa do gás
        API = 30 # Grau API
        RGL = 250 # Razão gás óleo
        R = 10.73 # Constante universal dos gases [ft^3.psi/°R/lb.mol]
        P_sc = 14.7 # Pressão na condição padrão [psia]
        T_sc = 60 # Temperatura na condição padrão [°F]
        
        # Pressão de Bolha - Correlação Standing
        a = 0.00091 * (T1 - 459.67) - 0.0125 * API
        Pb = 18.2 * ((RGL / dg)**0.83) * (10**a) - 1.4
        print('O ponto de bolha é', Pb)
        
        # Cálculo Mg (massa molecular do gás)
        Mg = (0.029 * 2.20462) * dg # [lb/lbmol]
        print(Mg)
        
        # Massa específica do gás
        me_g = ((P1 * Mg) / (Z * R * T1)) * 16.018 # [kg/m^3]
        print(me_g)
        
        # Fase Gás
        # Viscosidade da fase gás
        # Correlação Lee et al.
        xv = 3.448 + (986.4 / T1) + 0.01009 * Mg
        yv = 2.4 - (0.2 * xv)
        kv = ((9.379 + 0.0160 * Mg) * T1**1.5) / (209.2 + (19.26 * Mg) + T1)
        µg = (((10**-4) * kv * np.exp(xv * ((rho_g1 / 62.4)**yv)))) # [cP]
        
        # Fator volume-formação do gás
        Bg = ((P_sc / T_sc) * (Z * (T1 - 459.67) / P1)) # [m^3/scf]
        
        # Fator Expansão do gás
        Eg = 1 / Bg
        print(f'Viscosidade do Gás: {round(µg, 5)} Pa.s')
        print(f'Bg: {round(Bg, 5)} ft^3/scf')
        print(f'Eg: {round(Eg, 5)}')
        
        return Ppr, Tpr, Z, Pb, µg, Bg, Eg, dg

    def fase_oleo(do, dg, API, T1, P1, RGL, Pb, Bg):
        # Fase Óleo
        # Razão Solubilidade
        if P1 > Pb:
            # Correlação de Standing
            Rsb = (dg * (((Pb / 18.2) + 1.4) * 10**(0.0125 * API - 0.00091 * (T1 - 459.67)))**(1 / 0.83))
            print(f'Rs: {round(Rsb, 5)} scf/stb')
            Rs = Rsb
            # Para Rs em SI:
            Rsb_si = (dg * (((Pb / 18.2) + 1.4) * 10**(0.0125 * API - 0.00091 * (T1 - 459.67)))**(1 / 0.83)) * 0.17810760667903525
            print(f'Rs_si: {round(Rsb_si, 5)} sm^3/sm^3')
            Rs_si = Rsb_si
            
            # Fator Volume-Formação do óleo
            Bob = 0.9459 + (0.00012 * (((Rs * ((dg / do)**0.5)) + 1.25 * (T1 - 459.67))**1.2)) # [bbl / stb]
            me_ob = (62.4 * do + 0.0136 * Rs * dg) / Bob
            Co = (10**(-6)) * np.exp((me_ob + 0.004347 * (P1 - Pb) - 79.1) / ((0.0007141 * (P1 - Pb)) - 12.938))
            Bo = Bob * np.exp(-Co * (P1 - Pb))
            Bt = Bo
        else:
            Rsb = (dg * (((Pb / 18.2) + 1.4) * 10**(0.0125 * API - 0.00091 * (T1 - 459.67)))**(1 / 0.83))
            # Correlação de Standing:
            Rs = (dg * (((P1 / 18.2) + 1.4) * 10**(0.0125 * API - 0.00091 * (T1 - 459.67)))**(1 / 0.83))
            print(f'Rs: {round(Rs, 5)} scf/stb')
            # Para Rs em SI:
            Rs_si = (dg * (((P1 / 18.2) + 1.4) * 10**(0.0125 * API - 0.00091 * (T1 - 459.67)))**(1 / 0.83)) * 0.17810760667903525
            print(f'Rs_si: {round(Rs_si, 5)} sm^3/sm^3')
            
            # Fator Volume-Formação do óleo com P<Pb
            # Correlação de Standing
            Bo = 0.9759 + 0.00012 * (Rs * ((dg / do)**0.5) + 1.25 * (T1 - 459.67))
            Bt = Bo + (Rsb - Rs) * Bg
            
        # Viscosidade do Óleo Morto - Correlação de Beal-Standing
        A = 10**(0.43 + (8.33 / API))
        μod = 0.32 + (1.8 * 10**7 / (API**4.53)) * ((360 / (T1 - 260))**A)
        
        if P1 <= Pb:
            # Viscosidade do Óleo Saturado P<=Pb:
            a = 10**(((-7.4 * 10**(-4)) * Rs) + ((2.2 * 10**(-7)) * (Rs**2)))
            b = (0.68 / (10**((8.62 * 10**(-5)) * Rs))) + \
                (0.25 / (10**((1.1 * 10**(-3)) * Rs))) + \
                (0.062 / (10**((3.74 * 10**(-3)) * Rs)))
            μob = a * μod**b
            μo = μob #cp
        else:
            # Viscosidade do óleo Sub-Saturado P>Pb
            a = 10**(((-7.4 * 10**(-4)) * Rs) + ((2.2 * 10**(-7)) * (Rs**2)))
            b = (0.68 / (10**((8.62 * 10**(-5)) * Rs))) + \
                (0.25 / (10**((1.1 * 10**(-3)) * Rs))) + \
                (0.062 / (10**((3.74 * 10**(-3)) * Rs)))
            μob = a * μod**b
            μo = (μob + (0.001 * (P1 - Pb)) * (0.024 * μob**1.6 + 0.038 * μob**0.56))
            
        return Rs, Bo, μo, Bt

    def fase_agua(T1, P1, bsw):
        # Fase Água
        # Aproximando S=BSW
        rho_w = (62.368 + 0.438603 * bsw + 1.60074 * 10**(-3) * bsw**2) * 16.018
        print(f'Massa Específica da água: {round(rho_w, 5)} kg/m^3')
        
        print('COEFFICIENTES')
        a0 = 8.15839
        a1 = -6.12265 * 10**(-2)
        a2 = 1.91663 * 10**(-4)
        a3 = -2.1654 * 10**(-7)
        b0 = 1.01021 * 10**(-2)
        b1 = -7.44241 * 10**(-5)
        b2 = 3.05553 * 10**(-7)
        b3 = -2.94883 * 10**(-10)
        c0 = -9.02505
        c1 = 0.130237
        c2 = -8.53425 * 10**(-4)
        c3 = 2.34122 * 10**(-6)
        c4 = -2.37049 * 10**(-9)
        
        A = a0 + a1 * (T1 - 459.67) + a2 * (T1 - 459.67)**2 + a3 * (T1 - 459.67)**3
        B = b0 + b1 * (T1 - 459.67) + b2 * (T1 - 459.67)**2 + b3 * (T1 - 459.67)**3
        C = (c0 + c1 * (T1 - 459.67) + c2 * (T1 - 459.67)**2 + c3 * (T1 - 459.67)**3 + c4 * (T1 - 459.67)**4) * 10**(-7)
        
        print(f'A = {round(A, 5)}')
        print(f'B = {round(B, 5)}')
        print(f'C = {round(C, 8)}')
        print('')
        print('RSW')
        
        RSW = A + (B * P1) + (C * (P1)**2)
        RSW = max(RSW, 0) # Verificação para garantir que RSW não seja negativo
        print(f'Razão Solubilidade Gás - Água: {round(RSW, 5)}')
        
        dVwt = -1.0001 * (10**-2) + 1.33391 * (10**-4) * (T1 - 459.67) + 5.50654 * (10**-7) * (T1 - 459.67)**2
        dVwp = -1.95301 * (10**-9) * P1 * (T1 - 459.67) - 1.72834 * (10**-13) * (T1 - 459.67) * (P1**2) - 3.58922 * (10**-7) * P1 - 2.25341 * (10**-10) * (P1**2)
        
        Bw = ((1 + dVwt) * (1 + dVwp))
        print(f'Fator Volume Formação Água: {round(Bw, 5)} bbl/STB')
        
        # Viscosidade da Água
        a0 = 109.527
        a1 = -8.40564
        a2 = 0.313314
        a3 = 8.72213 * 10**(-3)
        b0 = -1.12166
        b1 = 2.63951 * 10**(-2)
        b2 = -6.79461 * 10**(-4)
        b3 = -5.47119 * 10**(-5)
        b4 = -1.55586 * 10**(-6)
        
        A = a0 + a1 * bsw + a2 * bsw**2 + a3 * bsw**3
        B = b0 + b1 * bsw + b2 * bsw**2 + b3 * bsw**3 + b4 * bsw**4
        u_a1 = A * (T1 - 459.67)**B
        print(f'Viscosidade da Água - Temperatura: {round(u_a1, 5)} Cp')
        
        # Efeito da Pressão na água - Viscosidade
        μw = (u_a1 * (0.994 + 4.0295 * 10**(-5) * P1 + 3.1062 * 10**(-9) * P1**2))
        print(f'Viscosidade da Água: {round(μw, 8)} cP')
        
        return rho_w, RSW, Bw, μw
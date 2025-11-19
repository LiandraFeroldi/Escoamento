# CÓDIGO PVT RESERVATÓRIO
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

class pvt_reser():
    def calculate_pvt_reser():
        dg = 0.7 # Densidade relativa do gás
        API = 30 # Grau API
        RGL = 250 # Razão gás óleo
        R = 10.73 # Constante universal dos gases [ft^3.psi/°R/lb.mol]
        M_ar = 0.029 # Massa molecular do ar [kg/mol]
        P_sc = 14.7 # Pressão na condição padrão [psia]
        T_sc = 60 # Temperatura na condição padrão [°F]
        TEC_marinho = 1 #[w/mk]
        TEC_poco = 2
        T_sup = 15
        L_bomba = 1050 #m
        L_manifold = 700 #m
        Z_poco = 450
        d = 3.5 * 0.0254 # diâmetro máxima para obter surgencia
        e = 3 * 10**(-6)
        rho_w = 1000
        rho_ar = 1.225 #kg/m**3
        P1 = 350 * 14.504 #psi
        P1_bar = 350 #bar
        P3 = 5 * 14.504 #psi
        P3_bar = 5 #bar
        bsw = 25
        g = 9.81 #m/s**2
        T1 = 80 * (9/5) + 491.67
        T2_C = 4 #°C
        T1_F = T1 - 459.67
        T2 = 4 * (9/5) + 491.67
        T3 = 15 * (9/5) + 491.67
        σwg = 0.004 #N/m
        σog = 0.00841 #N/m
        0_rad1 = math.radians(60)
        0_rad2 = math.radians(30)
        theta_3 = math.radians(90)
        print(T1, T2, T3)
        print(P1)
        print(T1_F)
        
        # Área de passagem
        Ap = np.pi * (d / 2)**2
        print(f'Área de Passagem: {round(Ap, 5)} m^2')
        
        # Densidade do óleo
        do = 141.5 / (API + 131.5)
        print(f'Densidade do óleo: {round(do, 5)} Kg/m^3')
        rho_o = rho_w * do
        print(f'Massa específica do óleo: {round(rho_o, 5)}')
        
        # Densidade do gás
        rho_g = dg * rho_ar
        print(rho_g, 'massa específica do gás')
        Q_I = 10000 / 86400 # [m^3/s]
        print(f'Vazão Volumétrica: {round(Q_I, 5)} kg/m^3')
        q_m = Q_I * do
        print(f'Vazão Mássica: {round(q_m, 5)} kg/m^3')
        Cp = ((2 * 10**-3) * T_sup - 1.429) * do + (2.67 * 10**-3) * T_sup + 3.049
        print(f'Calor Específico Óleo: {round(q_m, 5)} kg/m^3') # O print usa q_m, pode ser um typo no original
        
        L_poco = Z_poco / np.sin(teta_rad1)
        print(f'Poço - Manifold: {round(L_poco, 5)} m')
        v = Q_I / Ap
        print(v, 'velocidade')
        coef_ang = (T2 - T1) / 450
        print(f'Coeficiente: {round(coef_ang, 5)} m')
        
        P_pc = 677 + 15.0 * dg - 37.5 * dg**2
        T_pc = 168 + 325 * dg - 12.5 * dg**2
        print(T_pc, P_pc)
        
        # Inclinado
        n = 500
        d_L = L_poco / n
        print(f'Distância entre pontos Poço: {round(d_L, 5)} m')
        d_T = (T1 - T2) / n
        print(f'Distância entre pontos Temperatura: {round(d_T, 5)} °R')
        x = np.linspace(0, Z_poco, n)
        T_novo = np.full_like(x, T1)
        T_old = [T1]
        T_L = []
        for i in range(len(x)):
            T_novo[i] -= d_T * i
            T_L_value = T_novo[i] - ((q_m * g * np.sin(teta_rad1)) / TEC_poco) * (np.exp((-TEC_poco / (q_m * Cp)) * L_poco)) * (T_novo[i] - (q_m * g * np.sin(teta_rad1)) / TEC_poco - T_old[i])
            T_L.append(T_L_value)
            T_old_value = T_L_value
            T_old.append(T_old_value)
        print(T_L)
        
        # Reduzida
        T_pr_inc = []
        for i in range(len(x)):
            T_novo[i] -= d_T * i
            T_pr_ = T_L[i] / T_pc
            T_pr_inc.append(T_pr_)
        print(T_pr_inc)
        
        #Trecho Horizontal
        d_T = T2
        x_t = np.linspace(0, L_manifold, n)
        T_novo = np.full_like(x_t, T2)
        T_old = [T2]
        T_L2 = []
        for i in range(len(x_t)):
            T_novo[i] = d_T
            T_L2.append(T_novo[i])
        print("Valores de temperatura ao longo do trecho horizontal:", T_L2)
        
        # Reduzida
        T_pr_hor = []
        for i in range(len(x)):
            T_novo[i] = d_T
            T_pr_2 = T_L2[i] / T_pc
            T_pr_hor.append(T_pr_2)
        print(T_pr_hor)
        
        # VERTICAL
        d_T3 = (T3 - T2) / n
        y = np.linspace(Z_poco, L_bomba, n)
        T_novo = np.full_like(y, T2)
        T_old = [T2] # Inicializar T_old com T2
        T_L3 = []
        # Ajustar o cálculo da temperatura para garantir que termine em T1
        for i in range(len(y)):
            T_novo[i] += d_T3 * i
            T_L_value = T_novo[i] - ((q_m * g) / TEC_marinho) * (np.exp((-TEC_marinho / (q_m * Cp)) * L_bomba)) * (T_novo[i] - (q_m * g) / TEC_marinho - T_old[i])
            T_L3.append(T_L_value)
            T_old_value = T_L_value
            T_old.append(T_old_value)
        print("T_L3:", T_L3)
        
        # Temperatura reduzida
        T_pr_vert = []
        for T in T_L3:
            T_pr_ = T / T_pc
            T_pr_vert.append(T_pr_)
        print("T_pr_vert:", T_pr_vert)
        
        #gas
        # Cálculo do Z (Papay)
        Ppr = P1 / P_pc
        Tpr = T1 / T_pc
        print(Ppr, Tpr)
        Z = 1 - ((3.53 * Ppr) / (10**((0.9813) * Tpr))) + ((0.274 * Ppr**2) / (10**(0.8157 * Tpr)))
        print(Z) # primeiro z (reservatório)
        
        # Pressão de Bolha - Correlação Standing
        a = 0.00091 * T1_F - 0.0125 * API
        Pb = 18.2 * ((RGL / dg)**0.83) * (10**a) - 1.4
        print('O ponto de bolha é', Pb)
        
        # Cálculo Mg (massa molecular do gás)
        Mg = (0.029 * 2.20462) * dg # [lb/lbmol]
        print(Mg)
        
        # Massa específica do gás
        #me_g =((P1*Mg)/(Z*R*T1))* 16.018 # [kg/m^3]
        #print(me_g)
        
        # Fase Gás
        rho_g1 = rho_g / 16.01846337 # densidade do gás em lb/ft^3
        
        # Viscosidade da fase gás
        # Correlação Lee et al.
        xv = 3.448 + (986.4 / T1) + 0.01009 * Mg
        yv = 2.4 - (0.2 * xv)
        kv = ((9.379 + 0.0160 * Mg) * T1**1.5) / (209.2 + (19.26 * Mg) + T1)
        µg = (((10**-4) * kv * np.exp(xv * ((rho_g1 / 62.4)**yv)))) # [cP]
        # µg_si = (((10**-4) * kv * np.exp(xv* ((rho_g1/62.4) ** yv))))*0.001 # [Pa.s]
        
        # Fator volume-formação do gás
        Bg = ((P_sc / T_sc) * (Z * T1_F / P1)) # [m^3/scf]
        
        # Fator Expansão do gás
        Eg = 1 / Bg
        print(f'Viscosidade do Gás: {round(µg, 5)} cP')
        print(f'Bg: {round(Bg, 5)} ft^3/scf')
        print(f'Eg: {round(Eg, 5)}')
        
        # Fase Óleo
        # Razão Solubilidade
        # P > Pb
        # Rs = Rsb
        # Correlação de Standing
        Rs = (dg * (((Pb / 18.2) + 1.4) * 10**(0.0125 * API - 0.00091 * T1_F))**(1 / 0.83))
        print(f'Rs: {round(Rs, 5)} scf/stb')
        
        # Para Rs em SI:
        Rs_si = (dg * (((Pb / 18.2) + 1.4) * 10**(0.0125 * API - 0.00091 * T1_F))**(1 / 0.83)) * 0.17810760667903525
        print(f'Rs_si: {round(Rs_si, 5)} sm^3/sm^3')
        
        # Compressibilidade Isotérmica - Correlação Vasquez e Beggs
        Co = (-1433 + 5 * Rs + 17.2 * T1_F - 1180 * dg + 12.61 * API) / (10**5 * P1)
        print(f'Co: {round(Co, 8)} psi^(-1)')
        
        # Para SI:
        Co_si = Co * (1 / 6894.76)
        print(f'Co_si: {round(Co_si, 10)} Pa^(-1)')
        
        # Fator Volume-Formação do Gás (Óleo)
        Bob = 0.9759 + (0.00012 * (((Rs * ((dg / do)**0.5)) + 1.25 * T1_F)**1.2)) # [bbl/stb]
        me_ob = (62.4 * do + 0.0136 * Rs * dg) / Bob
        Co = (10**(-6)) * np.exp((me_ob + 0.004347 * (P1 - Pb) - 79.1) / ((0.0007141 * (P1 - Pb)) - 12.938))
        Bo = Bob * np.exp(-Co * (P1 - Pb))
        Bt = Bo # Bt=Bo (P>Pb)
        print(f'Bob: {round(Bob, 5)} bbl/stb')
        print(f'Co: {round(Co, 8)} psi^(-1)')
        print(f'Bo: {round(Bo, 5)} bbl/stb')
        
        # Viscosidade do Óleo Morto - Correlação de Beal-Standing
        A = 10**(0.43 + (8.33 / API))
        # Cálculo da viscosidade do óleo morto
        μod = 0.32 + (1.8 * 10**7 / (API**4.53)) * ((360 / (T1 - 260))**A)
        
        # Cálculo da variável 'a'
        a = 10**(((-7.4 * 10**(-4)) * Rs) + ((2.2 * 10**(-7)) * (Rs**2)))
        # Cálculo da variável 'b'
        b = (0.68 / (10**((8.62 * 10**(-5)) * Rs))) + (0.25 / (10**((1.1 * 10**(-3)) * Rs))) + (0.062 / (10**((3.74 * 10**(-3)) * Rs)))
        
        # Cálculo da viscosidade do óleo
        μob = a * μod**b
        # Cálculo da viscosidade do óleo
        μo = μob + (0.001 * (P1 - Pb)) * (0.024 * μob**1.6 + 0.038 * μob**0.56)
        # Cálculo da viscosidade do óleo
        μo_si = (μob + (0.001 * (P1 - Pb)) * (0.024 * μob**1.6 + 0.038 * μob**0.56)) * (10**(-3))
        
        print(f'Viscosidade do Óleo Morto: {round(μod, 5)} cP')
        print(f'Viscosidade do Óleo: {round(μo, 5)} cP')
        print(f'Viscosidade do Óleo: {round(μo_si, 5)} Pa.s')
        
        # Fase Água
        #aproximando S=BSW
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
        
        A = a0 + a1 * T1_F + a2 * T1_F**2 + a3 * T1_F**3
        B = b0 + b1 * T1_F + b2 * T1_F**2 + b3 * T1_F**3
        C = (c0 + c1 * T1_F + c2 * T1_F**2 + c3 * T1_F**3 + c4 * T1_F**4) * 10**(-7)
        print(f'A = {round(A, 5)}')
        print(f'B = {round(B, 5)}')
        print(f'C = {round(C, 8)}')
        print(' ')
        print('RSW')
        
        RSW = A + B * P1 + C * P1**2
        print(f'Razão Solubilidade Gás - Água: {round(RSW, 5)}')
        dVwt = -1.0001 * (10**-2) + 1.33391 * (10**-4) * T1_F + 5.50654 * (10**-7) * (T1_F**2)
        dVwp = -1.95301 * (10**-9) * P1 * T1_F - 1.72834 * (10**-13) * T1_F * (P1**2) - 3.58922 * (10**-7) * P1 - 2.25341 * (10**-10) * (P1**2)
        
        Bw = ((1 + dVwt) * (1 + dVwp)) * 0.158987294928 # [m^3/sm^3]
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
        u_a1 = A * T1_F**B
        print(f'Viscosidade da Água - Temperatura: {round(u_a1, 5)} Cp')
        
        # Efeito da Pressão na água - Viscosidade
        μw = (u_a1 * (0.994 + 4.0295 * 10**(-5) * P1 + 3.1062 * 10**(-9) * P1**2))
        print(f'Viscosidade da Água: {round(μw, 8)} cP')
        μw_si = (u_a1 * (0.994 + 4.0295 * 10**(-5) * P1 + 3.1062 * 10**(-9) * P1**2)) * 0.001
        print(f'Viscosidade da Água: {round(μw_si, 8)} Pa.s')
        
        return T_L, T_L2, T_L3, T_pr_hor, T_pr_vert, T_pr_inc, e, d, Q_I, Bo, Bw, RGL, Rs, \
               RSW, Bg, rho_w, rho_o, rho_g1, σwg, σog, μw_si, μo_si, g, teta_rad1, rho_g, µg, P1, L_poco, \
               theta_3, P3, L_bomba, L_manifold, P1_bar, n, P_pc, P3_bar, do, API, bsw, Z # Z_poco foi trocado por Z

    def calculate_bb(T_L, T_L2, T_L3, T_pr_hor, T_pr_vert, T_pr_inc, e, d, Q_I, Bo, Bw,
                     RGL, RS, RSW, Bg, rho_w, rho_o, σwg, σog, μw_si, μo_si, g, teta_rad1, rho_g, µg, P1, L_poco,
                     theta_3, P3, L_bomba, L_manifold, P1_bar, n, P_pc, P3_bar, Z_poco):
        
        # Modelo Beggs & Brill
        BSW_1 = 0.25
        # Rugosidade Relativa
        e_d = e / d
        # área de passagem
        Ap = (np.pi * d**2) / 4
        # Vazão de óleo na superfície
        Qosc = Q_I * (1 - BSW_1)
        # Vazão de água na superfície
        Qwsc = Q_I * BSW_1
        # Vazão insitu da água e do óleo
        Qo_insitu = Qosc * Bo
        Qw_insitu = Qwsc * Bw
        # Vazão do Líquido insitu
        Ql_insitu = Qosc * Bo + Qwsc * Bw
        # Vazão de gás na superfície
        Qgsc = RGL * Q_I
        # Vazão insitu gás
        Qg_insitu = (Qgsc - Qosc * RS - Qwsc * RSW) * Bg
        # Velocidade superfícial do líquido
        Vsl = Ql_insitu / Ap
        # Velocidade superfícial do gás
        Vsg = Qg_insitu / Ap
        # Velocidade da mistura
        Vm = Vsl + Vsg
        # Houldup do líquido (lambda_l)
        xl = Vsl / Vm
        # Fração volumétrica da Água
        Fwc = Qw_insitu / (Qo_insitu + Qw_insitu)
        # Fração volumétrica do Óleo
        Foc = 1 - Fwc
        # Densidade do líquido
        rhol = rho_w * Fwc + rho_o * Foc
        # Tensão interfacial do líquido
        σl = σwg * Fwc + σog * Foc
        # Viscosidade do líquido
        μl = μw_si * Fwc + μo_si * Foc
        # Número de Froude da mistura
        Frm = Vm**2 / (g * d)
        
        print(f'Vazão Insitu Óleo: {round(Qo_insitu, 5)} m^3/s')
        print(f'Vazão Insitu Gás: {round(Qg_insitu, 5)} m^3/s')
        print(f'Vazão Insitu Água: {round(Qw_insitu, 5)} m^3/s')
        print(f'Vazão Insitu Líquido: {round(Ql_insitu, 5)} m^3/s')
        print(f'Velocidade da Superfície Óleo: {round(Vsl, 5)} m/s')
        print(f'Velocidade da Superfície Gás: {round(Vsg, 5)} m/s')
        print(f'Velocidade da Mistura: {round(Vm, 5)} m/s')
        print(f'Houldup do Líquido: {round(xl, 5)}')
        print(f'Fração Volumétrica da Água: {round(Fwc, 5)}')
        print(f'Fração Volumétrica do Óleo: {round(Foc, 5)}') # Gás no original, mas Foc é Óleo
        print(f'Densidade do Líquido: {round(rhol, 5)} kg/m^3')
        print(f'Tensão Interfacial do Líquido: {round(σl, 5)} N/m')
        print(f'Viscosidade do Líquido: {round(μl, 5)} Pa.s')
        print(f'Número de Froude da Mistura: {round(Frm, 5)}')
        
        # O modelo para estes resultados é INTERMITENTE.
        # Parâmetros L
        L1 = 316 * (xl**0.302)
        L2 = 0.0009252 * (xl**(-2.4684))
        L3 = 0.10 * (xl**(-1.4516))
        L4 = 0.5 * (xl**(-6.738))
        
        # Resultados
        print(f'Parâmetros L1: {round(L1, 5)}')
        print(f'Parâmetros L2: {round(L2, 5)}')
        print(f'Parâmetros L3: {round(L3, 5)}')
        print(f'Parâmetros L4: {round(L4, 5)}')
        
        # Dados (para Intermitente)
        a = 0.845
        b = 0.5351
        c = 0.0173
        
        # Holdup de líquido
        Hlo = (a * xl**b) / (Frm)**c
        Nlv = Vsl * ((rhol / (g * σl))**(1 / 4))
        C = (1 - xl) * (np.log(2.960 * (xl**0.305) * (Nlv**(-0.4473)) * (Frm)**0.0978))
        ψ = 1 + C * (np.sin(1.8 * teta_rad1) - 0.333 * (np.sin(1.8 * teta_rad1))**3)
        Hl = Hlo * ψ
        
        print(f'Holdup Líquido: {round(Hl, 5)}')
        print(f'Holdup horizontal: {round(Hlo, 5)}')
        print(f'Número de velocidades do Líquido: {round(Nlv, 5)}')
        print(f'Parâmetro de correção da inclinação: {round(ψ, 5)}')
        
        rho_ns = rhol * xl + rho_g * (1 - xl)
        μ_ns = μl * xl + μg * (1 - xl)
        
        Re_ns = (rho_ns * Vm * d) / μ_ns
        print(f'Densidade: {round(rho_ns, 5)}')
        print(f'Viscosidade: {round(μ_ns, 5)} Pa.s')
        print(f'Número de Reynolds: {round(Re_ns, 5)}')
        
        # Cálculo de Fn pela correlação de Hall:
        f_n = 0.0055 * (1 + ((2 * 10**4) * (e_d)) + (10**6 / Re_ns))**(1 / 3)
        print("O fator de atrito é:", f_n)
        
        # Parâmetro Y:
        y = xl / (Hlo)**2
        print("Parâmetro Y:", y)
        
        # Para y > 1.2
        # Parâmetro S:
        s = (math.log(y)) / (-0.0523 + 3.18 * math.log(y) - 0.8725 * (math.log(y))**2 + 0.01853 * (math.log(y))**4)
        print("Parâmetro S:", s)
        
        f_tp = f_n * np.exp(s)
        print("Fator de Atrito:", f_tp)
        
        rho_slip = (rhol * Hlo) + (rho_g * (1 - Hlo))
        print("Densidade da Mistura Slip:", rho_slip, "kg/m^3")
        
        # Seção 1 - Inclinada
        # Gradiente de Pressão de Fricção:
        dp_dl_f1 = f_tp * ((rho_ns * Vm**2) / (2 * d))
        print("Gradiente de Pressão de Fricção:", dp_dl_f1, "N/m")
        # Gradiente de Pressão Gravitacional:
        dp_dl_g1 = rho_slip * g * np.sin(teta_rad1)
        print("Gradiente de Pressão Gravitacional:", dp_dl_g1, "N/m")
        # Parâmetro Ek:
        Ek1 = ((rho_slip * Vm * Vsg)) / (P1 * 6894.76)
        print("Parâmetro Ek:", Ek1)
        # Desconsideraremos o parâmetro Ek, pois ele é de ordem muito pequena.
        # Gradiente de Pressão Total:
        dp_dl_t1 = (dp_dl_f1 + dp_dl_g1)
        print("Gradiente de Pressão Total:", dp_dl_t1, "N/m")
        # Gradiente de Pressão por Aceleração:
        dp_dl_a1 = Ek1 * (dp_dl_t1)
        print("Gradiente de Pressão de Aceleração:", dp_dl_a1, "N/m")
        # O gradiente de pressão por aceleração será desconsiderado por ser um valor pequeno.
        
        dp1 = dp_dl_t1 * L_poco
        print("Diferença de Pressão:", dp1, "Pa")
        P_suc = (P1 * 6894.76) - dp1
        print("Pressão de Entrada na Bomba:", P_suc, "Pa")
        P_suc_bar = P_suc * (10**-5)
        
        # Seção 2 - Horizontal
        # Gradiente de Pressão de Fricção:
        dp_dl_f2 = f_tp * ((rho_ns * Vm**2) / (2 * d))
        print("Gradiente de Pressão de Fricção:", dp_dl_f2, "Pa/m")
        # Gradiente de Pressão Total:
        dp_dl_t2 = (dp_dl_f2)
        print("Gradiente de Pressão Total:", dp_dl_t2, "'Pa/m")
        
        # Seção 3 - Vertical
        d1 = 6.2 * 0.0254
        # Gradiente de Pressão de Fricção:
        dp_dl_f3 = f_tp * ((rho_ns * Vm**2) / (2 * d1))
        print("Gradiente de Pressão de Fricção:", dp_dl_f3, "Pa/m")
        # Gradiente de Pressão Gravitacional:
        dp_dl_g3 = rho_slip * g * np.sin(theta_3)
        print("Gradiente de Pressão Gravitacional:", dp_dl_g3, "Pa/m")
        # Parâmetro Ek:
        Ek3 = ((rho_slip * Vm * Vsg)) / (P3 * 6894.76)
        print("Parâmetro Ek:", Ek3)
        # Desconsideraremos o parâmetro Ek, pois ele é de ordem muito pequena.
        # Gradiente de Pressão Total:
        dp_dl_t3 = (dp_dl_f3 + dp_dl_g3)
        print("Gradiente de Pressão Total:", dp_dl_t3, "Pa/m")
        # Gradiente de Pressão por Aceleração:
        dp_dl_a3 = Ek3 * (dp_dl_t3)
        print("Gradiente de Pressão de Aceleração:", dp_dl_a3, "Pa/m")
        # O gradiente de pressão por aceleração será desconsiderado por ser um valor pequeno.
        
        dp3 = dp_dl_t3 * L_bomba
        print("Diferença de Pressão:", dp3, "Pa")
        P2 = (dp3 + (P3 * 6894.76))
        P2_bar = P2 * 10**(-5)
        print("Pressão no manifold:", P2_bar, "Pa")
        
        # Pressão da Bomba
        dp2 = dp_dl_f2 * L_manifold
        print(dp2)
        P_desc = dp2 + P2
        print("Pressão de Descarga na Bomba:", P_desc, "Pa")
        P_desc_bar = P_desc * (10**-5)
        P_bomba = P_desc - P_suc
        print("Pressão da Bomba:", P_bomba, "Pa")
        
        # Discretização da Malha para Pressão
        #Trecho inclinado
        d_P = (P1_bar - P_suc_bar) / n
        x_p = np.linspace(0, L_poco, n)
        P_novo = np.full_like(x_p, P1_bar)
        P_old = [P1_bar]
        P_L1 = []
        for i in range(len(x_p)):
            P_novo[i] -= d_P * i
            P_L1.append(P_novo[i])
        print("Valores de pressão ao longo do trecho inclinado:", P_L1)
        
        # Reduzida
        P_pr_inc = []
        Z1_inc = []
        for i in range(len(x_p)):
            # P_novo[i] -= d_P*i # Esta linha parece duplicada/errada no original
            P_pr_ = P_L1[i] / P_pc
            P_pr_inc.append(P_pr_)
            Z_inc = 1 - ((3.53 * P_pr_inc[i]) / (10**((0.9813) * T_pr_inc[i]))) + ((0.274 * P_pr_inc[i]**2) / (10**(0.8157 * T_pr_inc[i])))
            Z1_inc.append(Z_inc)
        print("Valores de Pressão Reduzida para o trecho inclinado: ", P_pr_inc)
        print("Valores de Z para o trecho inclinado:", Z1_inc)
        
        #Trecho Horizontal
        d_P2 = (P_desc_bar - P2_bar) / n
        x_p2 = np.linspace(Z_poco, L_manifold, n)
        P_novo = np.full_like(x_p2, P_desc_bar)
        P_old = [P_desc_bar]
        P_L2 = []
        for i in range(len(x_p2)):
            P_novo[i] -= d_P2 * i
            P_L2.append(P_novo[i])
        print("Valores de pressão ao longo do trecho horizontal:", P_L2)
        
        # Reduzida
        P_pr_hor = []
        Z2_hor = []
        for i in range(len(x_p2)):
            # P_novo[i] -= d_P2*i # Esta linha parece duplicada/errada no original
            P_pr_2 = P_L2[i] / P_pc
            P_pr_hor.append(P_pr_2)
            Z_hor = 1 - ((3.53 * P_pr_hor[i]) / (10**((0.9813) * T_pr_hor[i]))) + ((0.274 * P_pr_hor[i]**2) / (10**(0.8157 * T_pr_hor[i])))
            Z2_hor.append(Z_hor)
        print("Valores de Pressão Reduzida para o trecho horizontal:", P_pr_hor)
        print("Valores de Z para o trecho horizontal:", Z2_hor)
        
        # Trecho vertical
        d_P3 = (P2_bar - P3_bar) / n
        x_p1 = np.linspace(Z_poco, L_bomba, n)
        P_novo = np.full_like(x_p1, P2_bar)
        P_old = [P2_bar]
        P_L3 = []
        for i in range(len(x_p1)):
            P_novo[i] -= d_P3 * i
            P_L3.append(P_novo[i])
        print("Valores da Pressão ao longo do trecho vertical:", P_L3)
        
        # Reduzida
        P_pr_vert = []
        Z3_vert = []
        for i in range(len(x_p1)):
            # P_novo[i] -= d_P3*i # Esta linha parece duplicada/errada no original
            P_pr_1 = P_L3[i] / P_pc
            P_pr_vert.append(P_pr_1)
            Z_vert = 1 - ((3.53 * P_pr_vert[i]) / (10**((0.9813) * T_pr_vert[i]))) + ((0.274 * P_pr_vert[i]**2) / (10**(0.8157 * T_pr_vert[i])))
            Z3_vert.append(Z_vert)
        print("Valores de Pressão Reduzida para o trecho vertical:", P_pr_vert)
        print("Valores de Z para o trecho vertical:", Z3_vert)
        
        return P_L1, P_L2, P_L3, P_pr_hor, P_pr_inc, P_pr_vert, Z2_hor, Z3_vert, Z1_inc
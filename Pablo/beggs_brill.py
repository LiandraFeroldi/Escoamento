# beggs_brill.py
import math

g = 9.81   # gravidade

# ----------------------------------------------------------
#  TABELAS DE COEFICIENTES (do material / slide)
# ----------------------------------------------------------

ABC_TABLE = {
    "segregated":    (0.98,  0.4846, 0.0868),
    "intermittent":  (0.845, 0.5351, 0.0173),
    "distributed":   (1.065, 0.5824, 0.0609)
}

CEFG_TABLE_ASC = {
    "segregated":    (0.011, -3.768,  3.539,  -1.614),
    "intermittent":  (2.960, 0.305,  -0.4473, 0.0978),
    "distributed":   (0.011, -3.768,  3.539,  -1.614)   # mas AGORA não será usado (C=0)
}


# ----------------------------------------------------------
# 1. Velocidades superficiais
# ----------------------------------------------------------
def velocities_from_q(Q_st, D, Bo, Rs, RGL, Bg, BSW):
    A = math.pi * D**2 / 4.0

    Q_o_st = Q_st * (1 - BSW)
    Q_w_st = Q_st * BSW

    Ql = Q_o_st * Bo + Q_w_st * 1.0
    Qg = max((RGL * Q_st - Q_o_st * Rs) * Bg, 0.0)

    Vsl = Ql / (A + 1e-12)
    Vsg = Qg / (A + 1e-12)
    Vm  = Vsl + Vsg + 1e-12

    return Vsl, Vsg, Vm


# ----------------------------------------------------------
# 2. Parâmetros L1, L2, L3, L4
# ----------------------------------------------------------
def froude_mixture_squared(Vm, D):
    return Vm**2 / (g * D + 1e-12)

def lambda_no_slip(Vsl, Vm):
    return max(0.0, min(Vsl / (Vm + 1e-12), 1.0))

def L_params(lambda_L):
    L1 = 316.0      * (lambda_L**0.302)
    L2 = 0.0009252  * (lambda_L**-2.4684)
    L3 = 0.1        * (lambda_L**-1.4516)
    L4 = 0.5        * (lambda_L**-6.738)
    return L1, L2, L3, L4


# ----------------------------------------------------------
# 3. Seleção do regime (exatamente como o slide)
# ----------------------------------------------------------
def determine_regime(lambda_L, Frm2, L1, L2, L3, L4):
    if (lambda_L < 0.4 and Frm2 >= L1) or (lambda_L >= 0.4 and Frm2 > L4):
        return "distributed"

    if (lambda_L < 0.01 and Frm2 < L1) or (lambda_L >= 0.001 and Frm2 < L2):
        return "segregated"

    if (0.01 <= lambda_L < 0.4 and L3 <= Frm2 <= L1) or (lambda_L >= 0.4 and L3 <= Frm2 <= L4):
        return "transition"

    return "intermittent"


# ----------------------------------------------------------
# 4. Holdup básico HLo (a,b,c)
# ----------------------------------------------------------
def compute_HLo(lambda_L, Frm2, regime):
    a, b, c = ABC_TABLE[regime]
    HLo = (a * lambda_L**b) / (Frm2**c + 1e-12)

    # regra: HLo sempre > lambda_L
    if HLo <= lambda_L:
        HLo = lambda_L + 1e-6

    return max(min(HLo, 0.9999), 1e-6)


# ----------------------------------------------------------
# 5. C parameter (corrigido → C = 0 se distributed)
# ----------------------------------------------------------
def compute_C_parameter(lambda_L, Vsl, Frm2, rho_l, sigma_gl, regime, is_descending):
    """
    Correções aplicadas:
      ✔ se regime == distributed → C = 0 (psi = 1)
      ✔ se descendo → C = 0 (psi = 1)
    """
    if regime == "distributed" or is_descending:
        return 0.0

    dprime, e, f, gexp = CEFG_TABLE_ASC[regime]

    denom = g * sigma_gl + 1e-12
    Nlv = Vsl * (rho_l / denom)**0.25

    inside = max(dprime * (lambda_L**e) * (Nlv**f) * (Frm2**gexp), 1e-12)
    C = (1 - lambda_L) * math.log(inside)
    return C


def psi_inclination(C, theta_deg):
    th = math.radians(theta_deg)
    return 1.0 + C * (math.sin(1.8 * th) - 0.333 * (math.sin(1.8 * th)**3))


# ----------------------------------------------------------
# 6. HL para regime de transição (A * HL_seg + (1-A) HL_int)
# ----------------------------------------------------------
def compute_transition_HL(lambda_L, Frm2, Vsl, theta_deg, rho_l, sigma_gl):
    L1, L2, L3, L4 = L_params(lambda_L)

    A = (L3 - Frm2) / (L3 - L2 + 1e-12)
    A = max(min(A, 1.0), 0.0)

    # HL segregado
    HLo_seg = compute_HLo(lambda_L, Frm2, "segregated")
    C_seg = compute_C_parameter(lambda_L, Vsl, Frm2, rho_l, sigma_gl, "segregated", False)
    psi_seg = psi_inclination(C_seg, theta_deg)
    HL_seg = HLo_seg * psi_seg

    # HL intermitente
    HLo_int = compute_HLo(lambda_L, Frm2, "intermittent")
    C_int  = compute_C_parameter(lambda_L, Vsl, Frm2, rho_l, sigma_gl, "intermittent", False)
    psi_int = psi_inclination(C_int, theta_deg)
    HL_int  = HLo_int * psi_int

    HL = A * HL_seg + (1 - A) * HL_int
    return max(min(HL, 0.9999), 1e-6)


# ----------------------------------------------------------
# 7. Fator de atrito (Hall)
# ----------------------------------------------------------
def friction_factor_hall(eps, D, Re):
    term = 2e4 * (eps / (D + 1e-12)) + 1e6 / (Re + 1e-12)
    return 0.0055 *( 1 + (term ** (1/3)))


# ----------------------------------------------------------
# 8. Função principal do Beggs & Brill
# ----------------------------------------------------------
def beggs_brill(Q_st, D, eps, theta_deg, prop):
    """
    prop deve conter:
      Bo, Bg, Rs, RGL, BSW,
      rho_l, rho_g,
      mu_l, mu_g,
      sigma_gl
    """

    # Velocidades
    Vsl, Vsg, Vm = velocities_from_q(
        Q_st, D,
        prop['Bo'], prop['Rs'], prop['RGL'], prop['Bg'], prop['BSW']
    )

    lambda_L = lambda_no_slip(Vsl, Vm)
    Frm2 = froude_mixture_squared(Vm, D)
    L1, L2, L3, L4 = L_params(lambda_L)
    regime = determine_regime(lambda_L, Frm2, L1, L2, L3, L4)
    is_descending = (theta_deg < 0)

    # HL
    if regime == "transition":
        HL = compute_transition_HL(
            lambda_L, Frm2, Vsl, theta_deg,
            prop['rho_l'], prop['sigma_gl']
        )
    else:
        HLo = compute_HLo(lambda_L, Frm2, regime)
        C = compute_C_parameter(
            lambda_L, Vsl, Frm2,
            prop['rho_l'], prop['sigma_gl'],
            regime, is_descending
        )
        psi = psi_inclination(C, theta_deg)
        HL = max(min(HLo * psi, 0.9999), 1e-6)

    # Mistura
    rho_ns = prop['rho_l'] * lambda_L + prop['rho_g'] * (1 - lambda_L)
    rho_slip = prop['rho_l'] * HL + prop['rho_g'] * (1 - HL)

    mu_ns = prop['mu_l'] * lambda_L + prop['mu_g'] * (1 - lambda_L)

    Re_ns = max(rho_ns * Vm * D / (mu_ns + 1e-12), 1e-9)

    f_n = friction_factor_hall(eps, D, Re_ns)

    # y e s (dois ramos)
    y = lambda_L / (HL**2 + 1e-12)

    if 1 < y < 1.2:
        s = math.log(2.2 * y - 1.2 + 1e-12)
    else:
        ly = math.log(max(y, 1e-12))
        den = -0.0523 + 3.182 * ly - 0.8725 * ly**2 + 0.01853 * ly**4
        s = (ly / den) if abs(den) > 1e-12 else 0.0

    f_tp = f_n * math.exp(s)

    # Perdas por fricção e gravidade
    dp_fric = f_tp * rho_ns * Vm**2 / (2 * D)

    sin_theta = math.sin(math.radians(theta_deg))
    dp_grav = rho_slip * g * sin_theta

    # Função que retorna dp/dL
    def compute_dp(P_abs):
        Ek = (rho_slip * Vm * Vsg) / (P_abs + 1e-12)
        denom = max(1 - Ek, 1e-6)
        dp_total = -(dp_fric + dp_grav) / denom   # negativo ascendente
        return dp_total, HL, regime, {
            "dp_fric": dp_fric,
            "dp_grav": dp_grav,
            "Re_ns": Re_ns,
            "Vm": Vm,
            "f_n": f_n,
            "f_tp": f_tp,
            "y": y,
            "s": s
        }

    return compute_dp


prop = {
    "Bo": 1.2,
    "Bg": 0.004,
    "Rs": 60,
    "RGL": 150,
    "BSW": 0.30,
    "rho_l": 850,
    "rho_g": 12,
    "mu_l": 0.001,
    "mu_g": 1e-5,
    "sigma_gl": 0.02
}

compute_dp = beggs_brill(
    Q_st = 10000/86400,
    D = 0.089,
    eps = 3e-6,
    theta_deg = 37,
    prop = prop
)

dp_dl, HL, regime, info = compute_dp(500*1e5)

print("dp/dL:", dp_dl)
print("HL:", HL)
print("Regime:", regime)
print("Detalhes internos:", info)

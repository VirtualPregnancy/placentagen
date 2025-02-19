import numpy as np
import scipy as sp

def exchange_maternal_fetal_oxygen_no_vessel_resistance(Q_m, Q_f, C_ma, C_fa, N_co, N_sa, verbose=False):
    """
    Parameters
    ----------
    Q_m Maternal flow
    Q_f Fetal flow
    C_ma Concentration in the maternal arterial blood supply
    C_fa Concentration in fetal arterial blood
    N_co Number of cotyledons
    N_sa Number of Spiral arteries

    Returns
    -------
    C_fv Concentration in fetal venous circulation
    C_mv Concentration in maternal venous circulation

    """
    N_tv = N_co*N_sa
    
    D_t = 2 # Diffusion coefficient for terminal villous tissue in m^2/s
    L_tv = 16.5*10**-6 # Length scale for TV geom, roughly exchange area divided by thickness, 15 - 18 mm
    F = lambda theta: 1 - np.e**(-1*theta)

    Damkohler_fetal = (D_t*L_tv*N_tv)/Q_f
    Damkohler_maternal = ((D_t*L_tv*N_tv)/Q_m) * (F(Damkohler_fetal)/Damkohler_fetal)

    N_tot = (F(Damkohler_fetal)/Damkohler_fetal) * (F(Damkohler_maternal)/Damkohler_maternal) * D_t * L_tv * N_tv * (C_ma - C_fa)
    if verbose:
        print(f"{N_tot:.3} oxygen transferred, with fetal Damkohler: {Damkohler_fetal:.4}, maternal Damkohler: "
              f"{Damkohler_maternal:.4} and {N_tv} terminal villi")
    C_fv = (N_tot + Q_f*C_fa)/Q_f
    C_mv = (Q_m*C_ma - N_tot)/Q_m
    return C_fv, C_mv

def exchange_maternal_fetal_oxygen_with_vessel_resistance(P_m, P_f, C_ma, C_fa, N_co, N_sa, R_uta, R_utv, R_co, R_fpa, R_fpv, R_tv, verbose=False):
    """
    Parameters
    ----------
    P_m Pressure drop across maternal circulation from maternal arterial to maternal venous
    P_f Pressure drop across fetal circulation from fetal arterial to fetal venous
    C_ma Concentration in the maternal arterial blood supply
    C_fa Concentration in fetal arterial blood
    N_co Number of cotyledons
    N_sa Number of Spiral arteries
    R_uta Resistance to flow in a typical uterine arterial network
    R_utv Resistance to flow in a typical uterine venous network
    R_co Resistance to flow across the intervillous space in a single cotyledon from spiral artery to decidual vein
    R_fpa Resistance to flow in a typical fetal arterial network
    R_fpv Resistance to flow in a typical fetal venous network
    R_tv Resistance to flow of a single capillary network in a terminal villae

    Returns
    -------
    C_fv Concentration in fetal venous circulation
    C_mv Concentration in maternal venous circulation

    """
    N_tv = N_co*N_sa

    Q_m = P_m/(R_uta + R_co + R_utv)
    Q_f = P_f/(R_fpa + R_fpv + R_tv)

    D_t = 2 # Diffusion coefficient for terminal villous tissue in m^2/s
    L_tv = 16.5*10**-6 # Length scale for TV geom, roughly exchange area divided by thickness, 15 - 18 mm
    F = lambda theta: 1 - np.e**(-1*theta)

    Damkohler_fetal = (D_t*L_tv*N_tv)/Q_f
    Damkohler_maternal = ((D_t*L_tv*N_tv)/Q_m) * (F(Damkohler_fetal)/Damkohler_fetal)

    N_tot = (F(Damkohler_fetal)/Damkohler_fetal) * (F(Damkohler_maternal)/Damkohler_maternal) * D_t * L_tv * N_tv * (C_ma - C_fa)
    if verbose:
        print(f"{N_tot:.3} oxygen transferred, with fetal Damkohler: {Damkohler_fetal:.4}, maternal Damkohler: "
              f"{Damkohler_maternal} and {N_tv} terminal villi")
    C_fv = (N_tot + Q_f*C_fa)/Q_f
    C_mv = (Q_m*C_ma - N_tot)/Q_m
    return C_fv, C_mv

def consumption(t, C):
    """
    :param t: time
    :param C: Concentration of Oxygen in Fetal circulation
    :return: dodt - rate of change in fetal oxygen conctration with respect to time
    This function models the change in fetal oxygen concentration as a function of time using Michaelis-Menten mechanics
    TODD: properly parameterise Max_consumption and K_m
    """
    Max_consumption = 0.1
    K_m = 0.044#ml/ml
    dodt = -Max_consumption*C/(K_m + C) #ml/ml/aec
    if C<=0:
        dodt=0.0

    return dodt

def oxygen_consumption(C_fetal, cardiac_cyle_time):
    """
    :param C_fetal: Initial concentration of fetal oxygen
    :param cardiac_cyle_time: time it takes for the fetus to complete one cardiac cycle
    :return: two array-like objects with indexed values representing the time course of fetal oxygen concentration, the
    first array is time, and the second is oxygen concentration
    """
    sol= sp.integrate.solve_ivp(consumption, [0, cardiac_cyle_time], [C_fetal,])
    return sol.t, sol.y

def convert_po2_to_concentration(p_o2):
    """
    :param po2: partial pressure of oxygen in blood
    :return: C_o2, the concentration of oxygen in the blood in ml/ml
    """
    # concentration of oxygen in plasma given by 3 * 10^-5 * po2 (in mmhg) gives concentration in ml/ml
    C_o2_plasma = 3 * 10 ** -5 * p_o2
    # log pO2 = k1 − k2(pH − 7.4) + k3log(SHb/(100 − SHb)), solve for S_hb
    k1 = 1.445
    k2 = 0.456
    k3 = 0.371
    # o2_capacity = amout of oxygen that can be carried by haemoglobin (1.34 ml/g) * haemoglobin in the blood (~0.125 g/ml)
    Hb_cap = 1.34
    C_Hb = 0.125

    # Concentration of oxygen bound to haemoglobin ( C_Hb)
    # = (S_hb * o2_capacity)/100
    w = p_o2 ** (1 / k3) + np.e ** (-k1 / k3)
    S_Hb = (100 * w) / (1 + w)  # Rearrangement of modified hills equation from Mabelle Lins thesis to solve for SH

    O2_cap = Hb_cap * C_Hb

    C_o2_Hb = S_Hb * O2_cap * (1 / 100)
    C_o2 = C_o2_Hb + C_o2_plasma
    return C_o2
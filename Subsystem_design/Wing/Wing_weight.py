from Subsystem_design.common_constants import Constants
import math

def Torenbeek_method():
    c = Constants()
    kg_to_lbs = 2.20462
    m_to_ft = 3.28084

    W_MZF = c.MZFW_320neo       # [kg]
    # W_MZF = W_MZF*kg_to_lbs

    sweep = c.sweep_05          # [deg]
    sweep = sweep * math.pi/180

    b = c.b                     # [m]
    # b = b*m_to_ft

    n_ult = 2.5 * 1.5
    S = 112.4                           # [m^2]
    # S = S * m_to_ft**2

    cr = c.c_root                       # [m] root chord
    # cr = cr * m_to_ft

    tc_max = 0.15                       # NACA 23015 Amir taper ratio
    tr = cr*tc_max                      # maximum thickness of root chord
    # tr = tr * m_to_ft

    # Ww = 0.0017 * W_MZF * (b/math.cos(sweep))**(0.75) * ( 1 + (6.3* math.cos(sweep)/b)**(1/2))* n_ult**(0.55)* (b*S/tr*W_MZF*math.cos(sweep))**(0.30) #lbs

    Ww = 6.67*10**(-3) * W_MZF * (b / math.cos(sweep)) ** (0.75) * (1 + math.sqrt(1.905 * math.cos(sweep) / b)) * \
         n_ult ** (0.55) * (b * S / tr * W_MZF * math.cos(sweep)) ** (0.30)  # kg
    # Ww = Ww / kg_to_lbs
    print(Ww)

def ADSEE_method():
    c = Constants()
    W = c.MTOW_320neo
    n = 1.5 * 2.5
    Sw = 112.4
    A = 9.39
    t_c = 0.15
    taper = 0.24
    sweep = 25 * math.pi/180
    Sc = 45.12 #control surface area #Jane
    Ww = 0.0051*pow(W*n,0.557)*pow(Sw,0.649)*pow(A,0.5)* pow(t_c,-0.4)*pow(1+taper,0.1)*pow(math.cos(sweep),-1)*pow(Sc,0.1)
    print(Ww)

# Torenbeek_method()
ADSEE_method()
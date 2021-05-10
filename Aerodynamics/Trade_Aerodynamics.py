import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import math as m
import drag_coefficient_estimation_Roskam


def fus_wet_surface(l_cockpit, l_cabin, l_tail, df):
    """
    :param l_cockpit: length of the cockpit in [m]
    :param l_cabin: length of the cabin in [m]
    :param l_tail: length of the tail in [m]
    :param df: diameter length in [m]
    :return: Wet surface of the fuselage in [m^2]
    """

    S_w_fus = m.pi * df/4 * (1/3/ l_cockpit **2 * ((4*l_cockpit**2+df**2/4)-df**3/8) - df + 4 * l_cabin + 2 * m.sqrt(l_tail**2+ df**2/4))

    return S_w_fus

#######INSERT HERE DESIRED PARAMETERS

visc = 1.458 * 10**(-5) # N*s/m^2
air_d = 0.4135 # kg/m^3
u1 = 0.7 * np.sqrt(air_d * 287 * (273.15 - 50)) # m/s
M = 0.7

S = 122.6 # m^2
lf = 37.57 #fuselage length in [m]
df = 4.14  #fuselage maximum diameter in [m]
new_df = 5 #modified maximum diameter in [m]
l_cockpit = 5.04 #[m]
l_cabin = 29.53 - l_cockpit # [m]
l_tail = lf - 29.53 #[m]
S_b_fus = np.pi * 0.3/2 * 0.45/2

ct = 1.64 # [m] obtained from top view drawing
cr = 6.07 # [m] pbtained from top view drawing
sweep = 25 # [deg] wing sweep obtained from Elsevier data base



""" Concept 1: Cargo Hold"""
# Aerodynamic impact is 0


""" Concept 3: Longer Fuselage """
l_extra = 1 # m
S_fus = np.pi * (df/2)**2
S_wet_fus =
roskam = drag_coefficient_estimation_Roskam.Roskam_drag_coefficient(visc=visc, u1=u1, air_d=air_d, l_f=lf+l_extra, M=M,
                                                                    S=S)
R_wf, C_f_fus, C_D_o_fus = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin+l_extra,
                                                                        l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                        S_wet_fus=)

""" Concept 4: Flat Bottom """
S_fus = 0.5 * np.pi * (df/2)**2 + 2 * (2.07 * 2.17 - 0.6**2 + 0.25 * np.pi * 0.6**2)
S_wet_fus =
roskam = drag_coefficient_estimation_Roskam.Roskam_drag_coefficient(visc=visc, u1=u1, air_d=air_d, l_f=lf, M=M, S=S)
R_wf, C_f_fus, C_D_o_fus = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                        l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                        S_wet_fus=)

""" Concept 5: Wing Pods """

pods_number = 2
pods_cdo = 0.025
pods_interf = 1.3

pods_L1 = 0.5
pods_main_length = 5
pods_L3 = 0.5

S_w_pods = fus_wet_surface(pods_L1, pods_main_length, pods_L3)
pods_contrib_cdo = pods_cdo*pods_interf*S_w_pods #todo: divide by wet surface of the plane and by Cdo of A320

""" Concept 6: Beluga """
S_fus =
S_wet_fus =
roskam = drag_coefficient_estimation_Roskam.Roskam_drag_coefficient(visc=visc, u1=u1, air_d=air_d, l_f=lf, M=M, S=S)
R_wf, C_f_fus, C_D_o_fus = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                        l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                        S_wet_fus=)







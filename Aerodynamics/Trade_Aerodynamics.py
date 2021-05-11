import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import math as m
import drag_coefficient_estimation_Roskam
import S_wet_estimation


# def fus_wet_surface(l_cockpit, l_cabin, l_tail, df):
#     """
#     :param l_cockpit: length of the cockpit in [m]
#     :param l_cabin: length of the cabin in [m]
#     :param l_tail: length of the tail in [m]
#     :param df: diameter length in [m]
#     :return: Wet surface of the fuselage in [m^2]
#     """
#
#     S_w_fus = m.pi * df/4 * (1/3/ l_cockpit **2 * ((4*l_cockpit**2+df**2/4)-df**3/8) - df + 4 * l_cabin + 2 * m.sqrt(l_tail**2+ df**2/4))
#
#     return S_w_fus

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
S_fus = np.pi * (df/2)**2

S_wet_fus_st = S_wet_estimation.S_wet_estimation_standard(l_cockpit=l_cockpit,l_cabin=l_cabin,l_tail=l_tail,df1=df,
                                                          df2=df)
S_wet_fus_st.calculate_volume()
S_wet_fus_st.S_wet()
S_wet_fus = S_wet_fus_st.S_wet_fus

roskam = drag_coefficient_estimation_Roskam.Roskam_drag_coefficient(visc=visc, u1=u1, air_d=air_d, l_f=lf, M=M,
                                                                    S=S)
R_wf_st, C_f_fus_st, C_D_o_fus_st = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                        l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                        S_wet_fus=S_wet_fus)
print('\n Concept 1: Cargo Hold')
print(C_D_o_fus_st, '\n \n')


""" Concept 3: Longer Fuselage """
l_extra = 1.06 # m #todo: fill the right value
S_fus = np.pi * (df/2)**2

S_wet_fus_lg = S_wet_estimation.S_wet_estimation_standard(l_cockpit=l_cockpit,l_cabin=l_cabin+l_extra,l_tail=l_tail,
                                                          df1=df, df2=df)
S_wet_fus_lg.calculate_volume()
S_wet_fus_lg.S_wet()
S_wet_fus = S_wet_fus_lg.S_wet_fus

roskam = drag_coefficient_estimation_Roskam.Roskam_drag_coefficient(visc=visc, u1=u1, air_d=air_d, l_f=lf+l_extra, M=M,
                                                                    S=S)
R_wf_lg, C_f_fus_lg, C_D_o_fus_lg = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin+l_extra,
                                                                        l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                        S_wet_fus=S_wet_fus)
print('Concept 3: Longer Fuselage')
print(C_D_o_fus_lg, '\n \n')

""" Concept 4: Flat Bottom """

S_wet_fus_belly = S_wet_estimation.S_wet_estimation_belly(l_cockpit=l_cockpit,l_cabin=l_cabin,l_tail=l_tail,
                                                          df1=df, dh=0.1, rb=0.6)
S_wet_fus_belly.calculate_volume()
S_wet_fus_belly.S_wet()
S_fus = S_wet_fus_belly.cross_section
S_wet_fus = S_wet_fus_belly.S_wet_fus

roskam = drag_coefficient_estimation_Roskam.Roskam_drag_coefficient(visc=visc, u1=u1, air_d=air_d, l_f=lf, M=M, S=S)
R_wf_belly, C_f_fus_belly, C_D_o_fus_belly = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                        l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                        S_wet_fus=S_wet_fus)
print('Concept 4: Flat Bottom')
print(C_D_o_fus_belly, '\n \n')

""" Concept 5: Wing Pods """

N_pods = 2
C_D_o_pods_init = 0.025
Interf_pods = 1.3

L1_pods = 0.5
L2_pods= 5
L3_pods = 0.01
d_pods = 1.6 #[m]

S_wet_podded = S_wet_estimation.S_wet_estimation_standard(l_cockpit=L1_pods, l_cabin=L2_pods, l_tail=L3_pods, df1=d_pods,
                                                        df2=d_pods)
S_wet_podded.calculate_volume()
S_wet_podded.S_wet()
S_wet_pods = S_wet_podded.S_wet_fus

C_D_o_pods = N_pods * C_D_o_pods_init * Interf_pods * S_wet_pods/S_wet_fus
#todo: divide by wet surface of the plane and by Cdo of A320

print('Concept 5: Wing Pods')
print(C_D_o_pods)
print(C_D_o_pods/C_D_o_fus_st, ' %')
print('Plane CD_0 = ', C_D_o_pods+C_D_o_fus_st, '\n \n')

""" Concept 6: Beluga """

S_wet_fus_beluga = S_wet_estimation.S_wet_estimation_beluga(l_cockpit=5.04,l_cabin=24.49,l_tail=8.04,beluga=35,
                                                            df=4.14,dfb=4.14*0.5)
S_wet_fus_beluga.calculate_volume()
S_wet_fus_beluga.S_wet()
S_fus = S_wet_fus_beluga.S_beluga
S_wet_fus = S_wet_fus_beluga.S_wet_fus

roskam = drag_coefficient_estimation_Roskam.Roskam_drag_coefficient(visc=visc, u1=u1, air_d=air_d, l_f=lf, M=M, S=S)
R_wf_beluga, C_f_fus_beluga, C_D_o_fus_beluga = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit,
                                                                                             l_cabin=l_cabin,
                                                                        l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                        S_wet_fus=S_wet_fus)
print('Concept 6: Beluga')
print(C_D_o_fus_beluga, '\n \n')







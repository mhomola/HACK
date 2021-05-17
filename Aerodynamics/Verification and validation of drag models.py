import drag_coefficient_estimation_Roskam
import Trade_Aerodynamics
import S_wet_estimation
import numpy as np

# ===== INSERT HERE DESIRED PARAMETERS =====

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

# ===== Unit tests on dimensions and wetted area =====

print('\n===== Unit tests on dimensions and wetted area =====')

# 0. Standard configuration

S_wet_fus_st = S_wet_estimation.S_wet_estimation_standard(l_cockpit=l_cockpit, l_cabin=l_cabin, l_tail=l_tail, df1=df,
                                                          df2=df)
S_wet_fus_st.calculate_volume()
S_wet_fus_st.S_wet()
S_fus = np.pi * (df/2)**2
V_cpit = S_wet_fus_st.V_cockpit
V_cab = S_wet_fus_st.V_cabin
V_tail = S_wet_fus_st.V_tail
S_wet_fus = S_wet_fus_st.S_wet_fus
S_wet_fus_2 = Trade_Aerodynamics.fus_wet_surface(l_cockpit, l_cabin, l_tail, df)

roskam = drag_coefficient_estimation_Roskam.Roskam_drag_coefficient(visc=visc, u1=u1, air_d=air_d, l_f=lf, M=M,
                                                                    S=S)
_, _, C_D_o_fus_st = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                                 l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                                 S_wet_fus=S_wet_fus)
_, C_D_b, C_D_o_fus_st2 = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                                 l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                                 S_wet_fus=S_wet_fus_2)

print('\n For the standard configuration: '
      '\n The volume fo the cockpit is = ', V_cpit, ' m^3',
      '\n The volume of the cabin is = ', V_cab, ' m^3',
      '\n The volume of the tail is = ', V_tail, ' m^3',
      '\n Total volume of the fuselage is = ', V_cpit + V_cab + V_tail, ' m^3',
      '\n Wetted area of the fuselage according to regression = ', S_wet_fus, 'm^2',
      '\n Wetted area of the fuselage according to ADSEE simplification = ', S_wet_fus_2, 'm^2',
      '\n C_D_0 using wetted area from regression = ', C_D_o_fus_st,
      '\n C_D_0 using adsee wetted area = ', C_D_o_fus_st2,
      '\n C_D_b using adsee wetted area = ', C_D_b)
# 1. Elongated fuselage configuration

S_wet_fus_lg = S_wet_estimation.S_wet_estimation_standard(l_cockpit=l_cockpit, l_cabin=l_cabin + 2.76, l_tail=l_tail,
                                                          df1=df, df2=df)
S_wet_fus_lg.calculate_volume()
S_wet_fus_lg.S_wet()
S_fus = np.pi * (df/2)**2
V_cpit = S_wet_fus_lg.V_cockpit
V_cab = S_wet_fus_lg.V_cabin
V_tail = S_wet_fus_lg.V_tail
S_wet_fus = S_wet_fus_lg.S_wet_fus
S_wet_fus_2 = Trade_Aerodynamics.fus_wet_surface(l_cockpit, l_cabin + 2.76, l_tail, df)

roskam = drag_coefficient_estimation_Roskam.Roskam_drag_coefficient(visc=visc, u1=u1, air_d=air_d, l_f=lf+2.76, M=M,
                                                                    S=S)
_, _, C_D_o_fus_lg = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit,
                                                                                 l_cabin=l_cabin+2.76,
                                                                                 l_tail=l_tail, S_fus=S_fus,
                                                                                 S_b_fus=S_b_fus,
                                                                                 S_wet_fus=S_wet_fus)
_, _, C_D_o_fus_lg2 = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit,
                                                                                  l_cabin=l_cabin+2.76,
                                                                                  l_tail=l_tail, S_fus=S_fus,
                                                                                  S_b_fus=S_b_fus,
                                                                                  S_wet_fus=S_wet_fus_2)

print('\n For the longer configuration: '
      '\n The volume fo the cockpit is = ', V_cpit, ' m^3',
      '\n The volume of the cabin is = ', V_cab, ' m^3',
      '\n The volume of the tail is = ', V_tail, ' m^3',
      '\n Total volume of the fuselage is = ', V_cpit + V_cab + V_tail, ' m^3',
      '\n Wetted area of the fuselage according to regression = ', S_wet_fus, 'm^2',
      '\n Wetted area of the fuselage according to ADSEE simplification = ', S_wet_fus_2, 'm^2',
      '\n C_D_0 using wetted area from regression = ', C_D_o_fus_lg,
      '\n C_D_0 using adsee wetted area = ', C_D_o_fus_lg2,
      '\n', (C_D_o_fus_lg/C_D_o_fus_st - 1)*100, 'CD_0 increase in % for regression method',
      '\n', (C_D_o_fus_lg2/C_D_o_fus_st2 - 1)*100, 'CD_0 increase in % for adsee case', '\n \n')

# 2. Flat belly configuration

S_wet_fus_belly = S_wet_estimation.S_wet_estimation_belly(l_cockpit=l_cockpit,l_cabin=l_cabin,l_tail=l_tail,
                                                          df1=df, dh=0.1, rb=0.55, pc=12.25/l_cabin)
S_wet_fus_belly.calculate_volume()
S_wet_fus_belly.S_wet()
V_cab = S_wet_fus_belly.V_cabin
S_fus = S_wet_fus_belly.cross_section
S_wet_fus = S_wet_fus_belly.S_wet_fus
S_wet_fus_2 = Trade_Aerodynamics.fus_wet_surface(l_cockpit, l_cabin, l_tail, np.sqrt(4 * S_fus / np.pi))

_, _, C_D_o_fus_belly = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                                 l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                                 S_wet_fus=S_wet_fus)
_, _, C_D_o_fus_belly2 = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                                  l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                                  S_wet_fus=S_wet_fus_2)

print('\n For flat belly configuration: '
      '\n Surface are of the fuselage cross section = ', S_fus, ' m^2',
      '\n The volume of the cabin is = ', V_cab, ' m^3',
      '\n Total volume of the fuselage is = ', V_cpit + V_cab + V_tail, ' m^3',
      '\n Wetted area of the fuselage according to regression = ', S_wet_fus, 'm^2',
      '\n Wetted area of the fuselage according to ADSEE simplification = ', S_wet_fus_2, 'm^2',
      '\n C_D_0 using wetted area from regression = ', C_D_o_fus_belly,
      '\n C_D_0 using adsee wetted area = ', C_D_o_fus_belly2,
      '\n', (C_D_o_fus_belly/C_D_o_fus_st - 1)*100, 'CD_0 increase in % for regression method',
      '\n', (C_D_o_fus_belly2/C_D_o_fus_st2 - 1)*100, 'CD_0 increase in % for adsee case', '\n \n')

# 3. Wing podded configuration

N_pods = 2
C_D_o_pods_init = 0.025
Interf_pods = 1.3
L1_pods = 0.5
L2_pods = 4.5
L3_pods = 1
d_pods = 2.1  # [m]

S_wet_fus_pod = S_wet_estimation.S_wet_estimation_standard(l_cockpit=L1_pods, l_cabin=L2_pods, l_tail=L3_pods,
                                                           df1=d_pods, df2=d_pods)
S_wet_fus_pod.calculate_volume()
S_wet_fus_pod.S_wet()
V_1 = S_wet_fus_pod.V_cockpit
V_2 = S_wet_fus_pod.V_cabin
V_3 = S_wet_fus_pod.V_tail
S_wet_pods = S_wet_fus_pod.S_wet_fus
S_wet_pods_2 = Trade_Aerodynamics.fus_wet_surface(L1_pods, L2_pods, L3_pods, d_pods)

C_D_o_pods = N_pods * C_D_o_pods_init * Interf_pods * S_wet_pods_2/S_wet_fus
C_D_o_pods2 = N_pods * C_D_o_pods_init * Interf_pods * S_wet_pods_2/S_wet_fus_2
_, _, C_D_o_pods3 = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=L1_pods, l_cabin=L2_pods,
                                                                   l_tail=L3_pods, S_fus=np.pi/4*d_pods**2, S_b_fus=0,
                                                                   S_wet_fus=S_wet_pods_2)

print('\n For wing podded configuration: '
      '\n Total volume of one podded tank is = ', V_1 + V_2 + V_3, ' m^3',
      '\n Wetted area of one tank according to regression = ', S_wet_pods, 'm^2',
      '\n Wetted area of one tank according to ADSEE simplification = ', S_wet_pods_2, 'm^2'
      '\n Fuselage + tanks CD_0 using the regression method = ', C_D_o_pods + C_D_o_fus_st,
      '\n Fuselage + tanks CD_0 using the adsee and Torenbeek method = ', C_D_o_pods2, '+', C_D_o_fus_st2, '=',
      C_D_o_pods2 + C_D_o_fus_st2,
      '\n Fuselage + tanks CD_0 using the adsee method and Roskam = 2 x ', C_D_o_pods3, '+', C_D_o_fus_st2, '=',
      2 * C_D_o_pods3 + C_D_o_fus_st2,
      '\n', (C_D_o_pods/C_D_o_fus_st)*100, 'CD_0 increase in % for regression method',
      '\n', (C_D_o_pods2/C_D_o_fus_st2)*100, 'CD_0 increase in % for adsee adn Torenbeek method'
      '\n', (2*C_D_o_pods3/C_D_o_fus_st2)*100, 'CD_0 increase in % for adsee and Roskam method \n \n')

# 4. Beluga configuration

S_wet_fus_beluga = S_wet_estimation.S_wet_estimation_beluga(l_cockpit=5.04,l_cabin=24.49,l_tail=8.04,beluga=100,
                                                            df=4.14,dfb=1.968,h1 =5.283,h2=3.8187)
S_wet_fus_beluga.calculate_volume()
S_wet_fus_beluga.S_wet()
S_fus = S_wet_fus_beluga.S_beluga
V_cpit = S_wet_fus_beluga.V_cockpit
V_cab = S_wet_fus_beluga.V_beluga_cabin
V_tail = S_wet_fus_beluga.V_tail
S_wet_fus = S_wet_fus_beluga.S_wet_fus
S_wet_fus_2 = Trade_Aerodynamics.fus_wet_surface(l_cockpit, l_cabin, l_tail, np.sqrt(4 * S_fus / np.pi))

_, _, C_D_o_fus_beluga = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                                 l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                                 S_wet_fus=S_wet_fus)
_, _, C_D_o_fus_beluga2 = roskam.run_Roskam_drag_coefficient_functions(l_cockpit=l_cockpit, l_cabin=l_cabin,
                                                                                  l_tail=l_tail, S_fus=S_fus, S_b_fus=S_b_fus,
                                                                                  S_wet_fus=S_wet_fus_2)

print('\n For beluga configuration: '
      '\n Surface are of the fuselage cross section = ', S_fus, ' m^2',
      '\n The volume fo the cockpit is = ', V_cpit, ' m^3',
      '\n The volume of the cabin is = ', V_cab, ' m^3',
      '\n The volume of the tail is = ', V_tail, ' m^3',
      '\n Total volume of the fuselage is = ', V_cpit + V_cab + V_tail, ' m^3',
      '\n Wetted area of the fuselage according to regression = ', S_wet_fus, 'm^2',
      '\n Wetted area of the fuselage according to ADSEE simplification = ', S_wet_fus_2, 'm^2',
      '\n C_D_0 using wetted area from regression = ', C_D_o_fus_beluga,
      '\n C_D_0 using adsee wetted area = ', C_D_o_fus_beluga2,
      '\n', (C_D_o_fus_beluga/C_D_o_fus_st - 1)*100, 'CD_0 increase in % for regression method',
      '\n', (C_D_o_fus_beluga2/C_D_o_fus_st2 - 1)*100, 'CD_0 increase in % for adsee method')




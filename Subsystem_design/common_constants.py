import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from math import pi
from Subsystem_design.fuel_required import V_H2, V_k
from Subsystem_design.Tank_Design.Main_PreliminaryTank import mass_pod, mass_center_tank, volume_pod, volume_centre_tank
from Subsystem_design.Engine.EnergySplit import LHV_hack, ER_h2, ER_ker, MR_h2, MR_ker



"""
This file contains one class only, which is meant to contain the variables which are common to the entire subsystem
design. It may also contain some simple functions to compute constants derived from other constants (e.g. ISA).
"""


class Constants():
    def __init__(self):
        self.rho_0 = 1.225                                          # Sea level density                         [kg/m^3]
        self.p_0 = 101325                                           # Sea level pressure                           [Pa]
        self.T_0 = 288.15                                           # Sea level temperature                        [K]
        self.g_0 = 9.80665                                          # Gravity at sea level                       [m/s^2]
        self.R = 287.0                                              # Specific gas constant of air            [J/(kg*K)]
        self.R_univ = 8.314                                         # Universal gas constant                    [J/mol/K]
        self.gamma = 1.4                                            # Heat capacity ratio of air                   [-]
        self.a_0 = 340.294                                          # Sea level speed of sound                     [m/s]
        self.visc = 1.458 * 10**(-5)                                # Air viscosity                            [N*s/m^2]


        """Properties of H2"""
        self.H2_ed = 33.5                                       # Energy density of hydrogen                   [kWh/kg]
        self.LHV_h2 = 119.96                                    # Lower heating value hydrogen                 [MJ/kg]
        self.rho_h2 = 0.07111                                   # Density of hydrogen                          [kg/L]

        """Properties of kerosene"""
        self.k_ed = 12.0                                    # Energy density of kerosene                       [kWh/kg]
        self.LHV_ker = 43.2                                 # Lower heating value of kerosene                   [MJ/kg]
        self.rho_ker = 0.81                                 # Density of kerosene                               [kg/L]

        '''Performance'''
        self.cruise_altitude = 11280                                # Cruise altitude                               [m]
        self.max_altitude = 12130                                   # Max certified altitude                        [m]
        self.ISA_calculator(h_input=self.cruise_altitude)

        self.rho_c = self.rho                                       # Aird density at cruise latitude           [kg/m^3]
        self.a_c = self.a
        self.T_c = self.T
        self.p_c = self.p

        self.M = 0.78
        self.V_cruise = self.M * self.a_c

        '''Aerodynamics'''
        self.e = 0.992                                              # Oswald efficiency factor
        self.C_D_0_TO_neo = 0.078                                   # Zero-lift drag coefficient of A320neo - TO
        self.C_D_0_clean_neo = 0.023                                # Zero-lift drag coefficient of A320neo - cruise
        self.C_D_0_land_neo = 0.12                                  # Zero-lift drag coefficient of A320neo - land
        self.C_D_0_climb_neo = 0.043                                # Zero-lift drag coefficient of A320neo - climb
        self.C_D_0_approach_neo = 0.053                             # Zero-lift drag coefficient of A320neo - approach
        self.C_D_0_taxi_neo = 0.055                                 # Zero-lift drag coefficient of A320neo - taxi
        self.C_L_max_clean = 1.8                                    # Max lift coefficient clean (Roskam)
        self.C_L_max_TO = 2.2                                       # Max lift coefficient during TO (Roskam)
        self.C_L_max_land = 3.0                                     # Max lift coefficient during landing (Roskam)

        self.b_in = 2 * 6.4                                         # Span of inner wing trapezoid                  [m]
        self.b_out = 2 * 10.616                                     # Span of outer wing trapezoid                  [m]
        self.b = self.b_in + self.b_out                             # Span of the entire wing                       [m]
        self.c_fus = 6.0465                                         # Chord right where the fuselage ends           [m]
        self.c_root = 7.0465                                        # Root chord including the kink                 [m]
        self.c_kink_out = 3.72                                      # Chord at the point where the kink ends        [m]
        self.c_tip = 1.488                                          # Chord at the tip of the wing excl sharklet    [m]
        self.taper_in = self.c_kink_out / self.c_root               # Taper ratio of the inner trapezoid
        self.taper_out = self.c_tip / self.c_kink_out               # Taper ratio of the outer trapezoid
        self.h_sharklet = 2.43                                      # Height of the wing sharklets                  [m]
        self.sweep_025 = 25                                         # Sweep at the quarter chord                   [deg]
        self.sweep_05 = 22.4                                        # Sweep at the half chord                      [deg]
        self.x_mac = 3.295                                          # x distance from root chord to lemac           [m]
        self.x_LEMAC = 11.88 + self.x_mac                           # nose to lemac                                 [m]
        self.mac = 4.312                                            # length of MAC                                 [m]

        self.b_h = 2 * 6.12                                         # Span of the horizontal tail                   [m]
        self.c_r_h = 3.814                                          # Root chord of the horizontal tail             [m]
        self.c_t_h = 1.186                                          # Tip chord of the horizontal tail              [m]
        self.taper_h = self.c_t_h / self.c_r_h                      # Taper ratio of the horizontal tail
        self.sweep_LE_h = 33                                        # Sweep of the LE of the horizontal tail       [deg]
        self.sweep_025_h = 29                                       # Sweep of the 25c of the horizontal tail      [deg]
        self.S_h = 31                                               # Surface area of horizontal tail              [m^2]
        self.AR_h = self.b_h**2 / self.S_h                          # Aspect ratio of horizontal tail
        self.Vh_V = 0.85                                            # Ratio between V_h and V
        self.S_elevator = 7.63                                      # Surface area of the elevators                [m^2]
        self.max_elevator_deflection_nu = 30                        # Max elevator deflection nose up              [deg]
        self.max_elevator_deflection_nd = 17                        # Max elevator deflection nose down            [deg]
        self.x_LEMAC_h = 31.605 + 1.639                             # x distane from nose to lemac of h             [m]

        self.b_v = 2 * 5.87                                         # Double the span of the evrtical tail          [m]
        self.c_r_v = 5.595                                          # Root chord of the vertical tail               [m]
        self.c_t_v = 1.822                                          # Tip chord of the vertical tail                [m]
        self.taper_v = self.c_t_v / self.c_r_v                      # Taper ratio of the vertical tail
        self.sweep_LE_v = 41                                        # Sweep of the LE of the vertical tail         [deg]
        self.S_v = 21.5                                             # Surface area of vertical tail                [m^2]
        self.AR_v = self.b_v**2 / 4 / self.S_v                      # Aspect ratio of vertical tail

        ''' Dimensions of A320-HACK'''
        self.S = 122.6                                              # Wing surface area                            [m^2]
        self.l_f = 37.57                                            # Fuselage length                               [m]
        self.height_f = 4.14                                        # Fuselage height                               [m]
        self.width_f = 3.95                                         # Fuselage width                                [m]
        self.l_cockpit = 5.04                                       # Length of the cockpit                         [m]
        self.l_cabin = 29.53 - self.l_cockpit                       # Length of the cabin + H2 tank                 [m]
        self.l_tail = self.l_f - self.l_cabin - self.l_cockpit      # Length of the tail                            [m]
        self.S_b_fus = np.pi * 0.3/2 * 0.45/2                       # Base surface area                            [m^2]
        self.sweep_LE = 27                                          # Wing sweep                                   [deg]
        self.X_root_wing = 11.88                                    # Distance from nose to root of wing            [m]
        self.X_root_vtail = 30.11                                   # Distance from nose to root of vertical tail   [m]
        self.X_root_htail = 31.60                                   # Distance from nose to root of horizontal tail [m]
        self.D_fan = 78 * 0.0254                                    # Fan diameter, 78 [in]                         [m]
        self.A_fan = np.pi * self.D_fan**2 / 4                      # Area of the fan                               [m2]
        # self.D_fan_eff = 1.787                                      # Effective fan diameter for air intake in cruise [m]
        # self.A_fan_eff = np.pi * self.D_fan_eff**2 / 4              # Effective fan area for air intake in cruise   [m2]
        self.D_h = 2.3                                              # Diameter of cowling inlet                     [m]
        self.D_n = 2.5                                              # Diameter of the cowling                       [m]
        self.D_e = 2.2                                              # Diameter of the cowling exit                  [m]
        self.l_n = 3.75                                             # Length of the fan cowling                     [m]
        self.l_eng = 5.44                                           # Length of the complete engine                 [m]
        self.y_engine = 0.34*self.b/2                               # Y location of the engine on the wing          [m]

        """Fuel constant A320-HACK"""

        self.V_H2 = V_H2/1000 * 0.885                                 # Volume required of hydrogen                [m^3]
        self.V_k = V_k/1000                                           # Volume required of kerosene                [m^3]

        self.W1_Wto = 0.990
        self.W2_W1 = 0.990
        self.W3_W2 = 0.995
        self.W4_W3 = 0.980
        self.W6_W5 = 0.990
        self.W7_W6 = 0.992
        self.c_j_kerosene = 16.68 * 10 ** (-6)                     # Specific cruise fuel consumption of neo    [kg/N*s]
        self.c_j_k_H2_cruise = 11.83 * 10 ** (-6)                  # Specific cruise fuel consumption of HACK   [kg/N*s]

        """Fuel Contants A320neo"""
        self.fuel_capacity_a320neo_FUTURE = 23859*0.9
        self.fuel_capacity_320neo = 23.859                          # Maximum Fuel capacity of A320neo          [m^3]
        self.k_d = 810.0                                            # Mass density of kerosene                  [kg/m^3]

        """ STRUCTURES"""
        self.pylon_height = 0.38     # height of the pylon of the tank [m]
        self.pylon_weight = 19.126   # weight of each pylon            [kg]
        self.x_cg_pylon = 17.14      # From nose                       [m]

        """Tank design constants"""
        self.V_centre_tank = volume_centre_tank               # Volume of each centre tank                       [m^3]
        self.V_centre_pod = volume_pod                        # Volume of each wing pod                          [m^3]
        self.V_H2_centre = volume_centre_tank * 0.885         # Volume of H2 in each centre tank                 [m^3]
        self.V_H2_pod = volume_pod * 0.885                    # Volume of H2 in each wing pod                    [m^3]

        self.center_tank_mass = mass_center_tank              # Mass of each center tank (we have 2 tanks)       [kg]
        self.pod_tank_mass = mass_pod                         # Mass of each pod tank (we have 2 tank)           [kg]
        self.centre_H2_mass = self.V_H2_centre * 71.1         # Mass of H2 in each centre tank (we have 2 tank)  [kg]
        self.pod_H2_mass = self.V_H2_pod * 71.1               # Mass of H2 in each pod tank (we have 2 tank)     [kg]

        self.x_cg_pod = 18.23                                 # From the nose                                    [m]
        self.x_cg_centertank = 20.21                          # From the nose                                    [m]
        self.y_cg_pod = 0.55*self.b/2                         # Y location of the pods on the wing                [m]

        """Weights of HACK"""
        self.Fuel_idel_taxi_take_off_HACK = 262.88                # Fuel for before take -off                       [kg]
        self.MTOW_320hack = 73500
        self.MRW_320HACK = 73900
        self.OEW_320hack = 41755
        self.payload_320hack = 13257

        """Weights of A320neo"""
        self.MTOW_320neo = 73500                                    # Maximum Take-Off weight of A320neo            [kg]
        self.MLW_320neo = 66300                                     # Maximum Landing weight of A320neo             [kg]
        self.MRW_320neo = 73900                                     # Maximum Ramp weight of A320neo                [kg]
        self.MZFW_320neo = 62800                                    # Maximum Zero fuel weight of A320neo           [kg]
        self.MPLW_320neo = 18240                                    # Maximum Payload weight of A320neo             [kg]
        self.OEW_320neo = 44560                                     # Operational Empty weight of A320neo           [kg]
        self.Wing_Weight_320neo = 9150                              # Wing weight DOI:10.5139/IJASS.2014.15.4.38    [kg]
        self.W_engine = 2990                                        # Weight of one engine                          [kg]
        self.Fuel_idel_taxi_take_off_320neo = 400                   # Fuel for before take -off                     [kg]
        self.Max_fuel_mass_capacity_320neo = self.fuel_capacity_320neo * self.k_d # Maximum kerosene mass of A320neo[kg]
        self.x_cg_320neo_zf = 0.29
        self.x_cg_320neo_mtow = 0.275
        # self.x_cg_hack = self.x_cg_320neo_zf * self.MZFW_320neo + \
        #                 (self.x_cg_320neo_mtow -self.x_cg_320neo_zf)*self.Max_fuel_mass_capacity_320neo

        """Weights of A321neo"""
        self.MTOW_321neo = 89000                                    # Maximum Take-Off weight of A321neo            [kg]
        self.MLW_321neo = 77300                                     # Maximum Landing weight of A321neo             [kg]
        self.MRW_321neo = 89400                                     # Maximum Ramp weight of A321neo                [kg]
        self.MZFW_321neo = 73300                                    # Maximu Zero fuel weight of A321neo            [kg]
        self.MPLW_321neo = 22910                                    # Maximum Payload weight of A321neo             [kg]
        self.OEW_321neo = self.MZFW_321neo - self.MPLW_321neo       # Operational Empty weight of A320neo           [kg]

        """Dimensions of A320neo and A321neo"""
        self.l_f_321neo = 44.51                                     # Fuselage length of A321neo                    [m]
        self.l_f_320neo = 37.57                                     # Fuselage length of A320neo                    [m]
        self.l_cockpit_320neo = 5.04                                # Length of the cockpit of A320neo              [m]
        self.l_cabin_320neo = 29.53 - self.l_cockpit_320neo         # Length of the cabin of A320neo                [m]
        self.l_tail_320neo = self.l_f_320neo - 29.53                # Length of the tail of A320neo                 [m]

        """Propulsion"""
        self.cp_air = 1000                                          # Specific heat constant air                    [J/kg/K]
        self.cp_gas = 1150                                          # Specific heat constant gas                    [J/kg/K]
        self.k_air = 1.4                                            # Ration of specific heat for air
        self.k_gas = 1.3379776344421168                             # Ration of specific heat for air

        # self.N2_cp_data = np.array(np.genfromtxt('N2_cp.dat'))      # cp vs. T data for N2          T[K]; cp[kJ/(kg*K)]
        self.molarmass_N2 = 28.01340                                # Molar mass of N2                          [g/mol]

        # self.h2_cp_data = np.array(np.genfromtxt('h2_cp.dat'))      # cp vs. T data for h2          T[K]; cp[kJ/(kg*K)]
        self.molarmass_h2 = 2.01588                                 # Molar mass of h2                          [g/mol]

        # self.C12H26_cp_data = np.array(np.genfromtxt('C12H26_cp.dat'))  # cp vs. T data for dodecane                T[K]; cp[J/(mol*K)]
        self.h0_C12H26 = -290.90                                        # Zero enthalpy of dodecane                 [kJ/mol]         # https://www.chemeo.com/cid/34-125-5/n-Dodecane
        self.molarmass_C12H26 = 170.3348                                # Molar mass of dodecane                    [g/mol]

        self.stoich_ratio_ker = 1/15.66 #FAR
        self.stoich_ratio_h2 = 1/34.3 #FAR

        """" Altitude and speed --- USE DATAFRAME """
        # self.phases = np.array(['idle', 'taxi_out', 'takeoff', 'climb', 'cruise', 'approach', 'taxi_in'])
        # self.Thrust = np.array( [0.0993, 0.124, 0.45, 0.77, 0.29, 0.28, 0.068] ) * 35000*4.44822162
        # self.M0 = np.array([0.0, 0.01, 0, 0.5, 0.78, 0.5, 0.01])  # [-] Mach number
        # self.h = np.array([0, 0, 0, 3000, 11280, 300, 0])  # [m] altitude
        # self.T0, self.p0, self.rho0, self.a0 = np.array([]), np.array([]), np.array([]), np.array([])


        # for i in self.h:
        #     self.ISA_calculator(h_input=i)
        #     self.T0 = np.append(self.T0, self.T)
        #     self.p0 = np.append(self.p0, self.p)
        #     self.rho0 = np.append(self.rho0, self.rho)
        #     self.a0 = np.append(self.a0, self.a)
        #
        # self.v0 = self.M0 * self.a0

    def chord(self, x):
        if x <= 0.5*self.b_in:
            c = (self.c_kink_out - self.c_root) / (0.5 * self.b_in) * x + self.c_root
        else:
            c = (self.c_tip - self.c_kink_out) / (0.5 * self.b_out) * (x - 0.5 * self.b_in) + self.c_kink_out
        return c

    # def engine_data_neo(self, phase):
    #     if phase in ['climb','cruise','approach']:
    #         self.eta_inlet = 0.9208608681597723  # calculated
    #         self.PR_fan = 1.4206
    #         self.eta_fan = 0.90445
    #         self.BPR = 11.24426
    #         self.eta_LPC = 0.90019
    #         self.eta_HPC = 0.95469  # calculated # 0.91449 given
    #         self.PR_LPC = 2.69419
    #         self.PR_HPC = 9.73784
    #         self.eta_mech = 0.6 # 0.9  # not used
    #         self.eta_cc = 0.995  # that of Leap-1B, assumed
    #         self.PR_cc = 0.9395309126896629
    #         self.T04 = 1459.30433  # [K]
    #         self.eta_LPT = 0.9405
    #         self.eta_HPT = 0.9328  # computed # 0.91898 #(given)
    #         self.PR_LPT = 7.9204
    #         self.PR_HPT = 3.81933
    #         self.eta_nozzle = 0.981797  # 1.0737340755627587 (computed) # previous assumption: 0.98
    #         self.PR_noz_core = 0.985443  # Between stations 5 and 7
    #         self.PR_cr_noz_core = 1.653828  # computed
    #         self.PR_noz_fan = 0.987444  # Between stations 21 and 16
    #
    #         self.mf_air_init = self.rho0 * self.A_fan_eff * self.v0
    #
    #     elif phase in ['idle', 'taxi_out','takeoff', 'taxi_in']:
    #         self.eta_inlet = 0.92  # assumed
    #         self.PR_inlet = 0.9699975  # calculated
    #         self.PR_fan = 1.4
    #         self.eta_fan = 0.93
    #         self.BPR = 11.1
    #         self.eta_LPC = 0.92
    #         self.eta_HPC = 0.97735  # calculated # 0.92 given
    #         self.PR_LPC = 1.99241
    #         self.PR_HPC = 11.92998521
    #         self.eta_mech = 0.78  # not used
    #         self.eta_cc = 0.995  # that of Leap-1B, assumed
    #         self.PR_cc = 0.94
    #         self.T04 = 1630.39416  # [K]
    #         self.eta_LPT = 0.94
    #         self.eta_HPT = 0.94  # computed # 0.91898 #(given)
    #         self.PR_LPT = 6.36217496  # 4.835766
    #         self.PR_HPT = 3.82808
    #         self.eta_nozzle = 0.98  # assumed
    #         self.PR_noz_core = 0.99  # Between stations 5 and 7
    #         self.PR_cr_noz_core = 1.2380755  # computed
    #         self.PR_noz_fan = 0.99  # Between stations 21 and 16
    #         self.mf_air_init = 510  # [kg/s] from Arvind's document
    #
    #     self.mr_h2 = np.array([0, 0, 0, 0, 0, 0, 0])
    #     self.mr_ker = 1 - self.mr_h2
    #     self.ER_h2 = np.array([0, 0, 0, 0, 0, 0, 0])
    #     self.ER_ker = 1 - self.mr_h2
    #     self.LHV_f = np.array([self.LHV_ker]*7) # [MJ/kg]
    #
    #     self.ratio_air_cc = np.array(np.genfromtxt('mr_cc_neo.dat')) # percentage of core air that is used in combustion                                                                                   # [kg/s]
    #
    # def engine_data_hack(self, phase):
    #     if phase in ['climb','cruise','approach']:
    #         self.eta_inlet = 0.9208608681597723  # calculated
    #         self.PR_fan = 1.4206
    #         self.eta_fan = 0.90445
    #         self.BPR = 11.24426
    #         self.eta_LPC = 0.90019
    #         self.eta_HPC = 0.95469  # calculated # 0.91449 given
    #         self.PR_LPC = 2.69419
    #         self.PR_HPC = 9.73784
    #         # self.eta_mech_H =  0.644335665181638
    #         # self.eta_mech_L = 1
    #         self.eta_mech = 0.9  # not used
    #         self.eta_cc = 0.995  # that of Leap-1B, assumed
    #         self.PR_cc = 0.9395309126896629
    #         self.T04 = 1459.30433  # [K]
    #         self.eta_LPT = 0.9405
    #         self.eta_HPT = 0.9328  # computed # 0.91898 #(given)
    #         self.PR_LPT = 7.9204
    #         self.PR_HPT = 3.81933
    #         self.eta_nozzle = 0.981797  # 1.0737340755627587 (computed) # previous assumption: 0.98
    #         self.PR_noz_core = 0.985443  # Between stations 5 and 7
    #         self.PR_cr_noz_core = 1.653828  # computed
    #         self.PR_noz_fan = 0.987444  # Between stations 21 and 16
    #         self.mf_air_init = self.rho0[4] * self.A_fan_eff * self.v0[4]
    #
    #     elif phase in ['idle', 'taxi out', 'takeoff', 'taxi in']:
    #         self.eta_inlet = 0.92  # assumed
    #         self.PR_inlet = 0.9699975  # calculated
    #         self.PR_fan = 1.4
    #         self.eta_fan = 0.93
    #         self.BPR = 11.1
    #         self.eta_LPC = 0.92
    #         self.eta_HPC = 0.97735  # calculated # 0.92 given
    #         self.PR_LPC = 1.99241
    #         self.PR_HPC = 11.92998521
    #         # self.eta_mech_H =  0.644335665181638
    #         # self.eta_mech_L = 1
    #         self.eta_mech = 0.9  # not used
    #         self.eta_cc = 0.995  # that of Leap-1B, assumed
    #         self.PR_cc = 0.94
    #         self.T04 = 1630.39416  # [K]
    #         self.eta_LPT = 0.94
    #         self.eta_HPT = 0.94  # computed # 0.91898 #(given)
    #         self.PR_LPT = 6.36217496  # 4.835766
    #         self.PR_HPT = 3.82808
    #         self.eta_nozzle = 0.98  # assumed
    #         self.PR_noz_core = 0.99  # Between stations 5 and 7
    #         self.PR_cr_noz_core = 1.2380755  # computed
    #         self.PR_noz_fan = 0.99  # Between stations 21 and 16
    #         self.mf_air_init = 510  # [kg/s]
    #
    #     # Fuel properties
    #     self.mr_h2 = MR_h2
    #     self.mr_ker = MR_ker
    #     self.ER_h2 = ER_h2
    #     self.ER_ker = ER_ker
    #     self.LHV_f = LHV_hack # [MJ/kg]
    #     self.ratio_air_cc = np.array(np.genfromtxt('mr_cc_hack.dat'))



    # def fuselage_length(self,vol_eff, vol_fus):
    #     """
    #
    #     :param vol_eff: Volumetric efficiency of integral tanks (Ratio of usable tank volume-to-volume occupied
    #                     in the fuselage                                         [-]
    #     :param vol_fus: Volume available for tanks in the center wingbox        [m^3]
    #     :return:
    #     """
    #     self.V_H2_center_w = vol_eff * vol_fus  # Volume of hydrogen stored on center wing box [m^3]
    #     self.V_H2_ext_fus = self.V_H2 - self.V_H2_center_w  # Volume of hydrogen stored on extended fuselage[m^3]
    #     self.l_f = self.V_H2_ext_fus/(pi * self.width_f/2 * self.height_f/2)

    def speed_of_sound(self, T):
        """
        Compute the speed of sound in air for any temperature
        :param T: Temperature [K]
        :return: the speed of sound in air [m/s]
        """
        a = np.sqrt(self.gamma * self.R * T)
        return a

    def ISA_temp(self, gr, T_init, h_init, h):
        """
        International Standard Atmosphere equation for temperature
        :param gr: Temperature gradient [K/m]
        :param T_init: Temperature at the beginning of the atmospheric layer [K]
        :param h_init: Altitude at the beginning of the atmospheric layer [m]
        :param h: Altitude of interest [m]
        :return: Temperature at altitude h according to the ISA [K]
        """
        ISA_T = T_init + gr * (h - h_init)
        return ISA_T

    def ISA_grad_pressure(self, p1_grad, T1_grad, ISA_T, gr):
        """
        International Standard Atmosphere equation for pressure in a layer with changing temperature
        :param p1_grad: Pressure at the beginning of the atmospheric layer [Pa]
        :param T1_grad: Temperature at the beginning of the atmospheric layer [K]
        :param ISA_T: Temperature at the desired altitude [K]
        :param gr: Temperature gradient [K/m]
        :return: Pressure at altitude h according to the ISA [Pa]
        """
        ISA_p = p1_grad * (ISA_T / T1_grad)**(-self.g_0 / (gr * self.R))
        return ISA_p

    def ISA_iso_pressure(self, p1_iso, T1_iso, h1_iso, h):
        """
        International Standard Atmosphere equation for pressure in a layer with constant temperature
        :param p1_iso: Pressure at the beginning of the atmospheric layer [Pa]
        :param T1_iso: Temperature at the beginning of the atmospheric layer [K]
        :param h1_iso: Altitude at the beginning of the atmospheric layer [m]
        :param h: Altitude of interest [m]
        :return: Pressure at altitude h according to the ISA [Pa]
        """
        ISA_p = p1_iso * np.exp(-self.g_0 / (self.R * T1_iso) * (h - h1_iso))
        return ISA_p

    def ISA_density(self, ISA_p, ISA_T):
        """
        Ideal gas law for air
        :param ISA_p: Pressure at altitude h according to the ISA [Pa]
        :param ISA_T: Temperature at altitude h according to the ISA [K]
        :return: Air density at altitude h according to the ISA [kg/m^3]
        """
        rho = ISA_p / (ISA_T * self.R)
        return rho

    def ISA_calculator(self, h_input):
        """
        Compute the atmospheric conditions at altitude h_input
        :param h_input: Altitude of interest [m]
        :return: Nothing, it only defines attributes of the class
        """

        h_atmos = np.array([0., 11000., 20000., 32000., 47000., 51000., 71000., 86000])
        p_atmos = np.array([101325., 22625.79, 5471.9, 867.25, 110.766, 66.848, 3.94898, 0.301618])
        T_atmos = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 184.65])
        grad_atmos = np.array([-0.0065, 0., 0.001, 0.0028, 0., -0.0028, -0.002, 0.])

        loc_atmos = np.where(h_atmos <= h_input)[0][-1]
        h1 = h_atmos[loc_atmos]
        p1 = p_atmos[loc_atmos]
        T1 = T_atmos[loc_atmos]
        grad = grad_atmos[loc_atmos]

        self.T = self.ISA_temp(grad, T1, h1, h_input)
        self.a = self.speed_of_sound(self.T)

        if grad != 0.:
            self.p = self.ISA_grad_pressure(p1, T1, self.T, grad)

        else:
            self.p = self.ISA_iso_pressure(p1, T1, h1, h_input)

        self.rho = self.ISA_density(self.p, self.T)


# class engine_data_neo:
#     def __init__(self):
#         self.eta_inlet = 0.97
#         self.PR_fan = 1.6
#         self.eta_fan = 0.93
#         self.BR = 11.1
#         self.eta_LPC = 0.92
#         self.eta_HPC = 0.92
#         self.eta_LPT = 0.94
#         self.eta_HPT = 0.94
#         self.eta_mech = 0.9
#         self.eta_cc = 0.99
#         self.PR_LPC = 2
#         self.PR_HPC = 11.93
#         self.eta_nozzle = 0.98
#         self.PR_cc = 0.96
#         self.T04 = 1630 # [K]
#         self.LHV_f = 43.2 # [MJ/kg]


# class engine_data_hack:
#     def __init__(self):
        # self.eta_inlet = 0.97
        # self.PR_fan = 1
        # self.eta_fan = 0.93
        # self.BR = 12
        # self.eta_LPC = 0.92
        # self.eta_HPC = 0.92
        # self.eta_LPT = 0.94
        # self.eta_HPT = 0.94
        # self.eta_mech = 0.9
        # self.eta_cc = 0.99
        # self.PR_LPC = 2.3
        # self.PR_HPC = 13
        # self.eta_nozzle = 0.98
        # self.PR_cc = 0.96
        # self.T04 = 1630 # [K]
        #
        # # Fuel properties
        #
        # self.mr_h2 = np.array([ 1, 1, 0.1376, 0.1376, 0.1376, 0.1376, 1  ])
        # self.mr_ker = 1 - self.mr_h2
        #
        # self.ER_h2 = ( self.mr_h2*self.LHV_h2 ) / (  self.mr_h2*self.LHV_h2 + self.mr_ker*self.LHV_ker)
        # self.ER_ker = ( self.mr_ker*self.LHV_ker ) / (  self.mr_h2*self.LHV_h2 + self.mr_ker*self.LHV_ker)
        #
        # # find LHV_f for each phase, according to mass fractions
        # self.LHV_f = self.ER_h2*self.LHV_h2 + self.ER_ker*self.LHV_ker  # [MJ/kg]

# Try out the class

#Weigth_Centre_Tanks

if __name__ == '__main__':
    c = Constants()
    c.ISA_calculator(h_input=53000)
    print('\n T = ', c.T, ' K',
          '\n P = ', c.p, ' Pa',
          '\n rho = ', c.rho, ' kg/m^3')

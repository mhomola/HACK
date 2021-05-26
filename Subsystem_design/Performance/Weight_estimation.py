#Following file computes a first estimation of the weight of the A320-HACK

from tank_sizing import liquid_H_tanks
#from Subsystem_design.fuel_required import fuel_volume_calc
from Subsystem_design.common_constants import Constants
from fuel_constants import *
from math import pi,log
import matplotlib.pyplot as plt
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
import numpy as np



class Compute_weight(Constants):

    def __init__(self):
        super().__init__()
        self.m_k = 12883.86                                        # Kerosene mass                                  [kg]
        self.sigma_y_FS = 350 * 10 ** 6                            # Yield strendth of Steel (assumption)           [Pa]
        self.rho_FS = 7870                                         # Density of steel                           [kg/m^3]
        self.D_FS = 0.4                                            # Diameter of the pipelines                      [m]
        self.p_FS = 1.66 * p_tank                                  # Assumed pressure inside the tanks              [Pa]
        self.FS_length = 66.5                                      # Length of feeding system(roughly)              [m]
        self.t_FS = (self.p_FS * self.D_FS / 2) / (self.sigma_y_FS)  # Thickness of the feeding system                [m]




    def Tank_mass(self,h2_vol_center, h2_vol_f):
        """

        :param h2_vol_center: Volume of H2 contained in the center tanks                    [m^3]
        :param h2_vol_f: Volume of H2 contained in the extended fuselage tanks              [m^3]
        :return: m_H2_center_w, m_tank_center_w, m_H2_ext_f, m_tank_ext_f
        """
        self.m_H2_center_w,_, _ = liquid_H_tanks(h2_vol_center*1000)                    # Hydrogen mass             [kg]
        self.m_tank_center_w = self.center_tank_mass                                    # Center Tank's mass       [kg]
        self.m_H2_ext_f,_, _ = liquid_H_tanks(h2_vol_f * 1000)                          # Hydrogen mass             [kg]
        self.m_tank_ext_f = self.fuselage_tank_mass                                     # Fuselage tank's mass      [kg]

        return self.m_H2_center_w, self.m_tank_center_w,self.m_H2_ext_f, self.m_tank_ext_f

    def Feeding_sys_m(self):

        self.FS_mass = pi * self.D_FS * self.t_FS * self.FS_length * self.rho_FS        # Mass of feeding system     [kg]

        return self.FS_mass

    def Struc_m(self):
        # Added mass per meter of added fuselage                    [kg/m]
        extd_m_r = (self.OEW_321neo - self.OEW_320neo)/(self.l_f_321neo - self.l_f_320neo)
        self.struc_mass = extd_m_r * (self.l_f - self.l_f_320neo)      # Extra Structural mass               [kg]

        return self.struc_mass

    def weight_break_down_HACK(self,h2_vol_center, h2_vol_f):
        """

        :param h2_vol_center: Volume of H2 contained in the center tanks                    [m^3]
        :param h2_vol_f: Volume of H2 contained in the extended fuselage tanks              [m^3]
        :return: None
        """
        #Running Tank mass, feeding system mass and struc mass  function
        self.Tank_mass(h2_vol_center, h2_vol_f)
        self.Feeding_sys_m()
        self.Struc_m()

        # Operational Empty Weight of Hack    [kg]
        self.OEW_HACK = self.OEW_320neo + self.m_tank_center_w + self.m_tank_ext_f + self.FS_mass + self.struc_mass

        self.Max_fuel_mass_capacity_HACK = self.m_H2_center_w + self.m_H2_ext_f + self.m_k

        self.MPLW_HACK = self.MPLW_320neo                         # Maximum Payload weight of HACK                  [kg]
        self.MZFW_HACK = self.MPLW_HACK + self.OEW_HACK           # Maximum zero fuel weight of HACK                [kg]
        self.MTOW_HACK = (self.OEW_HACK - 16914.5)/0.376          # Maximum take-off weight of Hack                 [kg]

        self.Max_fuel_at_max_PL_HACK = self.MTOW_HACK - self.MPLW_HACK - self.OEW_HACK   #Maximum fuel at maximum payload[kg]

        self.MRW_HACK = self.OEW_HACK + self.MPLW_HACK + self.Max_fuel_at_max_PL_HACK + self.Fuel_idel_taxi_take_off_HACK


class performance(Compute_weight):
    #super().__init__()

    def flight_profile_weights(self,Mf,Mto):
        #Starts with MTOW, obtain Wf/Wto:
        Wf_Wto = Mf/Mto                                         # Fuel fraction between maximum fuel mass          [-]
                                                                # @ max. payload, and MTOW
        M_ff = 1 - Wf_Wto
        W5_W4 = M_ff  * (1/self.W4_W3) * (1/self.W6_W5) * (1/self.W7_W6)   # Fuel ratio between the end            [-]
                                                                           # and beginning of cruise
        return W5_W4

    def Range(self,L_D_ratio,cruise_f_ratio,SFC):
        R = (self.V_cruise/(SFC*self.g_0)) * (L_D_ratio) *log(1/cruise_f_ratio)
        return R

    def payload_range_dia_HACK(self,L_over_D,h2_vol_center, h2_vol_f,SFC):

        # Function to plot the payload-range diagram
        #Payload as a funstion of range
        #Get Range equation, get Wf and Wi from fuel fractions
        # R = V/cp * L/D * Wi/Wf
        # Running compute_weight class to get variables inside functions
        self.weight_break_down_HACK(h2_vol_center, h2_vol_f)


        #Point A: R = 0 and W_Pl = MPLW
        R_A = 0
        W_Pl_A = self.MPLW_HACK                                               # Payload mass at point A             [kg]

        #Point B: R = ?? and W_Pl = MPLW
        Mto = self.MTOW_HACK
        W_Pl_B = self.MPLW_HACK                                               # Payload mass at point B             [kg]
        Mf = self.Max_fuel_at_max_PL_HACK                                     # Maximum fuel at maximum payload    [kg]
        W5W4_B = self.flight_profile_weights(Mf,Mto)                          # Fuel ratio during cruise for point B [-]
        R_B = self.Range(L_over_D,W5W4_B,SFC)                                     # Range at point B                     [m]


        # Point C: R = ?? and Wf = MFW, Wpl = MTOW - OEM - Maxfuel
        Mf = self.Max_fuel_mass_capacity_HACK - self.Fuel_idel_taxi_take_off_HACK  # Maximum fuel capacity of HACK      [kg]
        W_Pl_C = Mto - Mf - self.OEW_HACK                                     #Payload mass at point C             [kg]
        W5W4_C = self.flight_profile_weights(Mf, Mto)                         # Fuel ratio during cruise for point C [-]
        R_C = self.Range(L_over_D, W5W4_C,SFC)                                    # Range at point C                     [m]


        # Point D: R = ?? and W_Pl = 0--> W_fuel = Maxfuel
        Mf = self.Max_fuel_mass_capacity_HACK - self.Fuel_idel_taxi_take_off_HACK  # Maximum fuel capacity of HACK      [kg]
        W_PL_D = 0                                                            # Payload mass at point D            [kg]
        Mto = self.OEW_HACK + Mf + W_PL_D                                     # Take off weight at point D         [kg]
        W5W4_D = self.flight_profile_weights(Mf, Mto)                         # Fuel ratio during cruise for point D [-]
        R_D = self.Range(L_over_D, W5W4_D,SFC)                                    # Range at point D                     [m]

        Range_array = np.array([R_A*0.001,R_B*0.001,R_C*0.001,R_D*0.001])
        Payload_array = np.array([W_Pl_A,W_Pl_B,W_Pl_C,W_PL_D])
        plt.plot(Range_array,Payload_array,marker = '*',color = 'tab:red')
        plt.xlabel('Range [km]')
        plt.ylabel('Payload Mass [kg]')
        plt.show()

    def payload_range_dia_320neo(self,L_over_D,SFC):

        # Function to plot the payload-range diagram
        #Payload as a funstion of range
        #Get Range equation, get Wf and Wi from fuel fractions
        # R = V/cp * L/D * Wi/Wf
        # Running compute_weight class to get variables inside functions


        #Point A: R = 0 and W_Pl = MPLW
        R_A = 0
        W_Pl_A = self.MPLW_320neo                                               # Payload mass at point A             [kg]

        #Point B: R = ?? and W_Pl = MPLW
        Mto = self.MTOW_320neo
        W_Pl_B = self.MPLW_320neo                                               # Payload mass at point B             [kg]
        Mf = self.MTOW_320neo - self.MPLW_320neo - self.OEW_320neo              # Maximum fuel at maximum payload    [kg]
        W5W4_B = self.flight_profile_weights(Mf,Mto)                            # Fuel ratio during cruise for point B [-]
        R_B = self.Range(L_over_D,W5W4_B,SFC)                                   # Range at point B                     [m]


        # Point C: R = ?? and Wf = MFW, Wpl = MTOW - OEM - Maxfuel
        Mf = self.Max_fuel_mass_capacity_320neo - self.Fuel_idel_taxi_take_off_320neo # Maximum fuel capacity of A320neo     [kg]
        W_Pl_C = Mto - Mf - self.OEW_320neo                                  # Payload mass at point C              [kg]
        W5W4_C = self.flight_profile_weights(Mf, Mto)                        # Fuel ratio during cruise for point C [-]
        R_C = self.Range(L_over_D, W5W4_C,SFC)                               # Range at point C                     [m]


        # Point D: R = ?? and W_Pl = 0--> W_fuel = Maxfuel
        Mf = self.Max_fuel_mass_capacity_320neo - self.Fuel_idel_taxi_take_off_320neo # Maximum fuel capacity of A320neo      [kg]
        W_PL_D = 0                                                            # Payload mass at point D            [kg]
        Mto = self.OEW_320neo + Mf + W_PL_D                                   # Take off weight at point D         [kg]
        W5W4_D = self.flight_profile_weights(Mf, Mto)                         # Fuel ratio during cruise for point D [-]
        R_D = self.Range(L_over_D, W5W4_D,SFC)                                # Range at point D                     [m]

        Range_array = np.array([R_A*0.001,R_B*0.001,R_C*0.001,R_D*0.001])
        Payload_array = np.array([W_Pl_A,W_Pl_B,W_Pl_C,W_PL_D])
        plt.plot(Range_array,Payload_array,marker = '*',color = 'tab:red')
        plt.xlabel('Range [km]')
        plt.ylabel('Payload Mass [kg]')
        plt.show()

if __name__ == '__main__':
    const = Constants()
    """Weights of A320-HACK"""

    AC_weights = Compute_weight()                                               # Initiallize class of weight estimation
    AC_weights.weight_break_down_HACK(h2_vol_center=11.34,h2_vol_f=30.69)       # Based on H2 volume estimates
    AC_weights.Struc_m()

    Aerodynamic_charac = AerodynamicCharacteristics()
    Aerodynamic_charac.L_over_D_cruise()

    Performance = performance()
    Performance.payload_range_dia_HACK(L_over_D=Aerodynamic_charac.L_D_ratio_HACK,h2_vol_center=11.34,h2_vol_f=30.69,SFC= const.c_j_k_H2_cruise)
    Performance.payload_range_dia_320neo(L_over_D=Aerodynamic_charac.L_D_ratio_HACK,SFC = const.c_j_kerosene)

    print('The OEW of the A320HACK is:',AC_weights.OEW_HACK)
    print('The MPLW of the A320HACK is:', AC_weights.MPLW_HACK)
    print('The max fuel mass of the A320HACK is:', AC_weights.Max_fuel_mass_capacity_HACK)
    print('The max fuel @ MPLW of the A320HACK is:', AC_weights.Max_fuel_at_max_PL_HACK)
    print('The MZFW of the A320HACK is:', AC_weights.MZFW_HACK)
    print('The MTOW of the A320HACK is:', AC_weights.MTOW_HACK)

    OEW = np.array([const.OEW_320neo, const.OEW_321neo])
    MTOW = np.array([const.MTOW_320neo, const.MTOW_321neo])
    plt.plot(MTOW,OEW,marker = '*')
    plt.ylabel('OEW [kg]')
    plt.xlabel('MTOW [kg]')
    plt.show()

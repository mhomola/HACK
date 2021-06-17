from Subsystem_design.common_constants import Constants
import matplotlib.pyplot as plt
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
import numpy as np
from math import log

class Compute_weight(Constants):
    def __init__(self):
        super().__init__()


    def Tank_mass(self):
        self.tank_mass = 2*self.pod_tank_mass
    def Feeding_sys_m(self):
        self.Feeding_mass = 22                          # Feeding system from tanks to engine + tank to tank
    def Struc_m(self):
        self.struc_mass = self.Wing_Weight_320HACK-self.Wing_Weight_320neo + self.pylon_weight
    def Fuel_cell_m(self):
        self.FC_m = 781                                 # Entire system [kg]
    def weight_break_down_HACK(self):
        self.Tank_mass()
        self.Feeding_sys_m()
        self.Struc_m()
        self.Fuel_cell_m()

        self.OEW_HACK = self.OEW_320neo +self.struc_mass +self.Feeding_mass + self.FC_m + self.tank_mass

        self.Max_fuel_mass_capacity_HACK = self.W_kerosene + 2* self.pod_H2_mass

        self.MPLW_HACK = self.MPLW_320neo                           # Maximum Payload weight of HACK                  [kg]
        self.MZFW_HACK = self.MPLW_HACK + self.OEW_HACK             # Maximum zero fuel weight of HACK                [kg]
        self.MTOW_HACK = self.MTOW_320neo                           # Maximum take-off weight of Hack                 [kg]

        self.Max_fuel_at_max_PL_HACK = self.MTOW_HACK - self.MPLW_HACK - self.OEW_HACK  # Maximum fuel at maximum payload[kg]

        self.MRW_HACK = self.MRW_320HACK

class performance(Compute_weight):
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
    def payload_range_dia_HACK(self,L_over_D,SFC):

        # Function to plot the payload-range diagram
        #Payload as a funstion of range
        #Get Range equation, get Wf and Wi from fuel fractions
        # R = V/cp * L/D * Wi/Wf
        # Running compute_weight class to get variables inside functions
        self.weight_break_down_HACK()


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

if __name__ == '__main__':
    const = Constants()
    """Weights of A320-HACK"""

    AC_weights = Compute_weight()                                               # Initiallize class of weight estimation
    AC_weights.weight_break_down_HACK()
    Aerodynamic_charac = AerodynamicCharacteristics()
    Aerodynamic_charac.L_over_D_cruise()

    Performance = performance()
    Performance.payload_range_dia_HACK(L_over_D=Aerodynamic_charac.L_D_ratio_HACK,SFC= const.c_j_k_H2_cruise)
    #Performance.payload_range_dia_320neo(L_over_D=Aerodynamic_charac.L_D_ratio_HACK,SFC = const.c_j_kerosene)
    #W4_W5 = (3200*10**3)/(((const.M * Aerodynamic_charac.a)/(const.c_j_k_H2_cruise*const.g_0)) * Aerodynamic_charac.L_D_ratio_HACK)

    print('The OEW of the A320HACK is:',AC_weights.OEW_HACK)
    print('The MPLW of the A320HACK is:', AC_weights.MPLW_HACK)
    print('The max fuel mass of the A320HACK is:', AC_weights.Max_fuel_mass_capacity_HACK)
    print('The max fuel @ MPLW of the A320HACK is:', AC_weights.Max_fuel_at_max_PL_HACK)
    print('The MZFW of the A320HACK is:', AC_weights.MZFW_HACK)
    print('The MTOW of the A320HACK is:', AC_weights.MTOW_HACK)
    print('The extra structural mass of the A320HACK is:', AC_weights.struc_mass)
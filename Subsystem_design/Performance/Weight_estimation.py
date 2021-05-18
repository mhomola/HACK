#Following file computes a first estimation of the weight of the A320-HACK

from tank_sizing import liquid_H_tanks
from Subsystem_design.fuel_required import fuel_volume_calc
from Subsystem_design.common_constants import *
from Subsystem_design.common_constants import Constants
from fuel_constants import *
from math import pi,log
import matplotlib.pyplot as plt



class compute_weight(Constants):

    def __init__(self):
        super().__init__()
        self.m_k = 11595.96                                        # Kerosene mass                                  [kg]
        self.sigma_y_FS = 350 * 10 ** 6                            # Yield strendth of Steel (assumption)           [Pa]
        self.rho_FS = 7750                                         # Density of steel                           [kg/m^3]
        self.D_FS = 0.4                                            # Diameter of the pipelines                      [m]
        self.p_FS = 1.66 * p_tank                                  # Assumed pressure inside the tanks              [Pa]
        self.FS_length = 66.5                                      # Length of feeding system(roughly)              [m]
        self.t_FS = (self.p_FS * self.D_FS / 2) / (self.sigma_y_FS)  # Thickness of the feeding system                [m]




    def Tank_mass(self,h2_vol_center, h2_vol_f):

        self.m_H2_center_w, self.m_tank_center_w, _ = liquid_H_tanks(h2_vol_center*1000) # Hydrogen mass and Tank mass[kg]
        self.m_H2_ext_f, self.m_tank_ext_f, _ = liquid_H_tanks(h2_vol_f * 1000)          # Hydrogen mass and Tank mass[kg]

    def Feeding_sys_m(self):

        self.FS_mass = pi * self.D_FS * self.t_FS * self.FS_length * self.rho_FS        # Mass of feeding system     [kg]

        return self.FS_mass

    def Struc_m(self):
        # Added mass per meter of added fuselage                    [kg/m]
        extd_m_r = (self.OEW_321neo - self.OEW_320neo)/(self.l_f_321neo - self.l_f_320neo)
        self.struc_mass = extd_m_r * (self.l_f - self.l_f_320neo)      # Extra Structural mass               [kg]

        return self.struc_mass

    def weight_break_down_HACK(self):

        # Operational Empty Weight of Hack    [kg]
        self.OEW_HACK = self.OEW_320neo + self.m_tank_center_w + self.m_tank_ext_f + self.FS_mass + self.struc_mass

        self.Max_fuel_mass = self.m_H2_center_w + self.m_H2_ext_f + self.m_k

        self.MPLW_HACK = self.MPLW_320neo                       # Maximum Payload weight of HACK                  [kg]
        self.MZFW_HACK = self.MPLW_HACK + self.OEW_HACK         # Maximum zero fuel weight of HACK                [kg]
        self.MTOW_HACK = self.MTOW_320neo                       # Maximum take-off weight of Hack                 [kg]
        self.Max_fuel_at_max_PL = self.MTOW_HACK - self.MPLW_HACK - self.OEW_HACK   #Maximum fuel at maximum payload[kg]

class performance(Constants,compute_weight):
    def __init__(self):
        super().__init__()
        self.W1_Wto = 0.990
        self.W2_W1 = 0.990
        self.W3_W2 = 0.995
        self.W4_W3 = 0.980
        self.W6_W5 = 0.990
        self.W7_W6 = 0.992

    def flight_profile_weights(self,Mf,Mto):
        #Starts with MTOW, obtain Wf/Wto:
        Wf_Wto = Mf/Mto                                         # Fuel fraction between maximum fuel mass          [-]
                                                                # @ max. payload, and MTOW
        M_ff = 1 - Wf_Wto
        W5_W4 = M_ff * (1/self.W1_Wto) * (1/self.W2_W1) * (1/self.W3_W2) \
                * (1/self.W4_W3) * (1/self.W6_W5) * (1/self.W7_W6)         # Fuel ratio between the end            [-]
                                                                           # and beginning of cruise
        return W5_W4

    def Range(self,L,D,cruise_f_ratio):
        R = (self.V_cruise/(self.c_j_kerosene*self.g_0)) * (L/D) *log(1/cruise_f_ratio)

    def payload_range_dia(self,L,D):

        # Function to plot the payload-range diagram
        #Payload as a funstion of range
        #Get Range equation, get Wf and Wi from fuel fractions
        # R = V/cp * L/D * Wi/Wf
        #PointA: R = 0 and W_Pl = MPLW
        #PointB: R = ?? and W_Pl = MPLW--> From MTOW calculate maxfuel for MPLW
        #PointC: R = ?? and Wf = MFW, Wpl = MTOW - OEM - Maxfuel
        #Point D: R = ?? and W_Pl = 0--> W_fuel = Maxfuel

        "Point A: R = 0 and W_Pl = MPLW"
        R_A = 0
        W_Pl_A = self.MPLW_HACK                                               # Payload mass at point A             [kg]

        #Point B: R = ?? and W_Pl = MPLW
        Mto = self.MTOW_HACK
        W_Pl_B = self.MPLW_HACK                                               # Payload mass at point B             [kg]
        Mf = self.MTOW_HACK - W_Pl_B - self.OEW_HACK                          # Maximum fuel at maximum payload    [kg]
        W5W4_B = self.flight_profile_weights(Mf,Mto)
        R_B = self.Range(L,D,W5W4_B)


        # Point C: R = ?? and Wf = MFW, Wpl = MTOW - OEM - Maxfuel
        Mf = self.Max_fuel_mass                                               # Maximum fuel capacity of HACK      [kg]
        W_Pl_C = Mto - Mf - self.OEW_HACK                                     #Payload mass at point C             [kg]
        W5W4_C = self.flight_profile_weights(Mf, Mto)
        R_C = self.Range(L, D, W5W4_C)


        # Point D: R = ?? and W_Pl = 0--> W_fuel = Maxfuel
        Mf = self.Max_fuel_mass                                               # Maximum fuel capacity of HACK      [kg]
        W_PL_D = 0                                                            # Payload mass at point D            [kg]
        Mto = self.OEW_HACK + Mf + W_PL_D                                     # Take off weight at point D         [kg]
        W5W4_D = self.flight_profile_weights(Mf, Mto)
        R_D = self.Range(L, D, W5W4_D)

        Range_array = np.array([R_A,R_B,R_C,R_D])
        Payload_array = np.array([W_Pl_A,W_Pl_B,W_Pl_C,W_PL_D])
        plt.plot(Range_array,Payload_array,marker = '*',color = 'tab:red')
        plt.xlabel('Range [m]')
        plt.ylabel('Payload Mass [kg]')
        plt.show()


const = Constants()
AC_weights = compute_weight()
Performance = performance()
performance.payload_range_dia()

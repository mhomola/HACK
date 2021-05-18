#Following file computes a first estimation of the weight of the A320-HACK

from tank_sizing import liquid_H_tanks
from Subsystem_design.fuel_required import fuel_volume_calc
from Subsystem_design.common_constants import *
from Subsystem_design.common_constants import Constants
from fuel_constants import *
from math import pi
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

        self.fuel_mass = self.m_H2_center_w + self.m_H2_ext_f + self.m_k

        self.MPLW_HACK = self.MPLW_320neo                       # Maximum Payload weight of HACK                  [kg]
        self.MZFW_HACK = self.MPLW_HACK + self.OEW_HACK         # Maximum zero fuel weight of HACK                [kg]
        self.MTOW_HACK = self.MTOW_320neo                       # Maximum take-off weight of Hack                 [kg]


class Performance(Constants):
    #Function to plot the payload-range diagram
    def range(self):

    def payload_range_dia(self):
        #Payload as a funstion of range
        #Get Range equation, get Wf and Wi from fuel fractions
        # R = V/cp * L/D * Wi/Wf
        #PointA: R = 0 and W_Pl = MPLW
        #PointB: R = ?? and W_Pl = MPLW--> From MTOW calculate maxfuel for MPLW
        #PointC: R = ?? and Wf = MFW, Wpl = MTOW - OEM - Maxfuel
        #Point D: R = ?? and W_Pl = 0--> W_fuel = Maxfuel


const = Constants()
AC_weights = compute_weight()
const.fuselage_length(vol_eff=0.93,vol_fus=8.250)
VH2_center = const.V_H2_center_w
VH2_f = const.V_H2_ext_fus
AC_weights.Tank_mass(h2_vol_center= VH2_center,h2_vol_f=VH2_f)
print(AC_weights.Feeding_sys_m())
print(AC_weights.m_tank_ext_f)
print(AC_weights.m_tank_center_w)

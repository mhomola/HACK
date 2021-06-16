from Subsystem_design.common_constants import Constants
import matplotlib.pyplot as plt
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
import numpy as np

class Compute_weight(Constants):
    def __init__(self):
        super().__init__()
        self.m_k = 12883.86                                        # Kerosene mass                                  [kg]

    # def Tank_mass(self, h2_vol_center, h2_vol_f):
    # def Feeding_sys_m(self):
    # def Struc_m(self):
    def weight_break_down_HACK(self, h2_vol_center, h2_vol_f):
        self.OEW_HACK = self.OEW_320neo + ...

        self.Max_fuel_mass_capacity_HACK = self.m_k + ...

        self.MPLW_HACK = self.MPLW_320neo                           # Maximum Payload weight of HACK                  [kg]
        self.MZFW_HACK = self.MPLW_HACK + self.OEW_HACK             # Maximum zero fuel weight of HACK                [kg]
        self.MTOW_HACK = self.MTOW_320neo                           # Maximum take-off weight of Hack                 [kg]

        self.Max_fuel_at_max_PL_HACK = self.MTOW_HACK - self.MPLW_HACK - self.OEW_HACK  # Maximum fuel at maximum payload[kg]

        self.MRW_HACK = self.OEW_HACK + self.MPLW_HACK + self.Max_fuel_at_max_PL_HACK + self.Fuel_idel_taxi_take_off_HACK

class performance(Compute_weight):

    def Range(self,L_D_ratio,cruise_f_ratio,SFC):

        R = (self.V_cruise/(SFC*self.g_0)) * (L_D_ratio) *log(1/cruise_f_ratio)
        return R

#Following file computes a first estimation of the weight of the A320-HACK

from tank_sizing import liquid_H_tanks
import fuel_constants
from Subsystem_design.fuel_required import fuel_volume_calc
from Subsystem_design.common_constants import *
from Subsystem_design.common_constants import Constants


class compute_weight(Constants)

    def __init__(self):
        super().__init__()
        self.m_H2, self.m_tank,_ = liquid_H_tanks(h2_vol)          # Hydrogen mass and Tank mass                    [kg]
        self.m_k = 11595.96                                        # Kerosene mass                                  [kg]


    def Feeding_sys_m(self):
        self.m_r = 0                                             # Mass of feeding system per meter               [kg/m]
        self.FS_length = 66.5                                    # Length of feeding system(roughly)              [m]
        self.FD_mass = self.m_r * self.FS_length

        return self.FD_mass

    def Struc_m(self):
        # Added mass per meter of added fuselage                    [kg/m]
        self.extd_m_r = (self.OEW_321neo - self.OEW_320neo)/(self.l_f_321neo - self.l_f_320neo)
        self.struc_mass = self.extd_m_r * (self.l_f - self.l_f_320neo)      # Extra Structural mass               [kg]

        return self.struc_mass

    def weight_break_down_HACK(self):
        # Operational Empty Weight of Hack    [kg]
        self.OEW_HACK = self.OEW_320neo + self.m_tank + self.FD_mass + self.struc_mass

        self.fuel_mass = self.m_H2 + self.m_k

        self.MPLW_HACK = self.MPLW_320neo

constants = Constants

AC_weights = compute_weight(constants)
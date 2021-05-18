# import fuel_constants
# from subsystem_design.fuel_required import fuel_volume_calc
from subsystem_design.common_constants import Constants
import numpy as np



class Climb_perdormance(Constants):
    def __init__(self):
        super().__init__()
        self.max_cont_thrust = 9.02*10**6
        self.takeoff_thrust = 9.17*10**6

    def diff_density(self, input_altitude):
        alt_arr = [input_altitude - 1, input_altitude, input_altitude + 1]
        rho_ratios = np.zeros(len(alt_arr))

        for i, alt in enumerate(alt_arr):
            self.ISA_calculator(h_input=alt)
            rho_ratios[i] = np.sqrt(self.rho_0 / self.rho)

        diff_rho_ratio = np.gradient(rho_ratios)

        return  diff_rho_ratio[1]

    def rc_steady(self, C_D, C_L, V, W):
        rc_s = (self.max_cont_thrust/W) - (C_D*V/C_L)
        return rc_s

    def rc_real(self, rc_steady, V_EAS):
        dvdh = V_EAS * self.diff_density(input_altitude=self.cruise_altitude)
        rc = rc_steady/(1 + self.V_cruise * dvdh / self.g_0)
        return rc

if __name__ == '__main__':
    p = Climb_perdormance()
    p.diff_density(input_altitude=11000)


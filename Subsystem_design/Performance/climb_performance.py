# import fuel_constants
# from subsystem_design.fuel_required import fuel_volume_calc
from Subsystem_design.common_constants import Constants
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
import numpy as np



class Climb_perdormance(Constants):
    def __init__(self):
        super().__init__()
        self.max_cont_thrust = 120.64*10**3
        self.takeoff_thrust = 118.68*10**3

    def diff_density(self, input_altitude):
        alt_arr = [input_altitude - 1, input_altitude, input_altitude + 1]
        rho_ratios = np.zeros(len(alt_arr))

        for i, alt in enumerate(alt_arr):
            self.ISA_calculator(h_input=alt)
            rho_ratios[i] = np.sqrt(self.rho_0 / self.rho)

        diff_rho_ratio = np.gradient(rho_ratios)

        return  diff_rho_ratio[1]

    def rc_steady(self, W):
        aerodynamics = AerodynamicCharacteristics()
        aerodynamics.L_over_D_cruise()
        rc_s = (self.max_cont_thrust - aerodynamics.C_D_start_cruise_HACK*self.V_cruise**2*0.5*self.rho_0*self.S)/W
        return rc_s

    def rc_real(self, rc_steady, V_EAS):
        dvdh = V_EAS * self.diff_density(input_altitude=self.cruise_altitude)
        rc = rc_steady/(1 + self.V_cruise * dvdh / self.g_0)
        return rc

    def max_rc(self,W):
        aerodynamics = AerodynamicCharacteristics()
        aerodynamics.L_over_D_cruise()
        aerodynamics.wing_AR()
        CD= 4 * aerodynamics.C_D_0_clean_neo
        CL=np.sqrt(3 * aerodynamics.C_D_0_clean_neo*np.pi*aerodynamics.AR*aerodynamics.e)
        print(CD, CL)
        maxrc = ((self.max_cont_thrust / W)-(CD/CL))*np.sqrt(2*W/(self.S*self.rho_0*CL))
        return maxrc


if __name__ == '__main__':
    p = Climb_perdormance()
    p.diff_density(input_altitude=11000)
    print('Steady rate of climb = ', p.rc_steady(W=702364.53))
    print('Max RC = ',p.max_rc(W=702364.53))



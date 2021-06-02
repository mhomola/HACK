from Subsystem_design.common_constants import Constants
from Subsystem_design.Performance.flight_envelope import FlightEnvelope
from Subsystem_design.Performance.Weight_estimation import Compute_weight
import numpy as np
import matplotlib.pyplot as plt

class FuselageLoads(Constants):
    def __init__(self):
        super().__init__()

    def tail_weight(self):
        """
        Compute the weight of the vertical and horizontal wings using the Torenbeek method.
        :return:
        """
        K_h = 1.1  # Factor for variable incidence stabilizers
        w = Compute_weight()
        w.weight_break_down_HACK(h2_vol_center=11.34, h2_vol_f=30.69)  # todo: Check if values are ok
        fe = FlightEnvelope(altitude=0, W=w.OEW_HACK*self.g_0)
        fe.plot_envelope()

        V_D = fe.V_D_eq * 1.94384  # Dive speed in KEAS (Knots Equivalent Air Speed)
        sweep_05_h_rad = np.arctan(np.tan(self.sweep_LE_h * np.pi / 180) +
                                   self.c_r_h / 2 / (self.b_h / 2) * (self.taper_h - 1))
        sweep_05_v_rad = np.arctan(np.tan(self.sweep_LE_v * np.pi / 180) +
                                   self.c_r_v / 2 / (self.b_v / 2) * (self.taper_v - 1))
        S_h_ft2 = self.S_h * 3.28084**2
        S_v_ft2 = self.S_v * 3.28084**2
        W_h_lbs = K_h * S_h_ft2 * (3.81 * (S_h_ft2**0.2 * V_D) / (1000 * (np.cos(sweep_05_h_rad))**0.5) - 0.287)
        self.W_h = W_h_lbs * 0.453592  # Weight of the horizontal tail in kg

        K_v = 1.0  # Factot for fuselage mounted horizontal tails
        W_v_lbs = K_v * S_v_ft2 * (3.81 * (S_v_ft2**0.2 * V_D) / (1000 * (np.cos(sweep_05_v_rad))**0.5) - 0.287)
        self.W_v = W_v_lbs * 0.453592  # Weight of the vertical tail in kg

    def tail_loads(self):
        w = Compute_weight()
        w.weight_break_down_HACK(h2_vol_center=11.34, h2_vol_f=30.69)  # todo: Check if values are ok
        fe = FlightEnvelope(altitude=0, W=w.OEW_HACK*self.g_0)
        fe.plot_envelope()
        fe.max_tail_loads()

        self.L_h_max = fe.L_H_up  # Maximum load which the horizontal tail has to carry [N]
        self.L_v = fe.L_v  # Maximum load which the vertical tail has to carry [N]
        print(self.L_v, self.L_h_max)

    def load_diagrams(self, x):
        self.tail_weight()
        self.tail_loads()

        self.S_y




if __name__ == '__main__':
    fl = FuselageLoads()
    fl.tail_weight()
    fl.tail_loads()
    print('\n The weight of the vertical tail is ', fl.W_v, ' kg',
          '\n The weight of the horizontal tail is ', fl.W_h, ' kg')
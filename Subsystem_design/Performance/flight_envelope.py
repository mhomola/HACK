from Subsystem_design.common_constants import Constants
from Subsystem_design.Performance.Weight_estimation import Compute_weight
import numpy as np

class FlightEnvelope(Constants):

    def __init__(self, altitude):
        super().__init__()

        # Atmosphere
        self.ISA_calculator(h_input=altitude)

        # Conversions
        self.conv1 = 3.28084  # 1 m = conv1 * ft
        self.conv2 = 0.00194032  # 1 kg/m^3 = conv2 * slug/ft^3
        self.conv3 = 2.20462  # 1 kg = conv3 * lbs
        self.conv4 = 1.94384  # 1 m/s = conv4 * kts

    def FE_functions(self):
        self.n_ultimate()
        self.speeds()

    def n_ultimate(self):
        w = Compute_weight()
        w.weight_break_down_HACK(h2_vol_center=11.34, h2_vol_f=30.69)  # todo: Check if values are ok
        WTO_lbs = w.MTOW_HACK * self.conv3

        n_max = 2.1 + (24000 / (WTO_lbs + 10000))

        if n_max > 3.8:
            n_max = 3.8
        elif n_max < 2.5:
            n_max = 2.5

        self.n_max = n_max

    def speeds(self):

        # Dive speed
        V_D1 = self.V_cruise / 0.8   # CS25 Design Dive speed V_c <= 0.8V_D
        V_D2 = self.V_cruise + 0.05 * self.a  # CS25: The margin may not be reduced to less than 0.05M
        V_D = np.array([V_D1, V_D2])
        V_D_eq = V_D * np.qrt(self.rho / self.rho_0)

        if V_D1/self.a_c <= 1:
            V_D_eq.append(V_D1)
        if V_D2/self.a_c <= 1:
            V_D_eq.append(V_D2)

        self.V_D_eq = max(V_D_eq)
        self.V_D = V_D[np.where(V_D_eq == self.V_D_eq)[0]]

        # Maneuver speed
        V_A = np.sqrt(self.n_max / )


if __name__ == '__main__':
    fe = FlightEnvelope(altitude=11280)
    fe.FE_functions()

    print('\n n_max = ', fe.n_max,
          '\n V_D = ', fe.V_D)




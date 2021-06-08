from Subsystem_design.common_constants import Constants
from Subsystem_design.Performance.Weight_estimation import Compute_weight
import numpy as np
import matplotlib.pyplot as plt

class FlightEnvelope(Constants):

    def __init__(self, altitude, W):
        super().__init__()

        # Weight of aircraft at desired moment
        self.W = W  # [N]

        # Atmosphere
        self.ISA_calculator(h_input=altitude)

        # Conversions
        self.conv1 = 3.28084  # 1 m = conv1 * ft
        self.conv2 = 0.00194032  # 1 kg/m^3 = conv2 * slug/ft^3
        self.conv3 = 2.20462  # 1 kg = conv3 * lbs
        self.conv4 = 1.94384  # 1 m/s = conv4 * kts

    def plot_maneuver_envelope(self):

        w = Compute_weight()
        w.weight_break_down_HACK(h2_vol_center=11.34, h2_vol_f=30.69)  # todo: Check if values are ok

        W_TO = w.MTOW_HACK * self.g_0
        W_land = self.MLW_320neo * self.g_0  # todo: Check if value is ok
        W_TO_lbs = W_TO * self.conv3 / self.g_0

        # n_max
        n_max = 2.1 + (24000 / (W_TO_lbs + 10000))
        if n_max > 3.8:
            n_max = 3.8
        elif n_max < 2.5:
            n_max = 2.5
        self.n_max = n_max

        # Cruise speed
        V_C = self.V_cruise
        V_C_eq = V_C * np.sqrt(self.rho / self.rho_0)
        print(V_C, V_C_eq)

        # Dive speed
        V_D1 = self.V_cruise / 0.8   # CS25 Design Dive speed V_c <= 0.8V_D
        V_D2 = self.V_cruise + 0.05 * self.a  # CS25: The margin may not be reduced to less than 0.05M
        V_D_lst = np.array([V_D1, V_D2])
        V_D_eq_lst = V_D_lst * np.sqrt(self.rho / self.rho_0)
        V_D_eq = np.max(V_D_eq_lst[np.where(V_D_eq_lst/self.a <= 1.0)])
        V_D = V_D_lst[np.where(V_D_eq_lst == V_D_eq)[0]]
        print(V_D_eq)

        # Maneuver speed
        V_A = np.sqrt(self.n_max / self.C_L_max_clean * self.W / (0.5 * self.rho * self.S))
        V_A_eq = V_A * np.sqrt(self.rho / self.rho_0)
        V_A_land = np.sqrt(2.0 / self.C_L_max_land * self.W / (0.5 * self.rho * self.S))
        V_A_land_eq = V_A_land * np.sqrt(self.rho / self.rho_0)
        V_A_TO = np.sqrt(2.0 / self.C_L_max_land * self.W / (0.5 * self.rho * self.S))
        V_A_TO_eq = V_A_TO * np.sqrt(self.rho / self.rho_0)
        print(V_A_eq, V_A_land_eq)

        # Stall speeds
        V_s1 = np.sqrt(1.0 / self.C_L_max_clean * self.W / (0.5 * self.rho * self.S))
        V_s1_eq = V_s1 * np.sqrt(self.rho / self.rho_0)
        V_s0_land = np.sqrt(self.W / (self.C_L_max_land * 0.5 * self.rho * self.S))
        V_s0_land_eq = V_s0_land * np.sqrt(self.rho / self.rho_0)
        V_s0_TO = np.sqrt(self.W / (self.C_L_max_TO * 0.5 * self.rho * self.S))
        V_s0_TO_eq = V_s0_TO * np.sqrt(self.rho / self.rho_0)
        print(V_s1_eq)

        # Design  wing-flap  speeds
        V_F1 = 1.6 * np.sqrt(W_TO / (self.C_L_max_TO * 0.5 * self.rho * self.S))
        V_F2 = 1.8 * np.sqrt(W_land / (self.C_L_max_land * 0.5 * self.rho * self.S))
        V_F = max(V_F1, V_F2)
        V_F_eq = V_F * np.sqrt(self.rho / self.rho_0)
        print(V_F_eq)

        # Negative loads
        V_H = np.sqrt(1 / self.C_L_max_clean * self.W / (self.rho * 0.5* self.S))
        V_H_eq = V_H * np.sqrt(self.rho / self.rho_0)
        print(V_H_eq)

        # cl_max curve
        v_clmax_curve = np.linspace(0, V_A_eq, 100)
        n_clmax_curve = self.C_L_max_clean * 0.5 * self.rho_0 * v_clmax_curve**2 * self.S / self.W

        # n_max curve
        v_max_curve = np.linspace(V_A_eq, V_D_eq, 5)
        n_max_curve = np.ones(len(v_max_curve)) * self.n_max

        # dive speed curve
        v_dive_curve = np.array([V_D_eq, V_D_eq])
        n_dive_curve = np.array([self.n_max, 0])

        # cl_max_land curve
        v_clmax_land_curve = np.linspace(0, V_A_land_eq, 100)
        n_clmax_land_curve = self.C_L_max_land * 0.5 * self.rho_0 * v_clmax_land_curve**2 * self.S / self.W
        v_clmax_land_curve = np.hstack((v_clmax_land_curve, np.linspace(V_A_land_eq, V_F_eq, 5)))
        n_clmax_land_curve = np.hstack((n_clmax_land_curve, np.ones(5) * 2.0))

        # negative cl_max curve
        v_neg_clmax_curve = np.linspace(0, V_H_eq, 100)
        n_neg_clmax_curve = - self.C_L_max_clean * 0.5 * self.rho_0 * v_neg_clmax_curve**2 * self.S / self.W

        # n_min curve
        v_min_curve = np.array([V_H_eq, V_C_eq, V_D_eq])
        n_min_curve = np.array([-1., -1., 0 ])

        fig, ax = plt.subplots(1,1)
        ax.plot(v_clmax_curve, n_clmax_curve, color='#5391d6')
        ax.plot(v_max_curve, n_max_curve, color='#5391d6')
        ax.plot(v_dive_curve, n_dive_curve, color='#5391d6')
        ax.plot(v_clmax_land_curve, n_clmax_land_curve, color='#5391d6')
        ax.plot(v_neg_clmax_curve, n_neg_clmax_curve, color='#5391d6')
        ax.plot(v_min_curve, n_min_curve, color='#5391d6')
        ax.axhline(y=0, color='k', linewidth=0.75)
        ax.axvline(x=0, color='k', linewidth=0.75)
        plt.show()


if __name__ == '__main__':
    fe = FlightEnvelope(altitude=11280, W=73500 * 9.80665)
    fe.plot_maneuver_envelope()
    print('\n n_max = ', fe.n_max,
          '\n V_D = ', fe.V_D)




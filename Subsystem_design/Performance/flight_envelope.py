from Subsystem_design.common_constants import Constants
from Subsystem_design.Performance.Weight_estimation import Compute_weight
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
import numpy as np
import matplotlib.pyplot as plt

class FlightEnvelope(Constants):

    def __init__(self, altitude, W, ):
        super().__init__()

        # Weight of aircraft at desired moment
        self.W = W  # [N]

        # Atmosphere
        self.alt = altitude
        self.ISA_calculator(h_input=altitude)

        # Conversions
        self.conv1 = 3.28084  # 1 m = conv1 * ft
        self.conv2 = 0.00194032  # 1 kg/m^3 = conv2 * slug/ft^3
        self.conv3 = 2.20462  # 1 kg = conv3 * lbs
        self.conv4 = 1.94384  # 1 m/s = conv4 * kts

    def create_maneuver_envelope(self):

        W_TO = self.MTOW_320neo * self.g_0
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
        self.V_C = self.V_cruise
        self.V_C_eq = self.V_C * np.sqrt(self.rho / self.rho_0)

        # Dive speed
        V_D1 = self.V_cruise / 0.8   # CS25 Design Dive speed V_c <= 0.8V_D
        V_D2 = self.V_cruise + 0.05 * self.a  # CS25: The margin may not be reduced to less than 0.05M
        V_D_lst = np.array([V_D1, V_D2])
        V_D_eq_lst = V_D_lst * np.sqrt(self.rho / self.rho_0)
        self.V_D_eq = np.max(V_D_eq_lst[np.where(V_D_eq_lst/self.a <= 1.0)])
        self.V_D = V_D_lst[np.where(V_D_eq_lst == self.V_D_eq)[0]]

        # Maneuver speed
        self.V_A = np.sqrt(self.n_max / self.C_L_max_clean * self.W / (0.5 * self.rho * self.S))
        V_A_eq = self.V_A * np.sqrt(self.rho / self.rho_0)
        V_A_land = np.sqrt(2.0 / self.C_L_max_land * self.W / (0.5 * self.rho * self.S))
        V_A_land_eq = V_A_land * np.sqrt(self.rho / self.rho_0)
        V_A_TO = np.sqrt(2.0 / self.C_L_max_land * self.W / (0.5 * self.rho * self.S))
        V_A_TO_eq = V_A_TO * np.sqrt(self.rho / self.rho_0)

        # Stall speeds
        self.V_s1 = np.sqrt(1.0 / self.C_L_max_clean * self.W / (0.5 * self.rho * self.S))
        self.V_s1_eq = self.V_s1 * np.sqrt(self.rho / self.rho_0)
        V_s0_land = np.sqrt(self.W / (self.C_L_max_land * 0.5 * self.rho * self.S))
        V_s0_land_eq = V_s0_land * np.sqrt(self.rho / self.rho_0)
        V_s0_TO = np.sqrt(self.W / (self.C_L_max_TO * 0.5 * self.rho * self.S))
        V_s0_TO_eq = V_s0_TO * np.sqrt(self.rho / self.rho_0)

        # Design  wing-flap  speeds
        V_F1 = 1.6 * np.sqrt(W_TO / (self.C_L_max_TO * 0.5 * self.rho * self.S))
        V_F2 = 1.8 * np.sqrt(W_land / (self.C_L_max_land * 0.5 * self.rho * self.S))
        V_F = max(V_F1, V_F2)
        V_F_eq = V_F * np.sqrt(self.rho / self.rho_0)

        # Negative loads
        V_H = np.sqrt(1 / self.C_L_max_clean * self.W / (self.rho * 0.5* self.S))
        V_H_eq = V_H * np.sqrt(self.rho / self.rho_0)

        # cl_max curve
        self.v_clmax_curve = np.linspace(0, V_A_eq, 100)
        self.n_clmax_curve = self.C_L_max_clean * 0.5 * self.rho_0 * self.v_clmax_curve**2 * self.S / self.W

        # n_max curve
        self.v_max_curve = np.linspace(V_A_eq, self.V_D_eq, 5)
        self.n_max_curve = np.ones(len(self.v_max_curve)) * self.n_max

        # dive speed curve
        self.v_dive_curve = np.array([self.V_D_eq, self.V_D_eq])
        self.n_dive_curve = np.array([self.n_max, 0])

        # cl_max_land curve
        self.v_clmax_land_curve = np.linspace(0, V_A_land_eq, 100)
        self.n_clmax_land_curve = self.C_L_max_land * 0.5 * self.rho_0 * self.v_clmax_land_curve**2 * self.S / self.W
        self.v_clmax_land_curve = np.hstack((self.v_clmax_land_curve, np.linspace(V_A_land_eq, V_F_eq, 5)))
        self.n_clmax_land_curve = np.hstack((self.n_clmax_land_curve, np.ones(5) * 2.0))

        # negative cl_max curve
        self.v_neg_clmax_curve = np.linspace(0, V_H_eq, 100)
        self.n_neg_clmax_curve = - self.C_L_max_clean * 0.5 * self.rho_0 * self.v_neg_clmax_curve**2 * self.S / self.W

        # n_min curve
        self.v_min_curve = np.array([V_H_eq, self.V_C_eq, self.V_D_eq])
        self.n_min_curve = np.array([-1., -1., 0 ])

    def create_gust_envelope_main_wing(self, ):
        ae = AerodynamicCharacteristics()
        ae.aero_functions(AoA_cruise=2)

        W_TO = self.MTOW_320neo * self.g_0
        W_land = self.MLW_320neo * self.g_0  # todo: Check if value is ok
        MZFW = self.MZFW_320neo * self.g_0
        W_TO_lbs = W_TO * self.conv3 / self.g_0

        Fgz_SL = 1 - self.max_altitude / 76200
        R1 = W_land / W_TO
        R2 = MZFW / W_TO
        Fgm_SL = np.sqrt(R2 * np.tan(np.pi * R1 / 4) )
        F_g_SL = 0.5 * (Fgz_SL + Fgm_SL)
        wing_load = self.W * self.conv3 / (self.S * self.conv1**2)
        mu_g = 2 * wing_load / (self.rho * self.conv2 * ae.mac * self.conv1 *
                                ae.CL_alpha_w * self.g_0 * self.conv1)
        K_g = 0.88 * mu_g / (5.3 + mu_g)

        def gust_profile(U_ref, F_g):

            H_arr = np.linspace(9, 107, 500)
            U_arr = np.zeros(len(H_arr))
            for i, H in enumerate(H_arr):
                Uds = U_ref * F_g * (H / 107)**(1/6)
                # s = np.linspace(0, 2*H, 500)
                # U_profile = (Uds / 2) * (1 - np.cos(np.pi * s / H))
                U_arr[i] = Uds
                # plt.plot(s, U_profile)
            # plt.show()
            U = np.max(U_arr)  # meters per second equivalent airspeed
            H = H_arr[np.where(U_arr==U)[0]]

            return U, H

        def ref_gust(alt):
            """

            :param alt: altitude to compute the gust
            :return: reference gust velocity in equivalent airspeed [m/s]
            """
            if alt < 4572:
                U_ref_BC = (13.41 - 17.07) / 4572 * alt + 17.07
            else:
                U_ref_BC = (6.36 - 13.41) / (18288 - 4572) * (alt - 4572) + 13.41
            U_ref_D = 0.5 * U_ref_BC

            return U_ref_BC, U_ref_D

        def F_g(alt, F_g_SL=F_g_SL, max_alt=self.max_altitude):

            F_g = (1 - F_g_SL) / max_alt * alt + F_g_SL
            return F_g

        def V_B_eq(U_ref, V_C=self.V_C_eq*self.conv4, K_g=K_g, a=ae.CL_alpha_w, wing_load=wing_load,
                   V_s1=self.V_s1_eq):
            """

            :param U_ref: in feet per second equivalent airspeed [ft/s]
            :param V_C:  in knots equivalent airspeed
            :param K_g: has no units
            :param a: wing lift curve slope in [1/rad]
            :param wing_load: in pounds per square foot
            :param V_s1: the unit of V_s1 defines the unit of V_B, use [m/s] EAS please
            :return: V_B_eq in equivalent airpseed [m/s]
            """
            return V_s1 * np.sqrt(1 + K_g * U_ref * V_C * a / (498 * wing_load))

        def n_g(V, U, H):

            CL_alpha_w = ae.lift_gradient_calc(M=V/self.a, AR=ae.AR, sweep_05=self.sweep_05)
            lam = self.W / self.S * 1 / CL_alpha_w * 2 / (self.rho * V * self.g_0)
            omega = np.pi * V / H
            t = np.linspace(0, 2*np.pi/omega, 500)
            Delta_n_s = U / (2 * self.g_0) * (omega * np.sin(omega * t) + 1 / (1 + (omega * lam)**(-2)) *
                                             (1 / lam * np.exp(-t/lam) - 1 / lam * np.cos(omega*t) -
                                              omega * np.sin(omega*t)))
            # plt.plot(t, Delta_n_s)
            # plt.show()
            return 1 + np.max(Delta_n_s), 1 - np.max(Delta_n_s)

        # Compute delta n for V_B, V_C and V_D
        U_ref_BC, U_ref_D = ref_gust(alt=self.alt)  # [m/s]
        self.U_BC, H_BC = gust_profile(U_ref=U_ref_BC*np.sqrt(self.rho_0 / self.rho), F_g=F_g(alt=self.alt))
        self.U_BC_eq = self.U_BC * np.sqrt(self.rho / self.rho_0)
        U_D, H_D = gust_profile(U_ref=U_ref_D*np.sqrt(self.rho_0 / self.rho), F_g=F_g(alt=self.alt))
        self.V_B_eq = V_B_eq(U_ref=U_ref_BC*self.conv1)  # [m/s]
        self.V_B = self.V_B_eq * np.sqrt(self.rho_0 / self.rho)
        self.n_B_pos, self.n_B_neg = n_g(V=self.V_B, U=self.U_BC, H=H_BC)
        self.n_C_pos, self.n_C_neg = n_g(V=self.V_C, U=self.U_BC, H=H_BC)
        self.n_D_pos, self.n_D_neg = n_g(V=self.V_D, U=U_D, H=H_D)

    def plot_envelope(self):

        self.create_maneuver_envelope()
        self.create_gust_envelope_main_wing()
        fig, ax = plt.subplots(1,1)
        ax.plot(self.v_clmax_curve, self.n_clmax_curve, color='#5391d6')
        ax.plot(self.v_max_curve, self.n_max_curve, color='#5391d6')
        ax.plot(self.v_dive_curve, self.n_dive_curve, color='#5391d6')
        ax.plot(self.v_clmax_land_curve, self.n_clmax_land_curve, color='#5391d6')
        ax.plot(self.v_neg_clmax_curve, self.n_neg_clmax_curve, color='#5391d6')
        ax.plot(self.v_min_curve, self.n_min_curve, color='#5391d6')
        ax.plot([self.V_B_eq, 0, self.V_B_eq], [self.n_B_neg, 1, self.n_B_pos], color='#f0067b')
        ax.plot([self.V_B_eq, self.V_C_eq], [self.n_B_pos, self.n_C_pos], color='#f0067b')
        ax.plot([self.V_B_eq, self.V_C_eq], [self.n_B_neg, self.n_C_neg], color='#f0067b')
        ax.plot([self.V_C_eq, self.V_D_eq], [self.n_C_pos, self.n_D_pos], color='#f0067b')
        ax.plot([self.V_C_eq, self.V_D_eq, self.V_D_eq], [self.n_C_neg, self.n_D_neg, self.n_D_pos], color='#f0067b')
        ax.axhline(y=0, color='k', linewidth=0.75)
        ax.axvline(x=0, color='k', linewidth=0.75)
        plt.show()

    def max_tail_loads(self):

        self.create_gust_envelope_main_wing()
        ae = AerodynamicCharacteristics()
        ae.aero_functions(AoA_cruise=2)

        # Vertical tail loads

        I_Z = 0.3 * 10**(12)  # mm^4 (assumed from Thesis of Ilhan
        CL_v_beta = 2 * np.pi / (1 + 3 / (self.AR_v * np.cos(self.sweep_LE_v * np.pi / 180)))  # Torenbeek
        l_v = (self.X_root_vtail + ae.x_mac_v) - (self.X_root_wing + ae.x_mac)
        mu_v = 2 * I_Z / (self.rho * ae.mac_v * CL_v_beta * self.S_v * l_v**2)
        k_g = 0.88 * mu_v / (5.3 + mu_v)

        self.L_v = k_g * 0.5 * self.rho_0 * self.U_BC_eq * self.V_C_eq * self.S_v * CL_v_beta

        # Horizontal tail load

        CL_h_delta = ae.CL_alpha_h * np.sqrt(self.S_elevator / self.S_h)
        print(CL_h_delta)
        print(self.V_C, self.rho)
        self.L_H_up = 0.9 * 0.5 * self.rho * self.V_A**2 * self.S_h * CL_h_delta * self.max_elevator_deflection_nu *\
                      np.pi / 180
        self.L_H_down = 0.9 * 0.5 * self.rho * self.V_A**2 * self.S_h * CL_h_delta * self.max_elevator_deflection_nd * \
                        np.pi / 180




if __name__ == '__main__':
    fe = FlightEnvelope(altitude=0, W=45000*9.80665)
    fe.plot_envelope()
    print('\n n_max = ', fe.n_max,
          '\n V_D = ', fe.V_D)




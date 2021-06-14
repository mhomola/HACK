import Subsystem_design.Tank_Design.Mechanical_Design as Mechanical_Design
import Subsystem_design.Tank_Design.Materials as Materials
from Subsystem_design.fuel_required import V_H2
import numpy as np
from Subsystem_design.common_constants import Constants
import matplotlib.pyplot as plt


class spacial_constraints():
    def __init__(self, volume, width, height):
        """
        The sizes of the current A320 neo will represent dimensional constraints for fitting in the LH2 tanks.
        :param length:  [m]
        :param width:   [m]
        :param height:  [m]
        """

        self.volume = volume
        self.width = width
        self.height = height


volume_req = V_H2/1000

volume_pod = volume_req/2

lengths = []
diameters = []
drags = []
masses = []

d_i = 1.85
for i in range(0, 170):

    pod = spacial_constraints(volume=volume_pod, width=d_i, height=d_i)

    pod_tank = Mechanical_Design.OnlyPods(constraints=pod, dp=Mechanical_Design.dp, s_a=Mechanical_Design.s_a,
                                          e_w=Mechanical_Design.e_w, material_insulation=Materials.MLI
                                          , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                          rho=Mechanical_Design.rho, t_tank=Mechanical_Design.t_tank,
                                          dt=Mechanical_Design.dt, p_tank=Mechanical_Design.p_tank)
    pod_tank.tank_design()

    weight_addition = pod_tank.mass_tank * 2
    volume_total = pod_tank.inner_vol_inner_wall * 2

    l_wing_pod = pod_tank.length
    d_wing_pod = pod_tank.r4 * 2


    class AerodynamicCharacteristics(Constants):

        def __init__(self):
            super().__init__()

        def aero_functions(self, AoA_cruise):
            self.wing_MAC()
            self.wing_AR()
            self.h_tail_MAC()
            self.v_tail_MAC()
            self.drag_increase_cruise(AoA_cruise=AoA_cruise)
            self.L_over_D_cruise()

        def wing_MAC(self):
            """
            The wing MAC is computed using the ADSEE-II slides. The wing surface is composed out of two trapezoids.
            First the mac of these two trapezoids is computed separately and after the global MAC is computed.
            :return: The length and position of the wing's MAC
            """
            # Inner trapezoid
            mac_in = (2/3) * self.c_root * (1 + self.taper_in + self.taper_in**2) / (1 + self.taper_in)  # MAC [m]
            y_mac_in = (self.b_in / 6) * (1 + 2 * self.taper_in) / (1 + self.taper_in)  # y position of MAC [m]
            x_mac_in = y_mac_in * np.tan(self.sweep_LE * np.pi / 180)  # x location of LEMAC from start of root chord [m]
            S_in = self.b_in * (self.c_root + self.c_kink_out) / 2  # Surface area [m^2]

            # Outer trapezoid
            mac_out = (2/3) * self.c_kink_out * (1 + self.taper_out + self.taper_out**2) / (1 + self.taper_out)
            y_mac_out = (self.b_out / 6) * (1 + 2 * self.taper_out) / (1 + self.taper_out)
            x_mac_out = y_mac_out * np.tan(self.sweep_LE * np.pi / 180)
            S_out = self.b_out * (self.c_kink_out + self.c_tip) / 2

            # Complete wing
            self.mac = (S_out * mac_out + S_in * mac_in) / (S_out + S_in)  # MAC of the complete wing [m]
            self.y_mac = (y_mac_in * S_in + (0.5 * self.b_in + y_mac_out) * S_out) / (S_out + S_in)  # y position of MAC [m]
            self.x_mac = self.y_mac * np.tan(self.sweep_LE * np.pi / 180)  # x location of LEMAC from root chord [m]

        def wing_AR(self):
            """
            Compute the aspect ratio of the wing including sharklets with ADSEE-II folrmula
            :return: The AR of the main wing
            """
            AR_no_sharklets = self.b**2 / self.S
            self.AR = AR_no_sharklets + 1.9 * (self.h_sharklet / self.b) * AR_no_sharklets

        def h_tail_MAC(self):
            """
            The MAC is computed using the ADSEE-II slides.
            :return: The length and position of the horizontal wing's MAC
            """
            self.mac_h = (2/3) * self.c_r_h * (1 + self.taper_h + self.taper_h**2) / (1 + self.taper_h)
            self.y_mac_h = (self.b_h / 6) * (1 + 2 * self.taper_h) / (1 + self.taper_h)
            self.x_mac_h = self.y_mac_h * np.tan(self.sweep_LE_h * np.pi / 180)


        def v_tail_MAC(self):
            """
            The MAC is computed using the ADSEE-II slides.
            :return: The length and position of the horizontal wing's MAC
            """
            self.mac_v = (2/3) * self.c_r_v * (1 + self.taper_v + self.taper_v**2) / (1 + self.taper_v)
            self.y_mac_v = (self.b_v / 6) * (1 + 2 * self.taper_v) / (1 + self.taper_v)
            self.x_mac_v = self.y_mac_v * np.tan(self.sweep_LE_v * np.pi / 180)

            # print('\n MAC of the vertical tail = ', self.mac_v, ' m',
            #       '\n y position of the MAC = ', self.y_mac_v, ' m',
            #       '\n x position of the LEMAC measured from the start of the root chord = ', self.x_mac_v, ' m')

        def Roskam_drag_prediction_cruise(self, rho, u1, l_cockpit, l_cabin, l_tail, AoA, S_b_fus, height_f, width_f):
            """
            This function computes the drag from the fuselage according to the procedure given by Roskam at transonic
            conditions.
            :param rho: Air density [kg/m^3]
            :param u1: Airspeed [m/s]
            :param l_f: Length of the entire fuselage [m]
            :param l_cockpit: Length of the cockpit [m]
            :param l_cabin: Length of the cabin [m]
            :param l_tail: Length of the tail [m]
            :param AoA: Angle of attack of the aircraft [deg]
            :return: The zero lift drag coefficient of the fuselage and the drag coefficient of the fusselage due to lift
            """
            # Wetted area as computed in ADSEE-II
            S_fus = np.pi * 0.25 * height_f * width_f
            self.d_f = np.sqrt(4 * S_fus / np.pi)
            S_wet_fus = np.pi * self.d_f / 4 * \
                        (1 / (3 * l_cockpit**2) * ((4 * l_cockpit**2 + self.d_f**2 / 4) - self.d_f**3 / 8)
                         - self.d_f + 4 * l_cabin + 2 * np.sqrt(l_tail**2 + self.d_f**2 / 4))

            # Compute zero lift drag for M = 0.6 for fuselage exclusive of base
            l_f = l_cockpit + l_cabin + l_tail
            R_n_fus = rho * u1 * l_f / self.visc

            if l_f > 20:
                R_wf = 1.015  # The wing/fuselage iterference factor from Figure 4.1 in Roskam-VI
                C_f_fus = 0.0016  # The turbulent flat plate skin friction coefficient from Figure 4.3 in Roskam-VI
            else:
                R_wf = 1
                C_f_fus = 0.0018


            ld = l_f / self.d_f
            C_D_o_fus_exc_base = R_wf * C_f_fus * (1 + 60 / ld**3 + 0.0025 * ld) * S_wet_fus / self.S
            # Compute the fuselage base drag coefficient
            bf = np.sqrt(4 / np.pi * S_b_fus) / self.d_f
            C_D_b_fus = 0.09 * bf**2

            # The zero lift drag coefficient of the fuselage becomes:
            C_D_0_fus = C_D_o_fus_exc_base + C_D_b_fus

            # The lift induced drag
            C_D_L_fus = (AoA * np.pi / 180)**2 * S_b_fus / self.S

            return C_D_0_fus, C_D_L_fus

        def drag_increase_cruise(self, AoA_cruise):
            """
            Compute the drag increase due to the new fuselage length
            :param AoA_cruise: Angle of attack during cruise for the lift induced drag [deg]
            :return: The most important parameter which is computed is C_D_0_HACK
            """
            # Compute density at cruise altitude
            self.ISA_calculator(h_input=self.cruise_altitude)

            # A320neo
            self.C_D_0_fus_neo, _ = self.Roskam_drag_prediction_cruise(rho=self.rho, u1=self.M*self.a,
                                                                       l_cockpit=self.l_cockpit_320neo,
                                                                       l_cabin=self.l_cabin_320neo,
                                                                       l_tail=self.l_tail_320neo, AoA=AoA_cruise,
                                                                       S_b_fus=self.S_b_fus, height_f=self.height_f,
                                                                       width_f=self.width_f)

            # A320-HACK
            self.l_tank_nose, self.l_tank_body, self.l_tank_tail = d_wing_pod / 2, l_wing_pod - d_wing_pod, 2*d_wing_pod
            self.d_tank = d_wing_pod
            self.C_D_0_tank_HACK, _ = self.Roskam_drag_prediction_cruise(rho=self.rho, u1=self.M*self.a,
                                                                         l_cockpit=self.l_tank_nose,
                                                                         l_cabin=self.l_tank_body,
                                                                         l_tail=self.l_tank_tail, AoA=AoA_cruise,
                                                                         S_b_fus=0, height_f=self.d_tank,
                                                                         width_f=self.d_tank)

            self.C_D_0_tank_sys_HACK = 2 * 1.3 * self.C_D_0_tank_HACK
            self.C_D_0_HACK = self.C_D_0_clean_neo + self.C_D_0_tank_sys_HACK
            self.C_D_0_HACK_TO = self.C_D_0_TO_neo + self.C_D_0_tank_sys_HACK
            self.C_D_0_HACK_land = self.C_D_0_land_neo + self.C_D_0_tank_sys_HACK
            self.C_D_O_HACK_approach = self.C_D_0_approach_neo + self.C_D_0_tank_sys_HACK
            self.C_D_0_HACK_approach = self.C_D_0_approach_neo + self.C_D_0_tank_sys_HACK
            self.C_D_0_HACK_taxi = self.C_D_0_taxi_neo + self.C_D_0_tank_sys_HACK

        def L_over_D_cruise(self):
            self.ISA_calculator(h_input=self.cruise_altitude)
            self.wing_AR()
            self.drag_increase_cruise(2)
            # self.DATCOM_podded_tanks(c_loc=4, x_fwd=1, d_tank=self.d_tank, alpha=2)

            W_start_cruise = self.MTOW_320neo * (0.995 * 0.98)
            V = self.M * self.a

            self.C_L_start_cruise = W_start_cruise * self.g_0 / (0.5 * self.rho * V ** 2 * self.S)

            self.C_D_start_cruise_neo = self.C_D_0_clean_neo + self.C_L_start_cruise ** 2 / (np.pi * self.AR * self.e)
            self.C_D_start_cruise_HACK = self.C_D_0_HACK + self.C_L_start_cruise ** 2 / (np.pi * self.AR * self.e)

            self.D_start_cruise_neo = self.C_D_start_cruise_neo * 0.5 * self.rho * V ** 2 * self.S
            self.D_start_cruise_HACK = self.C_D_start_cruise_HACK * 0.5 * self.rho * V ** 2 * self.S

            self.L_D_ratio_neo = self.C_L_start_cruise / self.C_D_start_cruise_neo
            self.L_D_ratio_HACK = self.C_L_start_cruise / self.C_D_start_cruise_HACK

    ae = AerodynamicCharacteristics()
    ae.aero_functions(AoA_cruise=2)

    lengths.append(l_wing_pod)
    diameters.append(d_wing_pod)
    # drags.append((ae.C_D_0_HACK / ae.C_D_0_clean_neo - 1) * 100)
    drags.append((ae.D_start_cruise_HACK / ae.D_start_cruise_neo - 1) * 100)
    masses.append(pod_tank.mass_tank)

    d_i = d_i + 0.002


if __name__ == '__main__':

    fig1, ax1 = plt.subplots(1, 1)
    ax1.plot(diameters, drags, label='Influence of size on drag')
    ax1.set_xlabel(r'$D$ in m', size=15)
    ax1.set_ylabel(r'$C_D$ Increase in %', size=15)
    plt.grid()


    fig2, ax2 = plt.subplots(1, 1)
    ax2.plot(diameters, lengths, label='Diameter - Length relationship')
    ax2.set_xlabel(r'$D$ in m', size=15)
    ax2.set_ylabel(r'$L$ in m', size=15)
    plt.grid()


    fig3, ax3 = plt.subplots(1, 1)
    ax3.plot(diameters, masses, label='Diameter - Mass relationship')
    ax3.set_xlabel(r'$D$ in m', size=15)
    ax3.set_ylabel(r'$Mass$ in kg', size=15)
    plt.grid()

    plt.show()


    # ae = AerodynamicCharacteristics()
    # ae.aero_functions(AoA_cruise=2)

    # print('\n Assuming the C_D_0 of the A320neo during cruise is = ', ae.C_D_0_clean_neo,
    #       '\n The C_D_0 of the A320-HACK now becomes = ', ae.C_D_0_HACK,
    #       '\n That is a ', (ae.C_D_0_HACK / ae.C_D_0_clean_neo - 1) * 100, '% increase')





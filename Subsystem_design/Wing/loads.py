from Subsystem_design.common_constants import Constants
from Subsystem_design.Performance.flight_envelope import FlightEnvelope
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import scipy.integrate as spint
from Subsystem_design.Wing.Kerosene_distirbution import kerosene_calc
from Subsystem_design.Tank_Design.Main_PreliminaryTank import d_wing_pod

class Loads_w(Constants):

    def __init__(self):
        super().__init__()
        self.dx = 0.01
        self.m1 = (self.c_kink_out - self.c_root) / (0.5 * self.b_in)
        self.m2 = (self.c_tip - self.c_kink_out) / (0.5 * self.b_out)
        self.xL_c = 0.163  # Distance from the quarter chord to the center of the wing box over the local chord length
        self.T_arm = 1.5  # Torque arm of the engine thrust and drag
        self.W_eng_arm = 3.5  # Moment arm from the weight of the engine
        self.max_T = 140000  # Max thrust [N]
        self.W_wing_arm_c = 0.085  # Moment arm of wing weight over the local chord length

    def compute_loads(self):

        altitudes = np.linspace(0, self.cruise_altitude, 2)
        W_TO = (self.MTOW_320neo) * self.g_0  # todo: Check if value is ok
        W_empty = (self.OEW_320neo + 1000) * self.g_0  # todo: Check if value is ok
        weights = np.linspace(W_empty, W_TO, 2)
        self.L_max = 0
        index_i, index_j = 0, 0

        for i, weight in enumerate(weights):
            for j, alti in enumerate(altitudes):
                fe = FlightEnvelope(altitude=alti, W=weight)
                fe.create_maneuver_envelope()
                fe.create_gust_envelope_main_wing()
                # fe.plot_envelope()
                n_max = max(2.5, fe.n_C_pos)
                L_tot = weight * n_max
                if L_tot > self.L_max:
                    index_i, index_j = i, j
                    self.L_max = L_tot

                # print('\n At an altitude of ', alti, ' m, a weight of ', weight, ' N, the max lift is ', L_tot)

        print('\n At an altitude of ', altitudes[index_j], ' m, a weight of ', weights[index_i],
              ' N, the max lift is ', self.L_max)


        self.ISA_calculator(h_input=altitudes[index_j])
        fe = FlightEnvelope(altitude=altitudes[index_j], W=weights[index_i])
        fe.create_maneuver_envelope()
        ae = AerodynamicCharacteristics()
        ae.aero_functions(AoA_cruise=2)
        self.q_crit = 0.5 * self.rho * fe.V_A**2

        x_arr = np.linspace(0.5*self.width_f, self.b/2, 1000)
        c_arr = np.zeros(len(x_arr))
        for i, x in enumerate(x_arr):
            c_arr[i] = self.chord(x=x)
        S_exp = 2 * spint.simps(y=c_arr, x=x_arr)  # Exposed surface ara [m^2]
        self.C_L_crit = self.L_max / (self.q_crit * S_exp)
        self.C_D_crit = ae.C_D_0_HACK - ae.C_D_0_fus_neo - ae.C_D_0_tank_sys_HACK - ae.C_D_o_engine - ae.C_D_0_Vtail - \
                        ae.C_D_0_Htail + self.C_L_crit / (np.pi * ae.AR * self.e)

        self.W_body = (self.MZFW_320neo - self.Wing_Weight_320neo - 2*self.W_engine - 2*self.pod_tank_mass) * self.g_0
        self.WbS = self.W_body / (self.S - S_exp)
        self.WwS = self.Wing_Weight_320neo * self.g_0 / self.S

    def Lift(self, x):
        c = self.chord(x=x)
        L_prime = self.C_L_crit * self.q_crit * c
        return L_prime

    def Drag_w(self, x):
        c = self.chord(x=x)
        D_prime = self.C_D_crit * self.q_crit * c
        return D_prime

    def component_drag(self):
        ae = AerodynamicCharacteristics()
        ae.aero_functions(AoA_cruise=2)
        self.D_tank_sys = ae.C_D_0_tank_sys_HACK * self.q_crit * self.S
        self.D_eng = ae.C_D_o_engine * self.q_crit * self.S

    def body_weight(self, x):
        Wb_prime = self.WbS * self.chord(x=x)
        return Wb_prime

    def wing_strcut_weight(self, x):
        Ww_prime = self.WwS * self.chord(x=x)
        return Ww_prime

    def kerosene_distribution(self):
        self.kerosene, dummy, self.kerosene_max = kerosene_calc()
        #self.kerosene = function of kerosene weight distribution from o to kerosene_max
        #self.kerosene_max = half spanwise location

    def mc_step(self, dist, i):
        return max(0, dist-self.dx/2)**(i + 1) / (dist-self.dx/2)

    def chord_integrals(self):
        def cint1(x):
            return self.m1 * (x**2/2 - (0.5*self.width_f)**2/2) + self.c_root * (x - 0.5*self.width_f)

        def cint2(x):
            return self.m1 * (x**2/2 - (self.b_in/2)**2/2) + self.c_root * (x - self.b_in/2)

        def cint3(x):
            return self.m2 * (x**2/2 - (self.b_in/2)**2/2 - self.b_in/2 * (x - self.b_in/2)) + \
                   self.c_kink_out * (x - self.b_in/2)

        def c2int1(x):
            return self.m1 * (x**3/6 - self.width_f**2/8*x + self.width_f**3/16 - self.width_f**3/48) + \
                   self.c_root * (x**2/2 - self.width_f/2*x + self.width_f**2/4 - self.width_f**2/8)

        def c2int2(x):
            return self.m1 * (x**3/6 - self.b_in**2/8*x + self.b_in**3/16 - self.b_in**3/48) + \
                   self.c_root * (x**2/2 - self.b_in/2*x + self.b_in**2/4 - self.b_in**2/8)

        def c2int3(x):
            return self.m2 * (x**3/6 - self.b_in**2/8*x + self.b_in**3/16 - self.b_in**3/48) + \
                   (self.c_kink_out - self.m2*self.b_in/2) * (x**2/2 - self.b_in/2*x + self.b_in**2/4 - self.b_in**2/8)

        def c_sq_int1(x):
            return self.m1**2 * (x**3/3 - self.width_f**3/24) + self.m1 * self.c_root * (x**2 - self.width_f**2/4) \
                   + self.c_root**2 * (x - self.width_f/2)

        def c_sq_int2(x):
            return self.m1**2 * (x**3/3 - self.b_in**3/24) + self.m1 * self.c_root * (x**2 - self.b_in**2/4) \
                   + self.c_root**2 * (x - self.b_in/2)

        def c_sq_int3(x):
            return self.m2**2 * (x**3/3 - self.b_in**3/24) + \
                   self.m2 * (self.c_kink_out - self.m2 * self.b_in/2) * (x**2 - self.b_in**2/4) \
                   + (self.c_kink_out - self.m2 * self.b_in/2)**2 * (x - self.b_in/2)

        self.cint1 = cint1
        self.cint2 = cint2
        self.cint3 = cint3
        self.c2int1 = c2int1
        self.c2int2 = c2int2
        self.c2int3 = c2int3
        self.c_sq_int1 = c_sq_int1
        self.c_sq_int2 = c_sq_int2
        self.c_sq_int3 = c_sq_int3

    def Reaction_forces(self):
        self.component_drag()

        x_arr = np.arange(self.width_f/2, self.b/2 + self.dx, self.dx)
        L_arr = np.zeros(len(x_arr))
        D_arr = np.zeros(len(x_arr))
        Ww_arr = np.zeros(len(x_arr))
        Wb_arr = np.zeros(len(x_arr))
        c_arr = np.zeros(len(x_arr))
        for i, x in enumerate(x_arr):
            L_arr[i] = self.Lift(x=x)
            D_arr[i] = self.Drag_w(x=x)
            Ww_arr[i] = self.wing_strcut_weight(x=x)
            c_arr[i] = self.chord(x=x)

        self.RF_y = spint.simps(y=(L_arr-Ww_arr), x=x_arr) - self.W_engine * self.g_0 \
                   - self.pod_tank_mass * self.g_0
        self.RF_x = spint.simps(y=D_arr, x=x_arr) + self.D_tank_sys + self.D_eng - self.max_T

    def S_y(self, x):
        Sy = - self.RF_y * self.mc_step(dist=x-0.5*self.width_f, i=0) \
             - self.WwS * (+ self.cint1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                           - self.cint2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                           + self.cint3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)) \
             + self.C_L_crit * self.q_crit * (+ self.cint1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                                              - self.cint2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                                              + self.cint3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)) \
             - self.W_engine * self.g_0 * self.mc_step(dist=x-self.y_engine, i=0) \
             - self.pod_tank_mass * self.g_0 * self.mc_step(dist=x-self.y_cg_pod, i=0)
        return Sy

    def S_x(self, x):
        self.component_drag()
        Sx = - self.RF_x * self.mc_step(dist=x-0.5*self.width_f, i=0) \
             + self.C_D_crit * self.q_crit * (+ self.cint1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                                              - self.cint2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                                              + self.cint3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)) \
             + self.D_tank_sys * self.mc_step(dist=x-self.y_cg_pod, i=0) \
             + (self.D_eng - self.max_T) * self.mc_step(dist=x-self.y_engine, i=0)
        return Sx

    def Reaction_moments(self):
        x_arr = np.arange(self.width_f/2, self.b/2 + self.dx, self.dx)
        Sy_arr = np.zeros(len(x_arr))
        Sx_arr = np.zeros(len(x_arr))
        for i, x in enumerate(x_arr):
            Sy_arr[i] = self.S_y(x=x)
            Sx_arr[i] = self.S_x(x=x)

        self.RM_x = + spint.simps(y=Sy_arr, x=x_arr)
        self.RM_y = - spint.simps(y=Sx_arr, x=x_arr)

    def M_x(self, x):
        Mx = - self.RM_x * self.mc_step(dist=x-0.5*self.width_f, i=0) \
             - self.RF_y * self.mc_step(dist=x-0.5*self.width_f, i=1) \
             - self.WwS * (+ self.c2int1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                           - self.c2int2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                           + self.c2int3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)) \
             + self.C_L_crit * self.q_crit * (+ self.c2int1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                                              - self.c2int2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                                              + self.c2int3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)) \
             - self.W_engine * self.g_0 * self.mc_step(dist=x-self.y_engine, i=1) \
             - self.pod_tank_mass * self.g_0 * self.mc_step(dist=x-self.y_cg_pod, i=1)
        return Mx

    def M_y(self, x):
        My = - self.RM_y * self.mc_step(dist=x-0.5*self.width_f, i=0) \
             + self.RF_x * self.mc_step(dist=x-0.5*self.width_f, i=1) \
             - self.C_D_crit * self.q_crit * (+ self.c2int1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                                              - self.c2int2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                                              + self.c2int3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)) \
             - self.D_tank_sys * self.mc_step(dist=x-self.y_cg_pod, i=1) \
             - (self.D_eng - self.max_T) * self.mc_step(dist=x-self.y_engine, i=1)
        return My

    def Reaction_torque(self):
        x = self.b/2
        self.RT = self.C_L_crit * self.q_crit * self.xL_c * \
                     (self.c_sq_int1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                      - self.c_sq_int2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                      + self.c_sq_int3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)) \
                     + (self.max_T - self.D_eng) * self.T_arm * self.mc_step(dist=x-self.y_engine, i=0) \
                     - self.W_engine * self.W_eng_arm * self.mc_step(dist=x-self.y_engine, i=0) \
                     + self.D_tank_sys * (self.pylon_height + d_wing_pod/2 + 0.091/2 * self.chord(x=self.y_cg_pod)) * \
                     self.mc_step(dist=x-self.y_cg_pod, i=0) \
                     + self.WwS * self.W_wing_arm_c * (self.c_sq_int1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                                                       - self.c_sq_int2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                                                       + self.c_sq_int3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0))

    def T_z(self, x):
        Tz = - self.RT * self.mc_step(dist=x-self.width_f/2, i=0) + \
             self.C_L_crit * self.q_crit * self.xL_c * (self.c_sq_int1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                                                        - self.c_sq_int2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                                                        + self.c_sq_int3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)) \
            + (self.max_T - self.D_eng) * self.T_arm * self.mc_step(dist=x-self.y_engine, i=0) \
            - self.W_engine * self.W_eng_arm * self.mc_step(dist=x-self.y_engine, i=0) \
            + self.D_tank_sys * (self.pylon_height + d_wing_pod/2 + 0.091/2 * self.chord(x=self.y_cg_pod)) * \
            self.mc_step(dist=x-self.y_cg_pod, i=0) \
            + self.WwS * self.W_wing_arm_c * (self.c_sq_int1(x=x) * self.mc_step(dist=x-0.5*self.width_f, i=0)
                                              - self.c_sq_int2(x=x) * self.mc_step(dist=x-self.b_in/2, i=0)
                                              + self.c_sq_int3(x=x) * self.mc_step(dist=x-self.b_in/2, i=0))
        return Tz

    def plot_loads(self):
        x_arr = np.arange(self.width_f/2, self.b/2 + self.dx, self.dx)
        Sy_arr = np.zeros(len(x_arr))
        Sx_arr = np.zeros(len(x_arr))
        My_arr = np.zeros(len(x_arr))
        Mx_arr = np.zeros(len(x_arr))
        Tz_arr = np.zeros(len(x_arr))

        for i, x in enumerate(x_arr):
            Sy_arr[i] = self.S_y(x=x)
            Sx_arr[i] = self.S_x(x=x)
            My_arr[i] = self.M_y(x=x)
            Mx_arr[i] = self.M_x(x=x)
            Tz_arr[i] = self.T_z(x=x)

        fig = plt.figure(constrained_layout=True, figsize=(15, 15))
        gs = GridSpec(3, 2, figure=fig)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])
        ax5 = fig.add_subplot(gs[2, :])

        ax1.plot(x_arr, Sy_arr/1000)
        ax1.set_ylabel(r'$S_y(z)$ [$kN$]', size=15)
        ax1.set_title(r'Shear force in $y$', size=20)

        ax2.plot(x_arr, Sx_arr/1000)
        ax2.set_ylabel(r'$S_x(z)$ [$kN$]', size=15)
        ax2.set_title(r'Shear force in $x$', size=20)

        ax3.plot(x_arr, Mx_arr/1000)
        ax3.set_ylabel(r'$M_x(z)$ [$kN/m$]', size=15)
        ax3.set_xlabel(r'$z$ [$m$]', size=15)
        ax3.set_title(r'Bending moment around $x$', size=20)

        ax4.plot(x_arr, My_arr/1000)
        ax4.set_ylabel(r'$M_y(z)$ [$kN/m$]', size=15)
        ax4.set_xlabel(r'$z$ [$m$]', size=15)
        ax4.set_title(r'Bending moment around $y$', size=20)

        ax5.plot(x_arr, Tz_arr/1000)
        ax5.set_ylabel(r'$T(z)$ [$kN/m$]', size=15)
        ax5.set_xlabel(r'$z$ [$m$]', size=15)
        ax5.set_title(r'Torque around $z$', size=20)

        for i in [ax1, ax2, ax3, ax4, ax5]:
            i.tick_params(axis='both', which='major', labelsize=12)

        plt.show()

if __name__ == '__main__':

    lw = Loads_w()
    lw.compute_loads()
    lw.chord_integrals()
    print(lw.c_sq_int1(x=lw.width_f/2))
    print(lw.c_sq_int2(x=lw.b_in/2))
    print(lw.c_sq_int3(x=lw.b_in/2))
    lw.Reaction_forces()
    lw.Reaction_moments()
    lw.Reaction_torque()

    lw.plot_loads()

    x_arr = np.arange(lw.width_f/2, lw.b/2 + lw.dx, lw.dx)
    Sy_arr = np.zeros(len(x_arr))
    Sx_arr = np.zeros(len(x_arr))
    My_arr = np.zeros(len(x_arr))
    Mx_arr = np.zeros(len(x_arr))

    for i, x in enumerate(x_arr):
        Sy_arr[i] = lw.S_y(x=x)
        Sx_arr[i] = lw.S_x(x=x)
        My_arr[i] = lw.M_y(x=x)
        Mx_arr[i] = lw.M_x(x=x)

    # print('Max lift = ', lw.L_max)
    # print('From the plot it is = ', 2 * spint.simps(y=L_arr, x=x_arr))
    # print('The integrated surface is = ', 2 * spint.simps(y=c_arr, x=x_arr))

    plt.plot(x_arr, Sy_arr)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('y', size=15)
    plt.ylabel('N', size=15)
    plt.show()



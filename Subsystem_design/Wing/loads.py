from Subsystem_design.common_constants import Constants
from Subsystem_design.Performance.flight_envelope import FlightEnvelope
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spint
import pandas as pd

class Loads_w(Constants):

    def __init__(self):
        super().__init__()
        self.dx = 0.01
        self.m1 = (self.c_kink_out - self.c_root) / (0.5 * self.b_in)
        self.m2 = (self.c_tip - self.c_kink_out) / (0.5 * self.b_out)

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

        x_arr = np.linspace(0.5*lw.width_f, lw.b/2, 1000)
        c_arr = np.zeros(len(x_arr))
        for i, x in enumerate(x_arr):
            c_arr[i] = lw.chord(x=x)
        S_exp = 2 * spint.simps(y=c_arr, x=x_arr)  # Exposed surface ara [m^2]
        self.C_L_crit = self.L_max / (self.q_crit * S_exp)
        self.C_D_crit = ae.C_D_0_HACK - ae.C_D_0_fus_neo - ae.C_D_0_tank_sys_HACK - ae.C_D_o_engine - ae.C_D_0_Vtail - \
                        ae.C_D_0_Htail + self.C_L_crit / (np.pi * ae.AR * self.e)


    def Lift(self, x):
        c = self.chord(x=x)
        L_prime = 1.5 * self.C_L_crit * self.q_crit * c  # With the safety factor of 1.5
        return L_prime

    def Drag_w(self, x):
        c = self.chord(x=x)
        D_prime = self.C_D_crit * self.q_crit * c

    def component_drag(self):
        ae = AerodynamicCharacteristics()
        ae.aero_functions(AoA_cruise=2)
        self.D_tank_sys = ae.C_D_0_tank_sys_HACK * self.q_crit * self.S
        self.D_eng = ae.C_D_o_engine * self.q_crit * self.S

    def mc_step(self, dist, i):
        return max(0, dist-self.dx/2)**(i + 1) / (dist-self.dx/2)

    def S_y(self, x):
        Sy = self.C_L_crit * self.q_crit * \
             (+ (self.m1 * (x**2/2 - (0.5*self.width_f)**2/2) + self.c_root * (x - 0.5*self.width_f)) *
              self.mc_step(dist=x-0.5*self.width_f, i=0)
              - (self.m1 * (x**2/2 - (self.b_in/2)**2/2) + self.c_root * (x - self.b_in/2)) *
              self.mc_step(dist=x-self.b_in/2, i=0)
              + (self.m2 * (x**2/2 - (self.b_in/2)**2/2 - self.b_in/2 * (x - self.b_in/2)) + self.c_kink_out *
                 (x - self.b_in/2)) * self.mc_step(dist=x-self.b_in/2, i=0))
        return Sy

    def S_x(self, x):
        self.component_drag()
        Sx = self.C_D_crit * self.q_crit * \
             (+ (self.m1 * (x**2/2 - (0.5*self.width_f)**2/2) + self.c_root * (x - 0.5*self.width_f)) *
              self.mc_step(dist=x-0.5*self.width_f, i=0)
              - (self.m1 * (x**2/2 - (self.b_in/2)**2/2) + self.c_root * (x - self.b_in/2)) *
              self.mc_step(dist=x-self.b_in/2, i=0)
              + (self.m2 * (x**2/2 - (self.b_in/2)**2/2 - self.b_in/2 * (x - self.b_in/2)) + self.c_kink_out *
                 (x - self.b_in/2)) * self.mc_step(dist=x-self.b_in/2, i=0)) \
             + self.D_tank_sys * self.mc_step(dist=x-) + self.D_eng * self.mc_step(dist=x-)






if __name__ == '__main__':

    lw = Loads_w()
    lw.compute_loads()

    x_arr = np.linspace(0.5*lw.width_f, lw.b/2, 500)
    L_arr = np.zeros(len(x_arr))
    Sy_arr = np.zeros(len(x_arr))
    c_arr = np.zeros(len(x_arr))
    for i, x in enumerate(x_arr):
        L_arr[i] = lw.Lift(x=x)
        Sy_arr[i] = lw.S_y(x=x)
        c_arr[i] = lw.chord(x=x)


    print('Max lift = ', lw.L_max)
    print('From the plot it is = ', 2 * spint.simps(y=L_arr, x=x_arr))
    print('The integrated surface is = ', 2 * spint.simps(y=c_arr, x=x_arr))

    plt.plot(x_arr, Sy_arr)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('y', size=15)
    plt.ylabel('N', size=15)
    plt.show()



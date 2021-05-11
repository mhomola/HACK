import numpy as np
import math as m
import matplotlib.pyplot as plt
import skimage.io as io
import sys



class Roskam_drag_coefficient():
    def __init__(self, visc, u1, air_d, l_f, M, S):
        '''
        :param visc: coefficient of viscosity of air
        :param u1: steady state airspeed
        :param air_d: density of air
        :param l_f: length of the fuselage
        :param M: Mach number
        :param S: Wing surface area
        '''
        self.visc = visc
        self.u1 = u1
        self.air_d = air_d
        self.l_f = l_f
        self.M = M
        self.S = S

    def run_Roskam_drag_coefficient_functions(self, l_cockpit, l_cabin, l_tail, S_fus, S_b_fus, S_wet_fus, ):
        self.wing_fus_interference()
        self.turbulent_flat_plate_skin_friction()
        self.zero_lift_drag_fus(S_fus, S_b_fus, S_wet_fus)

        return self.R_wf, self.C_f_fus, self.C_D_o_fus

    '''
    def exposed_mgc(self, ct, cr, sweep, d_fus, ):
        """
        This function calculates the mean geometrical chord of the exposed wing.
        :param ct: tip chord in [m]
        :param cr: root chord in [m]
        :param sweep: sweep angle in [deg]
        :param d_fus: fuselage diameter in [m]
        :return: exp_mgc, the mean geometric chord in [m]
        """
        delta_cr = d_fus/2/m.tan(m.radians(90-sweep)) # difference between root chord and exposed root chord
        exposed_cr = cr + delta_cr
        exposed_taper = ct/exposed_cr

        #Formula for mean geometric chord obtained from ADSEEII course
        self.exp_mgc = 2/3 * exposed_cr * (1 + exposed_taper + exposed_taper**2)/(1+exposed_taper)
    '''

    def wing_fus_interference(self, ):
        '''
        :return: R_wf, the wing/fuselage iterference factor as derived from Figure 4.1 in Roskam-VI
        '''
        self.R_n_fus = self.air_d * self.u1 * self.l_f / self.visc
        print('The Fuselage Reynolds Number R_f_fus is: ', self.R_n_fus, ' [-]')
        print('The Mach number M is: ', self.M)
        print('Next you will get a figure from which you can obtain the wing/fuselage iterference factor R_wf \n '
              'as a function of R_n_fus and M. \n'
              'Find this value and once the tab closes you can input the value.')

        # fus_RE = io.imread('fuselage_Re.PNG')
        # plt.imshow(fus_RE)
        # plt.show()
        self.R_wf = 1.015

    def turbulent_flat_plate_skin_friction(self, ):
        '''
        :return: C_f, the turbulent flat plate skin friction coefficient as derived from Figure 4.3 in Roska-VI
        '''
        print('The Fuselage Reynolds Number R_f_fus is: ', self.R_n_fus, ' [-]')
        print('The Mach number M is: ', self.M)
        print('Next you will get a figure from which you can obtain the turbulent flat plate skin \n '
              'friction coefficient C_f as a function of R_n_fus and M. \n'
              'Find this value and once the tab closes you can input the value.')

        # flat_plate = io.imread('flat_plate_skin_friction.PNG')
        # plt.imshow(flat_plate)
        # plt.show()
        self.C_f_fus = 0.002

    def zero_lift_drag_fus(self, S_b_fus, S_fus, S_wet_fus):
        '''
        :param S_b_fus: Base area (at the end of the fuselage)
        :param S_fus: Largest cross sectional area of the fuselage [m^2]
        :return: C_D_o_fus, the subsonic fuselage zero-lift drag coefficient
        '''
        self.d_f = np.sqrt(4 / np.pi * S_fus)
        ld = self.l_f / self.d_f
        self.C_D_o_fus_exc_base = self.R_wf * self.C_f_fus * ((1 + 60 / ld**3) + 0.0025 * ld) * S_wet_fus / self.S
        bf = np.sqrt(4 / np.pi * S_b_fus) / self.d_f
        self.C_D_b_fus = (0.029 * bf**3 / ((self.C_D_o_fus_exc_base * self.S / S_fus)**0.5)) * S_fus / self.S
        self.C_D_o_fus = self.C_D_o_fus_exc_base + self.C_D_b_fus













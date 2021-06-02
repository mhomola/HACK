"""
This files contains the classes and functions for determining the insulation properties of the tank.
"""
import numpy as np


T_h = 273.15 + 45
T_c = 19.75
DT = 0.25

r = 1.14/2

Cs = 7.3 * 10**-8
Cr = 7.07 * 10**-10
Cg = 1.46 * 10*8-4

epsl = 0.031

N_layers = 40*3         # number of layers
N_density = 22.7        #
P = 1.3*10**-3              # Pa


def heat_flux():

    q_tot = (Cr*epsl*(T_h**4.67-T_c**4.67)/N_layers) + (Cs*N_density**2.63*(T_h-T_c)*(T_h+T_c)/(2*N_layers+2)) + \
            (Cg*P*(T_h**0.52-T_c**0.52)/N_layers)
    Q_tot = (2*np.pi*r*(6.5-r*2) + 4*np.pi*r**2)*q_tot

    return q_tot, Q_tot


q_tank, Q_tank = heat_flux()

print(Q_tank)
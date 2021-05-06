import numpy as np

import matplotlib.pyplot as plt

import math

MAC = 4.29

total_length = 28.55



S = 77.3 #m^2

b = 26.21

AR = 8.89 * 1.15

taper_ratio = 0.356

qc_sweep = 15

x_wing = 10.34 #mat



"""

Preliminary calculations

"""

tan_le_sweep = math.tan(qc_sweep * math.pi/180)-(4/8.89)*(0-0.25)*((1-taper_ratio)/(1+taper_ratio))

root_chord = 2*S/((1+taper_ratio)*b)

tip_chord = taper_ratio*root_chord

y_lemac = (b/6)*((1+2*taper_ratio)/(1+taper_ratio))

x_lemac = y_lemac * tan_le_sweep



X_MAC = x_lemac + x_wing #leading edge

cg_OEW = X_MAC + MAC*0.25 - 0.6
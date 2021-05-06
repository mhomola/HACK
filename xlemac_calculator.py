import numpy as np

import matplotlib.pyplot as plt

import math

MAC = 4.29

total_length = 37.57



S = 122.4 #m^2

b = 33.91

AR = 9.39

taper_ratio = 0.24

qc_sweep = 25

x_wing = 11.65 #mat

# as a fraction of MAC
fwd_cg_limit = 0.17
aft_cg_limit = 0.368
cg_OEW_frac = 0.265



"""

Preliminary calculations

"""

tan_le_sweep = math.tan(qc_sweep * math.pi/180)-(4/9.39)*(0-0.25)*((1-taper_ratio)/(1+taper_ratio))

root_chord = 2*S/((1+taper_ratio)*b)

tip_chord = taper_ratio*root_chord

y_lemac = (b/6)*((1+2*taper_ratio)/(1+taper_ratio))

x_lemac = y_lemac * tan_le_sweep


X_MAC = x_lemac + x_wing #leading edge

print(X_MAC)

fwd_cg = X_MAC + fwd_cg_limit*MAC
aft_cg = X_MAC + aft_cg_limit*MAC
cg_OEW = X_MAC + cg_OEW_frac*MAC

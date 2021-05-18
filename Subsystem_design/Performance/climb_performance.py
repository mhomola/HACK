import fuel_constants
from Subsystem_design.fuel_required import fuel_volume_calc
from Subsystem_design.common_constants import *
from Subsystem_design.common_constants import Constants
import numpy as np

max_cont_thrust = 9.02*10**6
takeoff_thrust = 9.17*10**6
dvdh = V_EAS*  #V_EAS * d/dH(

def rc_steady(C_D, C_L, V, W):
    rc_s = (max_cont_thrust/W) - (C_D*V/C_L)
    return rc_s

def rc_real(rc_steady, V_EAS):
    rc = rc_steady/(1+V*dvdh/9.8)
    return



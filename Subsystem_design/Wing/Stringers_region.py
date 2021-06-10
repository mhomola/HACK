from Subsystem_design.Wing.inertia import Inertia
from Subsystem_design.common_constants import Constants
from Subsystem_design.Wing.loads import Loads_w
from Subsystem_design.Wing.Stresses import stresses
import numpy as np
import matplotlib.pyplot as plt
import math
###Initializing Loads Class
lw = Loads_w() #wing loads object
lw.compute_loads()
lw.chord_integrals()
lw.Reaction_forces()
lw.Reaction_moments()
lw.Reaction_torque()

def no_str(lw, x):
    """
    Functions for determining the number of stringers for supporting the yield strength.
    :return:
    """
    max_yield = 480  # [MPa]
    n_str = 0 # we start with 0 stringers

    MOI = Inertia(n_str=n_str)  # we start with 1 stringer
    MOI.compute_inertia(x=x)
    c = lw.chord(x)

    max_str = math.ceil(
        MOI.w_sk_c * c / MOI.h_str)  # check how many stringers can be fitted depending on their width and
    # width of the wingbox at that location

    wing_stress = stresses(Ixx=MOI.Ixx_no_str, Iyy=MOI.Iyy_no_str, Ixx_str=MOI.Ixx, Iyy_str=MOI.Iyy,
                  h=MOI.h_sp_c * c, L=MOI.w_sk_c * c, t_upper=MOI.t_sk, t_spar1=MOI.t_sp, t_spar2=MOI.t_sp,
                  t_lower=MOI.t_sk)

    wing_stress.shear_loads(Vx=lw.S_x(x=x), Vy=lw.S_y(x=x), T=lw.T_z(x=x))
    wing_stress.bending_loads(Mx=lw.M_x(x=x), My=lw.M_y(x=x))
    wing_stress.compute_stresses()
    wing_stress.vm_plotter(show=False)

    while wing_stress.vm_max / 10 ** 6 >= max_yield and n_str<max_str:
        n_str = n_str + 1
        MOI = Inertia(n_str=n_str)
        MOI.compute_inertia(x=x)
        wing_stress = stresses(Ixx=MOI.Ixx_no_str, Iyy=MOI.Iyy_no_str, Ixx_str=MOI.Ixx, Iyy_str=MOI.Iyy,
                      h=MOI.h_sp_c * c, L=MOI.w_sk_c * c, t_upper=MOI.t_sk, t_spar1=MOI.t_sp, t_spar2=MOI.t_sp,
                      t_lower=MOI.t_sk)
        wing_stress.shear_loads(Vx=lw.S_x(x=x), Vy=lw.S_y(x=x), T=lw.T_z(x=x))
        wing_stress.bending_loads(Mx=lw.M_x(x=x), My=lw.M_y(x=x))
        wing_stress.compute_stresses()
        wing_stress.vm_plotter(show=False)

    return n_str,wing_stress.vm_max/10**6

x_arr = np.arange(2.05, lw.b/2 + lw.dx,1)

for i, x in enumerate(x_arr):
    print(x,no_str(lw,x))
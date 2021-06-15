from Subsystem_design.Wing.inertia import Inertia_normal, Inertia_shear
from Subsystem_design.common_constants import Constants
from Subsystem_design.Wing.loads import Loads_w
from Subsystem_design.Wing.Stresses import stresses
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math

SF= 1.5 #safety factor
###Initializing Loads Class
lw = Loads_w() #wing loads object

lw.compute_loads()
lw.chord_integrals()
lw.Reaction_forces()
lw.Reaction_moments()
lw.Reaction_torque()

x_arr = np.arange(lw.width_f/2, lw.b/2 + lw.dx, lw.dx)
Sy_arr = np.zeros(len(x_arr))
Sx_arr = np.zeros(len(x_arr))
My_arr = np.zeros(len(x_arr))
Mx_arr = np.zeros(len(x_arr))
T_arr = np.zeros(len(x_arr))
vm_arr = np.zeros(len(x_arr))
sigma_arr = np.zeros(len(x_arr))
compression_arr = np.zeros(len(x_arr))
tension_arr = np.zeros(len(x_arr))
shear_stress_arr = np.zeros(len(x_arr))
vm = np.zeros((len(x_arr),1000))
q_arr = np.zeros((len(x_arr),1000))

Ixx_arr = np.zeros(len(x_arr))
Iyy_arr = np.zeros(len(x_arr))
Ixx_nostr_arr = np.zeros(len(x_arr))
Iyy_nostr_arr = np.zeros(len(x_arr))
### Initializing Inertia class
#MOI = Inertia(n_str = 10)
MOI_normal = Inertia_normal(n_str = 11)
MOI_shear = Inertia_shear(n_str = 11)
##We want to compute the inertia, the loads and the stresses at every span-wise location

for i, x in enumerate(x_arr):
    c = lw.chord(x=x) #chord dimension for scaling

    Sy_arr[i] = lw.S_y(x=x) * SF
    Sx_arr[i] = lw.S_x(x=x) * SF
    My_arr[i] = lw.M_y(x=x) * SF
    Mx_arr[i] = lw.M_x(x=x) * SF
    T_arr[i] = lw.T_z(x=x)  * SF
    MOI_normal.compute_inertia(x=x)
    MOI_shear.compute_inertia(x=x)



    ###we need to initialize the stresses class at every spanwise location

    ### NORMAL STRESSES
    wing_stress = stresses(Ixx=MOI_shear.Ixx_shear, Iyy=MOI_shear.Iyy_shear, Ixx_str=MOI_normal.Ixx_normal, Iyy_str=MOI_normal.Iyy_normal,
                           h=MOI_normal.h_sp_c*c, L=MOI_normal.w_sk_c*c, t_upper=MOI_normal.t_sk, t_spar1=MOI_normal.t_sp, t_spar2=MOI_normal.t_sp, t_lower=MOI_normal.t_sk)

    wing_stress.shear_loads(Vx=Sx_arr[i],Vy=Sy_arr[i],T=T_arr[i])
    wing_stress.bending_loads(Mx=Mx_arr[i],My=My_arr[i])
    wing_stress.compute_stresses()
    wing_stress.vm_plotter(show=False)
    wing_stress.sigma_plotter(show=False)
    sigma_arr[i] = wing_stress.sigma_max
    compression_arr[i] = wing_stress.sigma_compressive
    tension_arr[i] = wing_stress.sigma_tension
    sigma_arr[i] = wing_stress.sigma_max
    vm_arr[i] = wing_stress.vm_max

    ###SHEAR STRESSES
    wing_stress = stresses(Ixx=MOI_shear.Ixx_shear, Iyy=MOI_shear.Iyy_shear, Ixx_str=MOI_normal.Ixx_normal, Iyy_str=MOI_normal.Iyy_normal,
                           h=MOI_shear.h_sp_c * c, L=MOI_shear.w_sk_c * c, t_upper=MOI_shear.t_sk, t_spar1=MOI_shear.t_sp, t_spar2=MOI_shear.t_sp, t_lower=MOI_shear.t_sk)

    wing_stress.shear_loads(Vx=Sx_arr[i], Vy=Sy_arr[i], T=T_arr[i])
    wing_stress.shear_loads(Vx=Sx_arr[i], Vy=Sy_arr[i], T=T_arr[i])
    wing_stress.bending_loads(Mx=Mx_arr[i], My=My_arr[i])
    wing_stress.compute_stresses()
    wing_stress.shear_flow_plotter(type="total", show=False)
    q_arr[i] = wing_stress.q_max
    shear_stress_arr[i] = wing_stress.shear_stress_max



#lw.plot_loads()
#plt.figure()

###VON MISES STRESS PLOT
# plt.plot(x_arr,vm_arr/10**6)
# plt.xlabel("Span location[m]")
# plt.ylabel("Von Mises Stress [MPa]")
# plt.figure()

##COMPRESSION STRESS PLOT
print("Compression stress",min(compression_arr)/10**6)
plt.plot(x_arr,compression_arr/10**6)
plt.xlabel("Span location[m]")
plt.ylabel("Compression Stress [MPa]")
plt.figure()

###TENSION STRESS PLOT
# plt.plot(x_arr,tension_arr/10**6)
# plt.xlabel("Span location[m]")
# plt.ylabel("Tension Stress [MPa]")
# plt.figure()

###SHEAR FLOWS PLOT
# plt.plot(x_arr,q_arr/10**6)
# plt.xlabel("Span location[m]")
# plt.ylabel("Shear flows [MPa/m]")

###SHEAR STRESS PLOT
print("Max shear stress",max(shear_stress_arr)/10**6)
plt.plot(x_arr,shear_stress_arr/10**6)
plt.xlabel("Span location[m]")
plt.ylabel("Shear Stress [MPa]")
plt.show()

###MOI Graphs
# plt.plot(x_arr,Ixx_arr,"blue")
# plt.plot(x_arr,Iyy_arr,"red")
# plt.title("MOI with stringers")
# plt.figure()
# plt.plot(x_arr, Ixx_nostr_arr, "blue")
# plt.plot(x_arr, Iyy_nostr_arr, "red")
# plt.title("MOI with no stringers")

# sns.heatmap(vm.transpose(),cmap="magma",yticklabels=False,xticklabels=False) #cividis
# plt.show()


#Plot distribution at one span wise point

# ind = 180
# c = lw.chord(x = x_arr[ind])
# wing_stress_plot = stresses(Ixx=MOI.Ixx_no_str, Iyy=MOI.Iyy_no_str, Ixx_str=MOI.Ixx, Iyy_str=MOI.Iyy,
#                        h=MOI.h_sp_c * c, L=MOI.w_sk_c * c, t_upper=MOI.t_sk, t_spar1=MOI.t_sp, t_spar2=MOI.t_sp,
#                        t_lower=MOI.t_sk)
#
# wing_stress_plot.shear_loads(Vx=Sx_arr[ind], Vy=Sy_arr[ind], T=T_arr[ind])
# wing_stress_plot.bending_loads(Mx=Mx_arr[ind], My=My_arr[ind])
# wing_stress_plot.compute_stresses()
# print("Mx=",Mx_arr[ind],"My=",My_arr[ind])
# print("Vx=",Sx_arr[ind],"Vy=",Sy_arr[ind],"T=",T_arr[ind])
# wing_stress_plot.shear_flow_plotter(type = "total",show=True)
# wing_stress_plot.vm_plotter(show=True)

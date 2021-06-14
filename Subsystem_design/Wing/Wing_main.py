from Subsystem_design.Wing.inertia import Inertia, Inertia_initial
from Subsystem_design.common_constants import Constants
from Subsystem_design.Wing.loads import Loads_w
from Subsystem_design.Wing.Stresses import stresses
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math
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
vm = np.zeros((len(x_arr),1000))

Ixx_arr = np.zeros(len(x_arr))
Iyy_arr = np.zeros(len(x_arr))
Ixx_nostr_arr = np.zeros(len(x_arr))
Iyy_nostr_arr = np.zeros(len(x_arr))
### Initializing Inertia class
#MOI = Inertia(n_str = 10)
MOI = Inertia_initial(n_str = 10 )
##We want to compute the inertia, the loads and the stresses at every span-wise location

for i, x in enumerate(x_arr):
    c = lw.chord(x=x) #chord dimension for scaling

    Sy_arr[i] = lw.S_y(x=x)
    Sx_arr[i] = lw.S_x(x=x)
    My_arr[i] = lw.M_y(x=x)
    Mx_arr[i] = lw.M_x(x=x)
    T_arr[i] = lw.T_z(x=x)
    MOI.compute_inertia(x=x)



    ###we need to initialize the stresses class at every spanwise location
    wing_stress = stresses(Ixx=MOI.Ixx_no_str, Iyy=MOI.Iyy_no_str, Ixx_str=MOI.Ixx, Iyy_str=MOI.Iyy,
                           h=MOI.h_sp_c*c, L=MOI.w_sk_c*c, t_upper=MOI.t_sk, t_spar1=MOI.t_sp, t_spar2=MOI.t_sp, t_lower=MOI.t_sk)

    wing_stress.shear_loads(Vx=Sx_arr[i],Vy=Sy_arr[i],T=T_arr[i])
    wing_stress.bending_loads(Mx=Mx_arr[i],My=My_arr[i])
    wing_stress.compute_stresses()
    wing_stress.vm_plotter(show=False)
    wing_stress.sigma_plotter(show=False)
    vm_arr[i] = wing_stress.vm_max
    vm[i,:] = wing_stress.vm3

    sigma_arr[i] = wing_stress.sigma_max

    ###MOI verification
    Ixx_arr[i] = MOI.Ixx
    Iyy_arr[i] = MOI.Iyy
    Ixx_nostr_arr[i] = MOI.Ixx_no_str
    Iyy_nostr_arr[i] = MOI.Iyy_no_str

#lw.plot_loads()
#plt.figure()
plt.plot(x_arr,vm_arr/10**6)
plt.xlabel("Span location[m]")
plt.ylabel("Von Mises Stress [MPa]")
plt.figure()
plt.plot(x_arr,sigma_arr/10**6)
plt.xlabel("Span location[m]")
plt.ylabel("Bneding stress [MPa]")
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

ind = 180
c = lw.chord(x = x_arr[ind])
wing_stress_plot = stresses(Ixx=MOI.Ixx_no_str, Iyy=MOI.Iyy_no_str, Ixx_str=MOI.Ixx, Iyy_str=MOI.Iyy,
                       h=MOI.h_sp_c * c, L=MOI.w_sk_c * c, t_upper=MOI.t_sk, t_spar1=MOI.t_sp, t_spar2=MOI.t_sp,
                       t_lower=MOI.t_sk)

wing_stress_plot.shear_loads(Vx=Sx_arr[ind], Vy=Sy_arr[ind], T=T_arr[ind])
wing_stress_plot.bending_loads(Mx=Mx_arr[ind], My=My_arr[ind])
wing_stress_plot.compute_stresses()
print("Mx=",Mx_arr[ind],"My=",My_arr[ind])
print("Vx=",Sx_arr[ind],"Vy=",Sy_arr[ind],"T=",T_arr[ind])
wing_stress_plot.shear_flow_plotter(type = "total",show=True)
wing_stress_plot.vm_plotter(show=True)

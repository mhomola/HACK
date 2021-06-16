from Subsystem_design.Wing.Inertia_variable import Inertia_normal_var,Inertia_shear_var
from Subsystem_design.Wing.inertia import Inertia_normal,Inertia_shear
from Subsystem_design.common_constants import Constants
from Subsystem_design.Wing.loads import Loads_w
from Subsystem_design.Wing.Stresses import stresses
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math

SF = 1.5

#We will only look at one span-wise location for this analysis
##CHANGE THIS AS DESIRED
x = 3 #[m]
##SELECT A POINT OF REFERENCE
s1 = 0.05

# type_s = "Max"
type_s = "Point"
#We will look for a variation from -50 to 50
percentages = np.arange(-0.5,0.5,0.1)

#Define first the reference data

class Reference_Data():
    def __init__(self,x,s1):

        ###INITIALIZE LOADS
        SF = 1.5
        lw = Loads_w()
        lw.compute_loads()
        lw.chord_integrals()
        lw.Reaction_forces()
        lw.Reaction_moments()
        lw.Reaction_torque()

        ###INITIALIZE MOMENT OF INERTIA
        MOI_normal = Inertia_normal(n_str=11)
        MOI_shear = Inertia_shear(n_str=11)
        MOI_normal.compute_inertia(x=x)
        MOI_shear.compute_inertia(x=x)


        Vx = lw.S_x(x=x) * SF
        Vy = lw.S_y(x=x) * SF
        My = lw.M_y(x=x) * SF
        Mx = lw.M_x(x=x) * SF
        T  = lw.T_z(x=x) * SF

        ### NORMAL STRESSES
        wing_stress = stresses(Ixx=MOI_shear.Ixx_shear, Iyy=MOI_shear.Iyy_shear, Ixx_str=MOI_normal.Ixx_normal,
                               Iyy_str=MOI_normal.Iyy_normal,
                               h=MOI_normal.h_sp_c * MOI_normal.chord_inertia(x=x), L=MOI_normal.w_sk, t_upper=MOI_normal.t_sk,
                               t_spar1=MOI_normal.t_sp, t_spar2=MOI_normal.t_sp, t_lower=MOI_normal.t_sk)

        wing_stress.shear_loads(Vx=Vx, Vy=Vy, T=T)
        wing_stress.bending_loads(Mx=Mx, My=My)
        wing_stress.compute_stresses()
        wing_stress.sigma_plotter(show=False)
        max_sigma = wing_stress.sigma_max
        point_sigma = wing_stress.sigma1

        ###SHEAR STRESSES
        wing_stress = stresses(Ixx=MOI_shear.Ixx_shear, Iyy=MOI_shear.Iyy_shear, Ixx_str=MOI_normal.Ixx_normal,
                               Iyy_str=MOI_normal.Iyy_normal,
                               h=MOI_shear.h_sp_c * MOI_shear.chord_inertia(x=x), L=MOI_shear.w_sk, t_upper=MOI_shear.t_sk,
                               t_spar1=MOI_shear.t_sp, t_spar2=MOI_shear.t_sp, t_lower=MOI_shear.t_sk)

        wing_stress.shear_loads(Vx=Vx, Vy=Vy, T=T)
        wing_stress.bending_loads(Mx=Mx, My=My)
        wing_stress.compute_stresses()
        wing_stress.shear_flow_plotter(type="total", show=False)
        max_shear = wing_stress.shear_stress_max
        point_shear = wing_stress.q1_tot

        ###PARAMETERS NORMAL CROSS SECTION
        self.t_str_normal = MOI_normal.t_str
        self.t_sp_normal = MOI_normal.t_sp
        self.t_sk_normal = MOI_normal.t_sk
        self.h_str_normal = MOI_normal.h_str
        self.w_str_normal = MOI_normal.w_str
        self.n_str_normal = MOI_normal.n_str
        self.w_sk_normal = MOI_normal.w_sk
        self.h_sp_c_normal = MOI_normal.h_sp_c

        ###PARAMETER SHEAR CROSS SECTION
        self.t_sp_shear = MOI_shear.t_sp
        self.t_sk_shear = MOI_shear.t_sk
        self.w_sk_shear = MOI_shear.w_sk
        self.h_sp_c_shear = MOI_shear.h_sp_c

        ###STORING THE LOADS
        self.Vx = Vx
        self.Vy = Vy
        self.My = My
        self.Mx = Mx
        self.T  = T

        ###STORING THE STRESSES
        if type_s == "Max":
            self.max_sigma = max_sigma
            self.max_shear = max_shear
        if type_s == "Point":
            self.max_sigma = point_sigma(s1)
            self.max_shear = point_shear(s1)/self.t_sk_shear

        ###MOI elements for when only the loads are changed
        self.Ixx_shear = MOI_shear.Ixx_shear
        self.Ixx_normal = MOI_normal.Ixx_normal
        self.Iyy_shear = MOI_shear.Iyy_shear
        self.Iyy_normal = MOI_normal.Iyy_normal

        ###OTHER PARAMETERS
        self.c = MOI_shear.chord_inertia(x=x)

ref = Reference_Data(x=x,s1=s1)

###BENDING STRESS VARIATION

##Creating empty arrays for storing the results
sensi_bending_Mx = np.zeros(len(percentages))
sensi_bending_My = np.zeros(len(percentages))
sensi_bending_t_str = np.zeros(len(percentages))
sensi_bending_t_sp = np.zeros(len(percentages))
sensi_bending_t_sk = np.zeros(len(percentages))
sensi_bending_h_str = np.zeros(len(percentages))

###Changing moment about x-axis Mx
for i,p in enumerate(percentages):
    if ref.Mx > 0:
        Mx = (1+p) * ref.Mx
    else:
        Mx = (1 - p) * ref.Mx

    ws = stresses(Ixx=ref.Ixx_shear, Iyy=ref.Iyy_shear, Ixx_str=ref.Ixx_normal, Iyy_str=ref.Iyy_normal,
                h=ref.h_sp_c_normal*ref.c, L=ref.w_sk_normal, t_upper=ref.t_sk_normal,
                t_spar1=ref.t_sp_normal, t_spar2=ref.t_sp_normal, t_lower=ref.t_sk_normal)
    ws.bending_loads(Mx = Mx,My = ref.My)
    ws.shear_loads(Vx=ref.Vx,Vy=ref.Vy,T=ref.T)
    ws.compute_stresses()
    ws.sigma_plotter(show=False)
    if type_s == "Max":
        sensi_bending_Mx[i] = ws.sigma_max
    if type_s == "Point":
        sensi_bending_Mx[i] = ws.sigma1(s1)

###Changing moment about y-axis My
for i,p in enumerate(percentages):
    if ref.My>0:
        My = (1+p) * ref.My
    else:
        My = (1 - p) * ref.My
    ws = stresses(Ixx=ref.Ixx_shear, Iyy=ref.Iyy_shear, Ixx_str=ref.Ixx_normal, Iyy_str=ref.Iyy_normal,
                h=ref.h_sp_c_normal*ref.c, L=ref.w_sk_normal, t_upper=ref.t_sk_normal,
                t_spar1=ref.t_sp_normal, t_spar2=ref.t_sp_normal, t_lower=ref.t_sk_normal)

    ws.bending_loads(Mx = ref.Mx,My = My)
    ws.shear_loads(Vx=ref.Vx,Vy=ref.Vy,T=ref.T)
    ws.compute_stresses()
    ws.sigma_plotter(show=False)

    if type_s == "Max":
        sensi_bending_My[i] = ws.sigma_max
    if type_s == "Point":
        sensi_bending_My[i] = ws.sigma1(s1)

###Changing the thickness of the stringers
for i,p in enumerate(percentages):
    t_str = (1+p)*ref.t_str_normal
    MOI = Inertia_normal_var(n_str=ref.n_str_normal,t_str=t_str,t_sp=ref.t_sp_normal,
                             t_sk=ref.t_sk_normal,h_str=ref.h_str_normal,w_str=ref.w_str_normal)
    MOI.compute_inertia(x=x)

    ws = stresses(Ixx=ref.Ixx_shear, Iyy=ref.Iyy_shear, Ixx_str=MOI.Ixx_normal, Iyy_str=MOI.Iyy_normal,
                  h=ref.h_sp_c_normal * ref.c, L=ref.w_sk_normal, t_upper=ref.t_sk_normal,
                  t_spar1=ref.t_sp_normal, t_spar2=ref.t_sp_normal, t_lower=ref.t_sk_normal)
    ws.bending_loads(Mx=ref.Mx, My=ref.My)
    ws.shear_loads(Vx=ref.Vx, Vy=ref.Vy, T=ref.T)
    ws.compute_stresses()
    ws.sigma_plotter(show=False)

    if type_s == "Max":
        sensi_bending_t_str[i] = ws.sigma_max
    if type_s == "Point":
        sensi_bending_t_str[i] = ws.sigma1(s1)

###Changing the thickness of the spar
for i,p in enumerate(percentages):
    t_sp = (1+p)*ref.t_sp_normal
    MOI = Inertia_normal_var(n_str=ref.n_str_normal,t_str=ref.t_str_normal,t_sp=t_sp,
                             t_sk=ref.t_sk_normal,h_str=ref.h_str_normal,w_str=ref.w_str_normal)
    MOI.compute_inertia(x=x)

    ws = stresses(Ixx=ref.Ixx_shear, Iyy=ref.Iyy_shear, Ixx_str=MOI.Ixx_normal, Iyy_str=MOI.Iyy_normal,
                  h=ref.h_sp_c_normal * ref.c, L=ref.w_sk_normal, t_upper=ref.t_sk_normal,
                  t_spar1=t_sp, t_spar2=t_sp, t_lower=ref.t_sk_normal)
    ws.bending_loads(Mx=ref.Mx, My=ref.My)
    ws.shear_loads(Vx=ref.Vx, Vy=ref.Vy, T=ref.T)
    ws.compute_stresses()
    ws.sigma_plotter(show=False)

    if type_s == "Max":
        sensi_bending_t_sp[i] = ws.sigma_max
    if type_s == "Point":
        sensi_bending_t_sp[i] = ws.sigma1(s1)

###Changing the thickness of the skin
for i,p in enumerate(percentages):
    t_sk = (1+p)*ref.t_sk_normal
    MOI = Inertia_normal_var(n_str=ref.n_str_normal,t_str=ref.t_str_normal,t_sp=ref.t_sp_normal,
                             t_sk=t_sk,h_str=ref.h_str_normal,w_str=ref.w_str_normal)
    MOI.compute_inertia(x=x)

    ws = stresses(Ixx=ref.Ixx_shear, Iyy=ref.Iyy_shear, Ixx_str=MOI.Ixx_normal, Iyy_str=MOI.Iyy_normal,
                  h=ref.h_sp_c_normal * ref.c, L=ref.w_sk_normal, t_upper=t_sk,
                  t_spar1=ref.t_sp_normal, t_spar2=ref.t_sp_normal, t_lower=t_sk)
    ws.bending_loads(Mx=ref.Mx, My=ref.My)
    ws.shear_loads(Vx=ref.Vx, Vy=ref.Vy, T=ref.T)
    ws.compute_stresses()
    ws.sigma_plotter(show=False)

    if type_s == "Max":
        sensi_bending_t_sk[i] = ws.sigma_max
    if type_s == "Point":
        sensi_bending_t_sk[i] = ws.sigma1(s1)

###Changing the height of the stringer
for i,p in enumerate(percentages):
    h_str = (1+p)*ref.h_str_normal
    MOI = Inertia_normal_var(n_str=ref.n_str_normal,t_str=ref.t_str_normal,t_sp=ref.t_sp_normal,
                             t_sk=ref.t_sk_normal,h_str=h_str,w_str=ref.w_str_normal)
    MOI.compute_inertia(x=x)

    ws = stresses(Ixx=ref.Ixx_shear, Iyy=ref.Iyy_shear, Ixx_str=MOI.Ixx_normal, Iyy_str=MOI.Iyy_normal,
                  h=ref.h_sp_c_normal * ref.c, L=ref.w_sk_normal, t_upper=ref.t_sk_normal,
                  t_spar1=ref.t_sp_normal, t_spar2=ref.t_sp_normal, t_lower=ref.t_sk_normal)
    ws.bending_loads(Mx=ref.Mx, My=ref.My)
    ws.shear_loads(Vx=ref.Vx, Vy=ref.Vy, T=ref.T)
    ws.compute_stresses()
    ws.sigma_plotter(show=False)

    if type_s == "Max":
        sensi_bending_h_str[i] = ws.sigma_max
    if type_s == "Point":
        sensi_bending_h_str[i]= ws.sigma1(s1)

###SHEAR STRESS VARIATION

##Creating empty arrays for storing the results
sensi_shear_Vx = np.zeros(len(percentages))
sensi_shear_Vy = np.zeros(len(percentages))
sensi_shear_T = np.zeros(len(percentages))
sensi_shear_t_sp = np.zeros(len(percentages))
sensi_shear_t_sk = np.zeros(len(percentages))

###Changing the Vx load
for i,p in enumerate(percentages):
    if ref.Vx > 0:
        Vx = (1+p) * ref.Vx
    else:
        Vx = (1 - p) * ref.Vx
    ws = stresses(Ixx=ref.Ixx_shear, Iyy=ref.Iyy_shear, Ixx_str=ref.Ixx_normal, Iyy_str=ref.Iyy_normal,
                h=ref.h_sp_c_shear*ref.c, L=ref.w_sk_shear , t_upper=ref.t_sk_shear,
                t_spar1=ref.t_sp_shear, t_spar2=ref.t_sp_shear, t_lower=ref.t_sk_shear)

    ws.bending_loads(Mx = ref.Mx,My = ref.My)
    ws.shear_loads(Vx=Vx,Vy=ref.Vy,T=ref.T)
    ws.compute_stresses()
    ws.shear_flow_plotter(show=False,type="total")

    if type_s == "Max":
        sensi_shear_Vx[i] = ws.shear_stress_max
    if type_s == "Point":
        sensi_shear_Vx[i]= ws.q1_tot(s1)/ref.t_sk_shear

###Changing the Vy load
for i,p in enumerate(percentages):
    if ref.Vy > 0:
        Vy = (1+p) * ref.Vy
    else:
        Vy = (1-p) * ref.Vy

    ws = stresses(Ixx=ref.Ixx_shear, Iyy=ref.Iyy_shear, Ixx_str=ref.Ixx_normal, Iyy_str=ref.Iyy_normal,
                h=ref.h_sp_c_shear*ref.c, L=ref.w_sk_shear, t_upper=ref.t_sk_shear,
                t_spar1=ref.t_sp_shear, t_spar2=ref.t_sp_shear, t_lower=ref.t_sk_shear)

    ws.bending_loads(Mx = ref.Mx,My = ref.My)
    ws.shear_loads(Vx=ref.Vx,Vy=Vy,T=ref.T)
    ws.compute_stresses()
    ws.shear_flow_plotter(show=False,type="total")

    if type_s == "Max":
        sensi_shear_Vy[i] = ws.shear_stress_max
    if type_s == "Point":
        sensi_shear_Vy[i]= ws.q1_tot(s1)/ref.t_sk_shear

###Changing the T load
for i,p in enumerate(percentages):
    if ref.T>0:
        T = (1 + p) * ref.T
    else:
        T = (1 - p) * ref.T
    ws = stresses(Ixx=ref.Ixx_shear, Iyy=ref.Iyy_shear, Ixx_str=ref.Ixx_normal, Iyy_str=ref.Iyy_normal,
                h=ref.h_sp_c_shear*ref.c, L=ref.w_sk_shear, t_upper=ref.t_sk_shear,
                t_spar1=ref.t_sp_shear, t_spar2=ref.t_sp_shear, t_lower=ref.t_sk_shear)

    ws.bending_loads(Mx = ref.Mx,My = ref.My)
    ws.shear_loads(Vx=ref.Vx,Vy=ref.Vy,T=T)
    ws.compute_stresses()
    ws.shear_flow_plotter(show=False,type="total")

    if type_s == "Max":
        sensi_shear_T[i] = ws.shear_stress_max
    if type_s == "Point":
        sensi_shear_T[i] = ws.q1_tot(s1)/ref.t_sk_shear

###Changing the thickness of the spar
for i,p in enumerate(percentages):
    t_sp = (1+p) * ref.t_sp_shear

    MOI = Inertia_shear_var(t_sp=t_sp,t_sk=ref.t_sk_shear,n_str=ref.n_str_normal)
    MOI.compute_inertia(x=x)
    ws = stresses(Ixx=MOI.Ixx_shear, Iyy=MOI.Iyy_shear, Ixx_str=ref.Ixx_normal, Iyy_str=ref.Iyy_normal,
                h=ref.h_sp_c_shear*ref.c, L=ref.w_sk_shear, t_upper=ref.t_sk_shear,
                t_spar1=t_sp, t_spar2=t_sp, t_lower=ref.t_sk_shear)

    ws.bending_loads(Mx = ref.Mx,My = ref.My)
    ws.shear_loads(Vx=ref.Vx,Vy=ref.Vy,T=ref.T)
    ws.compute_stresses()
    ws.shear_flow_plotter(show=False,type="total")

    if type_s == "Max":
        sensi_shear_t_sp[i] = ws.shear_stress_max
    if type_s == "Point":
        sensi_shear_t_sp[i] = ws.q1_tot(s1)/ref.t_sk_shear

###Changing the thickness of the skin
for i, p in enumerate(percentages):
    t_sk = (1 + p) * ref.t_sk_shear

    MOI = Inertia_shear_var(t_sp=ref.t_sp_shear, t_sk=t_sk, n_str=ref.n_str_normal)
    MOI.compute_inertia(x=x)
    ws = stresses(Ixx=MOI.Ixx_shear, Iyy=MOI.Iyy_shear, Ixx_str=ref.Ixx_normal, Iyy_str=ref.Iyy_normal,
                  h=ref.h_sp_c_shear * ref.c, L=ref.w_sk_shear, t_upper=t_sk,
                  t_spar1=ref.t_sp_shear, t_spar2=ref.t_sp_shear, t_lower=t_sk)

    ws.bending_loads(Mx=ref.Mx, My=ref.My)
    ws.shear_loads(Vx=ref.Vx, Vy=ref.Vy, T=ref.T)
    ws.compute_stresses()
    ws.shear_flow_plotter(show=False, type="total")

    if type_s == "Max":
        sensi_shear_t_sk[i] = ws.shear_stress_max
    if type_s == "Point":
        sensi_shear_t_sk[i] = ws.q1_tot(s1)/t_sk

    print("Skin thickness",t_sk)
    print("Ixx shear",MOI.Ixx_shear)
    print("Iyy shear", MOI.Iyy_shear)

###Plotting
def bending_sensitivity_plot(Mx,My,t_str,t_sp,t_sk,h_str,imposed):
    if Mx == True:
        plt.plot(percentages, (sensi_bending_Mx - ref.max_sigma) / ref.max_sigma)
        plt.xlabel("Bending moment Mx percentage variation")
        plt.ylabel("Maximum Bending Stress Variation")
        plt.title("Change of maximum bending stress with variation of Mx ")
        plt.figure()

    if My == True:
        plt.plot(percentages, (sensi_bending_My - ref.max_sigma) / ref.max_sigma)
        plt.xlabel("Bending moment My percentage variation")
        plt.ylabel("Maximum Bending Stress Variation")
        plt.title("Change of maximum bending stress with variation of My ")
        plt.figure()

    if t_str == True:
        plt.plot(percentages, (sensi_bending_t_str - ref.max_sigma) / ref.max_sigma)
        plt.xlabel("Stringer thickness percentage variation")
        plt.ylabel("Maximum Bending Stress Variation")
        plt.title("Change of maximum bending stress with variation of t_str ")
        plt.figure()

    if t_sp == True:
        plt.plot(percentages, (sensi_bending_t_sp - ref.max_sigma) / ref.max_sigma)
        plt.xlabel("Spar thickness percentage variation")
        plt.ylabel("Maximum Bending Stress Variation")
        plt.title("Change of maximum bending stress with variation of t_sp ")
        plt.figure()

    if t_sk == True: #TODO CHECK THIS EFFECT
        plt.plot(percentages, (sensi_bending_t_sk - ref.max_sigma) / ref.max_sigma)
        plt.xlabel("Skin thickness percentage variation")
        plt.ylabel("Maximum Bending Stress Variation")
        plt.title("Change of maximum bending stress with variation of t_sk ")
        plt.figure()

    if h_str == True:
        plt.plot(percentages, (sensi_bending_h_str - ref.max_sigma) / ref.max_sigma)
        plt.xlabel("Stringer height percentage variation")
        plt.ylabel("Maximum Bending Stress Variation")
        plt.title("Change of maximum bending stress with variation of stringer height ")
        plt.figure()


    if imposed == True:
        plt.plot(percentages, (sensi_bending_Mx - ref.max_sigma) / ref.max_sigma,label="Mx")
        plt.plot(percentages, (sensi_bending_My - ref.max_sigma) / ref.max_sigma,label="My")
        plt.plot(percentages, (sensi_bending_t_str - ref.max_sigma) / ref.max_sigma,label="t_str")
        plt.plot(percentages, (sensi_bending_t_sp - ref.max_sigma) / ref.max_sigma,label="t_sp")
        plt.plot(percentages, (sensi_bending_t_sk - ref.max_sigma) / ref.max_sigma,label="t_sk")
        plt.plot(percentages, (sensi_bending_h_str - ref.max_sigma) / ref.max_sigma,label="h_str")
        plt.xlabel("Percentage variation of the parameter")
        plt.ylabel("Percentage variation of maximum bending stress")
        plt.legend()

    plt.show()

def shear_sensitivity_plot(Vx,Vy,T,t_sp,t_sk,imposed):

    if Vx==True:
        plt.plot(percentages, (sensi_shear_Vx - ref.max_shear) / ref.max_shear)
        plt.xlabel("Shear load Vx percentage variation")
        plt.ylabel("Maximum Shear Stress Variation")
        plt.title("Change of maximum shear stress with variation of Vx ")
        plt.figure()

    if Vy==True:
        plt.plot(percentages, (sensi_shear_Vy - ref.max_shear) / ref.max_shear)
        plt.xlabel("Shear load Vy percentage variation")
        plt.ylabel("Maximum Shear Stress Variation")
        plt.title("Change of maximum shear stress with variation of Vy ")
        plt.figure()

    if T==True:
        plt.plot(percentages, (sensi_shear_T- ref.max_shear) / ref.max_shear)
        plt.xlabel("Torque T percentage variation")
        plt.ylabel("Maximum Shear Stress Variation")
        plt.title("Change of maximum shear stress with variation of T")
        plt.figure()

    if t_sp==True:
        plt.plot(percentages, (sensi_shear_t_sp - ref.max_shear) / ref.max_shear)
        plt.xlabel("Spar thickness percentage variation")
        plt.ylabel("Maximum Shear Stress Variation")
        plt.title("Change of maximum shear stress with spar thickness variation")
        plt.figure()

    if t_sk == True:
        plt.plot(percentages, (sensi_shear_t_sk - ref.max_shear) / ref.max_shear)
        plt.xlabel("Skin thickness percentage variation")
        plt.ylabel("Maximum Shear Stress Variation")
        plt.title("Change of maximum shear stress with skin thickness variation")
        plt.figure()

    if imposed == True:
        plt.plot(percentages, (sensi_shear_Vx - ref.max_shear) / ref.max_shear,label="Vx")
        plt.plot(percentages, (sensi_shear_Vy - ref.max_shear) / ref.max_shear,label="Vy")
        plt.plot(percentages, (sensi_shear_T - ref.max_shear) / ref.max_shear,label="T")
        plt.plot(percentages, (sensi_shear_t_sp - ref.max_shear) / ref.max_shear,label="t_sp")
        plt.plot(percentages, (sensi_shear_t_sk - ref.max_shear) / ref.max_shear,label="t_sk")
        plt.xlabel("Percentage variation of the parameter")
        plt.ylabel("Percentage variation of maximum shear stress")
        plt.legend()

    plt.show()

#bending_sensitivity_plot(Mx=False,My=False,t_str=False,t_sp = False,t_sk=False,h_str=True,imposed=True)
#bending_sensitivity_plot(Mx=True,My=True,t_str=True,t_sp = True,t_sk=True,h_str=True,imposed=True)

shear_sensitivity_plot(Vx=True,Vy = True,T=True,t_sp=True,t_sk=True,imposed=True)





import numpy as np
import matplotlib.pyplot as plt
from Subsystem_design.Tank_Design import Mechanical_Design,Materials
from Subsystem_design.fuel_required import V_H2
import math as m

total_vol = V_H2 / 1000  # [m^3]

p_outside = 0.19028      # [bar] at 12100 [m] of altitude
p_tank = 3.5 * 1.5       # [bar] venting pressure
t_outside = 273.15 + 45  # [K]
t_tank = 19.75           # [K]
e_w = 0.8                # weld efficiency from Barron, 1985
s_a = 200                # [MPa] allowable stress dummy value for now!!
rho = 71.1               # [kg/m^3] liquid hydrogen density

dp = np.abs(p_tank-p_outside)
dt = np.abs(t_tank-t_outside)

class spacial_constraints():
    def __init__(self, length, width, height):
        """
        The sizes of the current A320 neo will represent dimensional constraints for fitting in the LH2 tanks.
        :param length:  [m]
        :param width:   [m]
        :param height:  [m]
        """

        self.length = length
        self.width = width
        self.height = height

class spacial_constraints_micha():
    def __init__(self, volume, width, height):
        """
        The sizes of the current A320 neo will represent dimensional constraints for fitting in the LH2 tanks.
        :param length:  [m]
        :param width:   [m]
        :param height:  [m]
        """

        self.volume = volume
        self.width = width
        self.height = height

def tank_length(vol,diameter1):
    """
    Outputs the length of the tank for a fixed diameter such that the required volume is attained.
    """
    length1 = 1
    pod_constaint1 = spacial_constraints(length=length1, width=diameter1, height=diameter1)
    pod_tank1 = Mechanical_Design.PodTank(constraints=pod_constaint1, dp=dp, s_a=s_a,
                                         e_w=e_w, material_insulation=Materials.MLI
                                         , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                         rho=rho, t_tank=t_tank,
                                         dt=dt, p_tank=(1 + p) * p_tank)
    pod_tank1.tank_design()
    dl = 0.01
    while pod_tank1.inner_vol_inner_wall < vol:
        length1 = length1 + dl
        pod_constaint1 = spacial_constraints(length=length1, width=diameter1, height=diameter1)
        pod_tank1 = Mechanical_Design.PodTank(constraints=pod_constaint1, dp=dp, s_a=s_a,
                                              e_w=e_w, material_insulation=Materials.MLI
                                              , material_inner=Materials.Al_2090_T81,
                                              material_outer=Materials.Al_2090_T81,
                                              rho=rho, t_tank=t_tank,
                                              dt=dt, p_tank=(1 + p) * p_tank)
        pod_tank1.tank_design()

    return length1

def tank_daiemet(vol,length1):
    """
        Outputs the diameter of the tank for a fixed length such that the required volume is attained.
       """
    diameter1 = 0.25
    pod_constaint1 = spacial_constraints(length=length1, width=diameter1, height=diameter1)
    pod_tank1 = Mechanical_Design.PodTank(constraints=pod_constaint1, dp=dp, s_a=s_a,
                                          e_w=e_w, material_insulation=Materials.MLI
                                          , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                          rho=rho, t_tank=t_tank,
                                          dt=dt, p_tank=(1 + p) * p_tank)
    pod_tank1.tank_design()
    dl = 0.01
    while pod_tank1.inner_vol_inner_wall < vol:
        diameter1 = diameter1 + dl
        pod_constaint1 = spacial_constraints(length=length1, width=diameter1, height=diameter1)
        pod_tank1 = Mechanical_Design.PodTank(constraints=pod_constaint1, dp=dp, s_a=s_a,
                                              e_w=e_w, material_insulation=Materials.MLI
                                              , material_inner=Materials.Al_2090_T81,
                                              material_outer=Materials.Al_2090_T81,
                                              rho=rho, t_tank=t_tank,
                                              dt=dt, p_tank=(1 + p) * p_tank)
        pod_tank1.tank_design()

    return diameter1


percentages = [-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5] #list of percentage increase
percentages = np.arange(-0.5,0.6,0.1)

###INSERT VALUES FOR CURRENTLY DESIGNED TANK
length_ref = 5.73 #[m]
diameter_ref = 2.23   #[m]

pod_reference = spacial_constraints(length=length_ref, width=diameter_ref, height=diameter_ref)
pod_tank_reference = Mechanical_Design.PodTank(constraints=pod_reference, dp=dp, s_a=s_a,
                                      e_w=e_w, material_insulation=Materials.MLI
                                      , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                      rho=Mechanical_Design.rho, t_tank=Mechanical_Design.t_tank,
                                      dt=Mechanical_Design.dt, p_tank=Mechanical_Design.p_tank)
pod_tank_reference.tank_design()

###CHANGING PRESSURE

sensi_pressure = np.zeros(len(percentages)) #list for storing the results of the sensitivity analysis caused by the pressure variation

for i,p in enumerate(percentages):

    pod_tank = Mechanical_Design.PodTank(constraints=pod_reference, dp=dp, s_a=s_a,
                                      e_w=e_w, material_insulation=Materials.MLI
                                      , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                      rho=rho, t_tank=t_tank,
                                      dt=dt, p_tank = (1 + p) * p_tank)
    pod_tank.tank_design()
    if p ==0:
        sensi_pressure[i] = pod_tank_reference.mass_tank
    else:
        sensi_pressure[i] = pod_tank.mass_tank

###CHANGING DIAMETER

sensi_diameter = np.zeros(len(percentages)) #list for storing the results of the sensitivity analysis caused by the diameter variation
sensi_diameter_micha = np.zeros(len(percentages))
for i, p in enumerate(percentages):
    #For each iteration we need to create new constraints
    diameter = (1+p)*diameter_ref
    #length = (total_vol/2 - 4/3 * m.pi * (diameter/2)**2)/(m.pi * pow(diameter/2,2))
    length = tank_length(vol = total_vol/2,diameter1=diameter)
    pod_constaint = spacial_constraints(length=length,width=diameter,height=diameter)
    pod_tank = Mechanical_Design.PodTank(constraints=pod_constaint, dp=dp, s_a=s_a,
                                         e_w=e_w, material_insulation=Materials.MLI
                                         , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                         rho=rho, t_tank=t_tank,
                                         dt=dt, p_tank=p_tank)
    pod_tank.tank_design()

    pod_constraint_micha = spacial_constraints_micha(volume=total_vol/2,width=diameter,height=diameter)
    pod_tank_micha = Mechanical_Design.OnlyPods(constraints=pod_constraint_micha,dp=dp, s_a=s_a, e_w=e_w,
                                                material_insulation=Materials.MLI,material_inner=Materials.Al_2090_T81,
                                                material_outer=Materials.Al_2090_T81,rho=rho,t_tank=t_tank,dt=dt,p_tank=p_tank)
    pod_tank_micha.tank_design()

    if p == 0:
        sensi_diameter[i] = pod_tank_reference.mass_tank
        sensi_diameter_micha[i]  = pod_tank_reference.mass_tank
    else:
        sensi_diameter[i] = pod_tank.mass_tank
        sensi_diameter_micha[i] = pod_tank.mass_tank


###CHANGING LENGTH

sensi_length = np.zeros(len(percentages)) #list for storing the results of the sensitivity analysis caused by the diameter variation
for i,p in enumerate(percentages):
    #For each iteration we need to create new constraints
    length = (1+p) * length_ref
    #diameter = m.sqrt(total_vol/(m.pi*length + m.pi * 4/3))
    diameter = tank_daiemet(vol=total_vol/2,length1=length)
    pod_constaint = spacial_constraints(length=length, width=diameter, height=diameter)
    pod_tank = Mechanical_Design.PodTank(constraints=pod_constaint, dp=dp, s_a=s_a,
                                         e_w=e_w, material_insulation=Materials.MLI
                                         , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                         rho=rho, t_tank=t_tank,
                                         dt=dt, p_tank=p_tank)
    pod_tank.tank_design()

    if p ==0:
        sensi_length[i] = pod_tank_reference.mass_tank
    else:
        sensi_length[i] = pod_tank.mass_tank


# plt.plot(percentages,(sensi_pressure-pod_tank_reference.mass_tank)/pod_tank_reference.mass_tank,label="Venting pressure")
# plt.title("Variation of tank mass with variation of venting pressure")
# plt.xlabel("Percentage change in venting pressure")
# plt.ylabel("Percentage change in mass")
# plt.figure()

plt.plot(percentages,(sensi_diameter-pod_tank_reference.mass_tank)/pod_tank_reference.mass_tank,label = "Diameter")
plt.title("Variation of tank mass with variation of diameter")
plt.xlabel("Percentage change in diameter")
plt.ylabel("Percentage change in mass")
# plt.figure()

# plt.plot(percentages,(sensi_diameter_micha-pod_tank_reference.mass_tank)/pod_tank_reference.mass_tank,label = "Diameter")
# plt.title("Variation of tank mass with variation of diameter(Micha program)")
# plt.xlabel("Percentage change in diameter")
# plt.ylabel("Percentage change in mass")
# plt.figure()

# plt.plot((1+np.array(percentages))*diameter_ref,sensi_diameter_micha,label = "Diameter")
# plt.title("Variation of tank mass with variation of diameter(Micha program)")
# plt.xlabel("Diameter[m]")
# plt.ylabel("Mass[kg]")

plt.plot(percentages,(sensi_length-pod_tank_reference.mass_tank)/pod_tank_reference.mass_tank,label = "Length")
plt.title("Variation of tank mass with variation of tank length")
plt.xlabel("Percentage change in length")
plt.ylabel("Percentage change in mass")
plt.xlim(xmin=-0.6,xmax=0.6)
plt.xticks(ticks=percentages)
plt.legend()
plt.grid()
plt.show()

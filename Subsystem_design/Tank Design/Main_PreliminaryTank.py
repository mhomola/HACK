import numpy as np
import matplotlib.pyplot as plt

import Mechanical_Design
import Materials
import math as m

central_H2_vol = 15 # [m^3] the volume of LH2 in the central tanks area
aft_H2_vol = 30     # [m^3] the volume of LH2 that must be stored behind the passenger cabin

class spacial_constraints():
    def __init__(self,length,width,height):
        """
        The sizes of the current A320 neo will represent dimensional constraints for fitting in the LH2 tanks.
        :param length:  [m]
        :param width:   [m]
        :param height:  [m]
        """

        self.length = length
        self.width = width
        self.height = height

#Define the areas of the constraints

central = spacial_constraints(length=3.5,width=1.14,height=2.67)

central_iter = spacial_constraints(length=6.5,width=1.14,height=2.67)
central_iter.width = central_iter.width *0.7 #we start from a smaller diameter

if __name__ == '__main__':
    iter_cond = 1 #will be sitched to 0 when the outerdiameter of the tank will surpass the central.width

    while central_iter.width < central.width and iter_cond==1:
        central_tank = Mechanical_Design.Tank(constraints=central_iter,dp=Mechanical_Design.dp, s_a=Mechanical_Design.s_a, e_w=Mechanical_Design.e_w,material_insulation=Materials.MLI
                          ,material_inner = Materials.Al_2090_T81,material_outer=Materials.Al_2090_T81,rho=Mechanical_Design.rho,t_tank=Mechanical_Design.t_tank,dt=Mechanical_Design.dt,p_tank=Mechanical_Design.p_tank)
        central_tank.tank_design()
        if (central.width*0.99 <= central_tank.r4*2): #1% margin
            iter_cond = 0
        print(central_iter.width)
        central_iter.width = central_iter.width + 0.01

    print("#######CENTRAL TANKS######")
    print("Mass of one tank",central_tank.mass_tank,"[kg]")
    print("Insulation of the tank thickness",central_tank.t_insulation*100,"[cm]")
    print("Thickness inner wall",central_tank.t_wall_inner*1000,"[mm]")
    print("Thickness outer wall",central_tank.t_wall_outer *1000,"[mm]")
    print("H2 volume that is stored per tank",central_tank.inner_vol_inner_wall,"[m^3]")
    print("Outer Diameter of the tank[including insulation]",central_tank.r4*2,"[m]")
    print("Number of tanks tht can be accommodated in the central sector",m.floor(central.height/(central_tank.r4*2)))

    central_H2_vol = central_tank.inner_vol_inner_wall * m.floor(central.height/(central_tank.r4*2)) #[m^3]
    #aft_H2_vol = -central_H2_vol

#ACT1 = spacial_constraint(length=,width=,height=)
#ACT1 = spacial_constraint(length=,width=,height=)
#aft = spacial_constraint(length=,width=,height=)
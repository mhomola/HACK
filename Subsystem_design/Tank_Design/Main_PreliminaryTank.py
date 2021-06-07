import Mechanical_Design
import Materials
import numpy as np
import math as m
from Subsystem_design.fuel_required import V_H2

total_vol = V_H2 / 1000 #[m^3]

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

#Define the areas of the constraints


central = spacial_constraints(length=3.14, width=1.24, height=2.67)

central_iter = spacial_constraints(length=3.14, width=1.24, height=2.67)
central_iter.width = central_iter.width * 0.7 #we start from a smaller diameter

iter_cond = 1 #will be sitched to 0 when the outerdiameter of the tank will surpass the central.width

#Iterations for the central tanks
while central_iter.width < central.width and iter_cond==1:
    central_tank = Mechanical_Design.Tank(constraints=central_iter,dp=Mechanical_Design.dp, s_a=Mechanical_Design.s_a, e_w=Mechanical_Design.e_w,material_insulation=Materials.MLI
                      ,material_inner = Materials.Al_2090_T81,material_outer=Materials.Al_2090_T81,rho=Mechanical_Design.rho,t_tank=Mechanical_Design.t_tank,dt=Mechanical_Design.dt,p_tank=Mechanical_Design.p_tank)
    central_tank.tank_design()
    if (central.width*0.99 <= central_tank.r4*2): #1% margin
        iter_cond = 0
    central_iter.width = central_iter.width + 0.01

central_H2_vol = central_tank.inner_vol_inner_wall * m.floor(central.height/(central_tank.r4*2)) #[m^3]the volume of LH2 in the central tanks area
aft_H2_vol = total_vol - central_H2_vol#[m^3] the volume of LH2 that must be stored behind the passenger cabin

# Iteration for the tank behind the cabin
aft = spacial_constraints(length=aft_H2_vol / (1.78 ** 2 * m.pi), width=3.56, height=3.56)
aft_tank = Mechanical_Design.Cyl_Tank(constraints=aft, dp=Mechanical_Design.dp, s_a=Mechanical_Design.s_a,
                                      e_w=Mechanical_Design.e_w, material_insulation=Materials.MLI
                                      , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                      rho=Mechanical_Design.rho, t_tank=Mechanical_Design.t_tank,
                                      dt=Mechanical_Design.dt, p_tank=Mechanical_Design.p_tank)
aft_tank.tank_design()

weight_addition = aft_tank.mass_tank + 2 * central_tank.mass_tank

ijk = 2 * central_tank.mass_tank


if __name__ == '__main__':

    print("#######CENTRAL TANKS######")
    print("Mass of one tank", central_tank.mass_tank, "[kg]")
    print("Insulation of the tank thickness", central_tank.t_insulation*100, "[cm]")
    print("Thickness inner wall", central_tank.t_wall_inner*1000, "[mm]")
    print("Thickness outer wall", central_tank.t_wall_outer *1000, "[mm]")
    print("H2 volume that is stored per tank", central_tank.inner_vol_inner_wall, "[m^3]")
    print("Outer Diameter of the tank[including insulation]", central_tank.r4*2, "[m]")
    print("Number of tanks tht can be accommodated in the central sector", m.floor(central.height/(central_tank.r4*2)))

    #Iteration for the tank behind the cabin

    print("#######AFT TANKS######")
    print("Mass of one tank", aft_tank.mass_tank, "[kg]")
    print("Insulation of the tank thickness", aft_tank.t_insulation * 100, "[cm]")
    print("Thickness inner wall", aft_tank.t_wall_inner * 1000, "[mm]")
    print("Thickness outer wall", aft_tank.t_wall_outer * 1000, "[mm]")
    print("H2 volume that is stored per tank", aft_tank.inner_vol_inner_wall, "[m^3]")
    print("Outer Diameter of the tank[including insulation]", aft_tank.r4 * 2, "[m]")
    print("Length", aft.length, "#m")

    print("Added Weight = ", weight_addition)

    Q_req = [central_tank.Q_req,aft_tank.Q_req]
    Q_req = np.array(Q_req)
    np.save("Q_required.npy",Q_req)



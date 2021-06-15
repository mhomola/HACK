import Subsystem_design.Tank_Design.Mechanical_Design as Mechanical_Design
import Subsystem_design.Tank_Design.Materials as Materials
import numpy as np
import math as m
from Subsystem_design.fuel_required import V_H2


class spacial_constraints():
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


""" ------ Pod Tank Design ------- """
total_vol = V_H2 / 1000  # [m^3]
pod_H2_vol = total_vol/2  # [m^3] the volume of LH2 in each pod

pod = spacial_constraints(volume=pod_H2_vol, width=2.2, height=2.2)
pod_tank = Mechanical_Design.OnlyPods(constraints=pod, dp=Mechanical_Design.dp, s_a=Mechanical_Design.s_a,
                                      e_w=Mechanical_Design.e_w, material_insulation=Materials.MLI
                                      , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                      rho=Mechanical_Design.rho, t_tank=Mechanical_Design.t_tank,
                                      dt=Mechanical_Design.dt, p_tank=Mechanical_Design.p_tank)
pod_tank.tank_design()

weight_addition = pod_tank.mass_tank * 2
volume_all = pod_tank.inner_vol_inner_wall * 2

""" ---- Important Parameters of Tanks ----"""
l_wing_pod = pod_tank.length
d_wing_pod = pod_tank.r4 * 2
mass_pod = pod_tank.mass_tank
volume_pod = pod_tank.inner_vol_inner_wall

if __name__ == '__main__':

    print("#######PODDED TANKS#######")
    print("Mass of one tank", mass_pod, "[kg]")
    print("Insulation of the tank thickness", pod_tank.t_insulation * 100, "[cm]")
    print("Thickness inner wall", pod_tank.t_wall_inner * 1000, "[mm]")
    print("Thickness outer wall", pod_tank.t_wall_outer * 1000, "[mm]")
    print("H2 volume that is stored per tank", volume_pod, "[m^3]")
    print("Outer Diameter of the tank[including insulation]", d_wing_pod, "[m]")
    print("Length", l_wing_pod, "[m]")
    print("Required Heat = ", pod_tank.Q_req, "vs Designed Heat =  ", pod_tank.Q_tot)

    print("#######OVERALL#######")
    print("Added Weight = ", weight_addition, "[kg]")
    print("Volume of all Tanks = ", volume_all, "[m^3]", "complying with the required ", V_H2/1000, "[m^3]")



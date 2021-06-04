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


""" ------ Central Tank Design ------- """
central = spacial_constraints(length=3.14, width=1.24, height=2.67)
central_iter = spacial_constraints(length=3.14, width=1.24, height=2.67)
central_iter.width = central_iter.width * 0.8                           # start from a smaller diameter

iter_cond = 1                 # will be switched to 0 when the outerdiameter of the tank will surpass the central.width

# Iterations for the central tanks
while central_iter.width < central.width and iter_cond == 1:
    central_tank = Mechanical_Design.Tank(constraints=central_iter, dp=Mechanical_Design.dp, s_a=Mechanical_Design.s_a, e_w=Mechanical_Design.e_w, material_insulation=Materials.MLI,
                                          material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81, rho=Mechanical_Design.rho, t_tank=Mechanical_Design.t_tank, dt=Mechanical_Design.dt, p_tank=Mechanical_Design.p_tank)

    central_tank.tank_design()
    if central.width*0.99 <= central_tank.r4*2:            # 1% margin
        iter_cond = 0
    central_iter.width = central_iter.width + 0.01

# [m^3]the volume of LH2 in the central tanks area
central_H2_vol = central_tank.inner_vol_inner_wall * m.floor(central.height/(central_tank.r4*2))

""" ------ Pod Tank Design ------- """

pod_H2_vol = (total_vol - central_H2_vol)/2  # [m^3] the volume of LH2 in each pod

pod = spacial_constraints(length=5.765, width=2, height=2)
pod_tank = Mechanical_Design.PodTank(constraints=pod, dp=Mechanical_Design.dp, s_a=Mechanical_Design.s_a,
                                      e_w=Mechanical_Design.e_w, material_insulation=Materials.MLI
                                      , material_inner=Materials.Al_2090_T81, material_outer=Materials.Al_2090_T81,
                                      rho=Mechanical_Design.rho, t_tank=Mechanical_Design.t_tank,
                                      dt=Mechanical_Design.dt, p_tank=Mechanical_Design.p_tank)
pod_tank.tank_design()

weight_addition = central_tank.mass_tank * 2 + pod_tank.mass_tank * 2
volume_all = central_tank.inner_vol_inner_wall * 2 + pod_tank.inner_vol_inner_wall * 2

if __name__ == '__main__':

    print("#######CENTRAL TANKS######")
    print("Mass of one tank", central_tank.mass_tank, "[kg]")
    print("Insulation of the tank thickness", central_tank.t_insulation*100, "[cm]")
    print("Thickness inner wall", central_tank.t_wall_inner*1000, "[mm]")
    print("Thickness outer wall", central_tank.t_wall_outer *1000, "[mm]")
    print("H2 volume that is stored per tank", central_tank.inner_vol_inner_wall, "[m^3]")
    print("Outer Diameter of the tank[including insulation]", central_tank.r4*2, "[m]")
    print("Number of tanks tht can be accommodated in the central sector", m.floor(central.height/(central_tank.r4*2)))
    print("Required Heat = ", central_tank.Q_req, "vs Designed Heat =  ", central_tank.Q_tot)

    print("#######PODDED TANKS######")
    print("Mass of one tank", pod_tank.mass_tank, "[kg]")
    print("Insulation of the tank thickness", pod_tank.t_insulation * 100, "[cm]")
    print("Thickness inner wall", pod_tank.t_wall_inner * 1000, "[mm]")
    print("Thickness outer wall", pod_tank.t_wall_outer * 1000, "[mm]")
    print("H2 volume that is stored per tank", pod_tank.inner_vol_inner_wall, "[m^3]")
    print("Outer Diameter of the tank[including insulation]", pod_tank.r4 * 2, "[m]")
    print("Length", pod.length, "[m]")
    print("Required Heat = ", pod_tank.Q_req, "vs Designed Heat =  ", pod_tank.Q_tot)

    print("Added Weight = ", weight_addition, "[kg]")
    print("Volume of all Tanks = ", volume_all, "[m^3]", "complying with the required ", V_H2/1000, "[m^3]")


    Q_req = [central_tank.Q_req, pod_tank.Q_req]
    Q_req = np.array(Q_req)
    np.save("Q_required.npy", Q_req)



"""
This files contains the classes and functions for determining the mechanical design of the tank (thickness of the wall).
"""
import numpy as np
import matplotlib.pyplot as plt
import math as m


# Hydrogen tank storage data
p_outside = 0.19028      # [bar] at 12100 [m] of altitude
p_tank = 3.5 * 1.5       # [bar] venting pressure
t_outside = 273.15 + 45  # [K]
t_tank = 19.75           # [K]
e_w = 0.8                # weld efficiency from Barron, 1985
s_a = 200                # [MPa] allowable stress dummy value for now!!
rho = 71.1                # [kg/m^3] liquid hydrogen density

# Determining the heat transfer rate that can be accepted based on the user requirement of no boil-off for 36h.
# Method from "Passive zero-boil-off storage of liquid hydrogen for long-time space missions"
dp = np.abs(p_tank-p_outside)
dt = np.abs(t_tank-t_outside)

# choose tank dimensions depending on constraints
# import the constraints for the Main file
# we maximize the volume of central tanks
# we assume circular cylinders with two hemisphere heads at ends


class Tank():
  def __init__(self,constraints,dp, s_a, e_w,material_insulation,material_inner,material_outer,rho,t_tank,dt,p_tank):
    """
    :param constraints: object of class constraints
    :param dp: pressure difference [bar]
    :param s_a: [Pa] allowable stress
    :param e_w: weld efficiency
    :param material_insulation: [material class object]
    :param material_inner: [material class object]
    :param material_outer: [material class object]
    :param rho: liquid hydrogen mass density [kg/m^3]
    :param t_tank: temperature inside the tank [K]
    :param dt:difference between temp inside th tank and outside temp
    """
    self.d_0 = min(constraints.width,constraints.height) #[m]
    self.length = constraints.length #including the hemisphere ends
    self.dp = dp * 10**5 #[Pa]
    self.s_a = s_a
    self.e_w = e_w
    self.K = 1/6 * (2 + 1) #circular condition
    self.material_insulation = material_insulation
    self.material_inner = material_inner
    self.material_outer = material_outer
    self.rho = rho
    self.t_tank = t_tank
    self.t_hot = 273.15 + 45        # [K]
    self.dt = dt
    self.p_tank = p_tank * 10**5    # [Pa]
    self.safety_factor = 1.5
    self.outer_vol_inner_wall = m.pi * (self.d_0/2)**2 * (self.length-self.d_0) + 4/3* m.pi * (self.d_0/2)**3 #cylinder volume + sphere volume

  def inner_wall(self):       # D_Verstraete_Thesis_2009
    t_wall_inner = self.p_tank * self.d_0 / (2 * self.material_inner.yield_strength *10**6 * self.e_w + 0.8 * self.p_tank) #10^6 factor for MPa ->Pa
    t_caps_inner = self.p_tank * self.d_0 * self.K / (2 * self.material_inner.yield_strength*10**6 * self.e_w + 2 * self.p_tank * (self.K - 0.1)) #we use directly p_tank cause p_tank - p_vac = p_tank

    self.t_wall_inner = t_wall_inner * self.safety_factor
    self.t_caps_inner = t_caps_inner * self.safety_factor

    self.inner_vol_inner_wall = m.pi * (self.d_0/2-self.t_wall_inner)**2 * (self.length-self.d_0) + 4/3* m.pi * (self.d_0/2-self.t_caps_inner)**3 #volume occupied by LH2
    self.mass_H2 = self.inner_vol_inner_wall * self.rho #calculates the mass of hydrogen that is encapsulated in the tank
    self.vol_wall_inner_wall = self.outer_vol_inner_wall - self.inner_vol_inner_wall #this will determine the mass of the inner wall of the tank
    self.mass_inner_wall = self.vol_wall_inner_wall*self.material_inner.density

  def insulation(self):
    self.r1 = self.d_0 / 2 - self.t_wall_inner
    self.r2 = self.d_0 / 2

    # Passive zero-boil-off storage of liquid hydrogen for long-time space missions
    t_boil = 20  # [K]
    Cp = 10310                            # [J/Kg*K] H2 specific heat capacity at 20[K]
    time = 36                             # [h] no boil-off should occur within 36h
    Q_req = Cp * self.mass_H2 * (t_boil - self.t_tank) / (time*3600)
    self.Q_req = Q_req

    Cs = 7.3 * 10**-8
    Cr = 7.07 * 10**-10
    Cg = 1.46 * 10**-4
    epsl = 0.021  # Effective Emissivity

    P = 1.3 * 10**-3 * 0.0075006168  # 1.3 mPa -> Torr
    N_density = 22.7  # [layers/cm]

    N_layers = 60  # number of layers
    q_tot = (Cr * epsl * (self.t_hot ** 4.67 - self.t_tank ** 4.67) / N_layers) + \
            (Cs * N_density ** 2.63 * (self.t_hot - self.t_tank) * (self.t_hot + self.t_tank) / (2 * N_layers + 2)) + \
            (Cg * P * (self.t_hot ** 0.52 - self.t_tank ** 0.52) / N_layers)
    self.Q_tot = (2 * np.pi * self.d_0/2 * (self.length - self.d_0) + 4 * np.pi * (self.d_0/2) ** 2) * q_tot

    self.t_insulation = N_layers / N_density * 0.01
    self.r3 = self.r2 + self.t_insulation
    self.outer_vol_insulation = m.pi * self.r3**2 * (self.length-self.d_0) + 4/3 * m.pi * (self.r3)**3
    self.vol_insulation = self.outer_vol_insulation - self.outer_vol_inner_wall # the volume of the insulation that will be used for computing the mass
    self.mass_insulation = self.vol_insulation * self.material_insulation.density

  def outer_wall(self):
    self.t_wall_outer = self.t_wall_inner      # we assume them equal for the moment
    self.r4 = self.r3 + self.t_wall_outer
    self.outer_vol_outer_wall = m.pi * self.r4**2 * (self.length-self.d_0) + 4/3* m.pi * self.r4**3
    self.vol_outer_wall = self.outer_vol_outer_wall - self.outer_vol_insulation
    self.mass_outer_wall = self.vol_outer_wall * self.material_outer.density

  def tank_design(self):
    self.inner_wall()
    self.insulation()
    self.outer_wall()
    self.mass_tank = self.mass_inner_wall + self.mass_insulation + self.mass_outer_wall


class PodTank():
  def __init__(self,constraints,dp, s_a, e_w,material_insulation,material_inner,material_outer,rho,t_tank,dt,p_tank):
    """
    :param constraints: object of class constraints
    :param dp: pressure difference [bar]
    :param s_a: [Pa] allowable stress
    :param e_w: weld efficiency
    :param material_insulation: [material class object]
    :param material_inner: [material class object]
    :param material_outer: [material class object]
    :param rho: liquid hydrogen mass density [kg/m^3]
    :param t_tank: temperature inside the tank [K]
    :param dt:difference between temp inside th tank and outside temp
    """
    self.d_0 = min(constraints.width, constraints.height) #[m]
    self.length = constraints.length # including the hemisphere ends
    self.dp = dp * 10**5 # [Pa]
    self.s_a = s_a
    self.e_w = e_w
    self.K = 1/6 * (2 + 1) # circular condition
    self.material_insulation = material_insulation
    self.material_inner = material_inner
    self.material_outer = material_outer
    self.rho = rho
    self.t_tank = t_tank
    self.t_hot = 273.15 + 45  # [K]
    self.dt = dt
    self.p_tank = p_tank * 10**5 # [Pa]
    self.safety_factor = 1.5
    self.outer_vol_inner_wall = m.pi * (self.d_0/2)**2 * (self.length-self.d_0) + 4/3* m.pi * (self.d_0/2)**3 #cylinder volume + sphere volume

  def inner_wall(self):       # D_Verstraete_Thesis_2009
    t_wall_inner = self.p_tank * self.d_0 / (2 * self.material_inner.yield_strength *10**6 * self.e_w + 0.8 * self.p_tank) #10^6 factor for MPa ->Pa
    t_caps_inner = self.p_tank * self.d_0 * self.K / (2 * self.material_inner.yield_strength*10**6 * self.e_w + 2 * self.p_tank * (self.K - 0.1)) #we use directly p_tank cause p_tank - p_vac = p_tank

    self.t_wall_inner = t_wall_inner * self.safety_factor
    self.t_caps_inner = t_caps_inner * self.safety_factor

    self.inner_vol_inner_wall = m.pi * (self.d_0/2-self.t_wall_inner)**2 * (self.length-self.d_0) + 4/3* m.pi * (self.d_0/2-self.t_caps_inner)**3 #volume occupied by LH2
    self.mass_H2 = self.inner_vol_inner_wall * self.rho #calculates the mass of hydrogen that is encapsulated in the tank
    self.vol_wall_inner_wall = self.outer_vol_inner_wall - self.inner_vol_inner_wall #this will determine the mass of the inner wall of the tank
    self.mass_inner_wall = self.vol_wall_inner_wall*self.material_inner.density

  def insulation(self):
    self.r1 = self.d_0 / 2 - self.t_wall_inner
    self.r2 = self.d_0 / 2

    #Passive zero-boil-off storage of liquid hydrogen for long-time space missions
    t_boil = 20  # [K]
    Cp = 10310  # [J/Kg*K] H2 specific heat capacity at 20[K]
    time = 36  # [h] no boil-off should occur within 36h
    Q_req = Cp * self.mass_H2 * (t_boil - self.t_tank) / (time*3600)
    self.Q_req = Q_req

    Cs = 7.3 * 10 ** -8
    Cr = 7.07 * 10 ** -10
    Cg = 1.46 * 10 ** -4
    epsl = 0.021  # Effective Emissivity

    P = 1.3 * 10 ** -3 * 0.0075006168  # 1.3 mPa -> Torr
    N_density = 22.7  # [layers/cm]

    N_layers = 35  # number of layers
    q_tot = (Cr * epsl * (self.t_hot ** 4.67 - self.t_tank ** 4.67) / N_layers) + \
            (Cs * N_density ** 2.63 * (self.t_hot - self.t_tank) * (self.t_hot + self.t_tank) / (2 * N_layers + 2)) + \
            (Cg * P * (self.t_hot ** 0.52 - self.t_tank ** 0.52) / N_layers)
    self.Q_tot = (2 * np.pi * self.d_0 / 2 * (self.length - self.d_0) + 4 * np.pi * (self.d_0 / 2) ** 2) * q_tot

    self.t_insulation = N_layers / N_density * 0.01
    self.r3 = self.r2 + self.t_insulation
    self.outer_vol_insulation = m.pi * self.r3**2 * (self.length-self.d_0) + 4/3 * m.pi * (self.r3)**3
    self.vol_insulation = self.outer_vol_insulation - self.outer_vol_inner_wall # the volume of the insulation that will be used for computing the mass
    self.mass_insulation = self.vol_insulation * self.material_insulation.density

  def outer_wall(self):
    self.t_wall_outer = self.t_wall_inner      # we assume them equal for the moment
    self.r4 = self.r3 + self.t_wall_outer
    self.outer_vol_outer_wall = m.pi * self.r4**2 * (self.length-self.d_0) + 4/3* m.pi * self.r4**3
    self.vol_outer_wall = self.outer_vol_outer_wall - self.outer_vol_insulation
    self.mass_outer_wall = self.vol_outer_wall * self.material_outer.density

  def tank_design(self):
    self.inner_wall()
    self.insulation()
    self.outer_wall()
    self.mass_tank = self.mass_inner_wall + self.mass_insulation + self.mass_outer_wall


class OnlyPods():
  def __init__(self,constraints,dp, s_a, e_w,material_insulation,material_inner,material_outer,rho,t_tank,dt,p_tank):
    """
    :param constraints: object of class constraints
    :param dp: pressure difference [bar]
    :param s_a: [Pa] allowable stress
    :param e_w: weld efficiency
    :param material_insulation: [material class object]
    :param material_inner: [material class object]
    :param material_outer: [material class object]
    :param rho: liquid hydrogen mass density [kg/m^3]
    :param t_tank: temperature inside the tank [K]
    :param dt:difference between temp inside th tank and outside temp
    """
    self.d_0 = min(constraints.width, constraints.height) # [m]
    self.volume = constraints.volume   # including the hemisphere ends
    self.dp = dp * 10**5               # [Pa]
    self.s_a = s_a
    self.e_w = e_w
    self.K = 1/6 * (2 + 1) # circular condition
    self.material_insulation = material_insulation
    self.material_inner = material_inner
    self.material_outer = material_outer
    self.rho = rho
    self.t_tank = t_tank
    self.t_hot = 273.15 + 45        # [K]
    self.dt = dt
    self.p_tank = p_tank * 10**5    # [Pa]
    self.safety_factor = 1.5

  def inner_wall(self):       # D_Verstraete_Thesis_2009
    t_wall_inner = self.p_tank * self.d_0 / (2 * self.material_inner.yield_strength *10**6 * self.e_w + 0.8 * self.p_tank) #10^6 factor for MPa ->Pa
    t_caps_inner = self.p_tank * self.d_0 * self.K / (2 * self.material_inner.yield_strength*10**6 * self.e_w + 2 * self.p_tank * (self.K - 0.1)) #we use directly p_tank cause p_tank - p_vac = p_tank

    self.t_wall_inner = t_wall_inner * self.safety_factor
    self.t_caps_inner = t_caps_inner * self.safety_factor

    self.inner_vol_inner_wall = self.volume
    self.length = (self.inner_vol_inner_wall - 4/3 * m.pi * (self.d_0/2-self.t_caps_inner)**3) /\
                  (m.pi * (self.d_0/2-self.t_wall_inner)**2) + self.d_0

    self.outer_vol_inner_wall = m.pi * (self.d_0/2)**2 * (self.length-self.d_0) + 4/3 * m.pi * (self.d_0/2)**3
    self.mass_H2 = self.inner_vol_inner_wall * self.rho #calculates the mass of hydrogen that is encapsulated in the tank
    self.vol_wall_inner_wall = self.outer_vol_inner_wall - self.inner_vol_inner_wall #this will determine the mass of the inner wall of the tank
    self.mass_inner_wall = self.vol_wall_inner_wall*self.material_inner.density

  def insulation(self):
    self.r1 = self.d_0 / 2 - self.t_wall_inner
    self.r2 = self.d_0 / 2

    # Passive zero-boil-off storage of liquid hydrogen for long-time space missions
    t_boil = 20  # [K]
    Cp = 10310                            # [J/Kg*K] H2 specific heat capacity at 20[K]
    time = 36                             # [h] no boil-off should occur within 36h
    self.Q_req = Cp * self.mass_H2 * (t_boil - self.t_tank) / (time*3600)


    Cs = 7.3 * 10**-8
    Cr = 7.07 * 10**-10
    Cg = 1.46 * 10**-4
    epsl = 6.8 * 10**-4 * (self.t_hot ** 0.67 - self.t_tank ** 0.67)    # Effective Emissivity

    P = 1.3 * 10**-3 * 0.0075006168  # 1.3 mPa -> Torr
    N_density = 22.7  # [layers/cm]

    N_layers = 34 # number of layers
    q_tot = (Cr * epsl * (self.t_hot ** 4.67 - self.t_tank ** 4.67) / N_layers) + \
            (Cs * N_density ** 2.63 * (self.t_hot - self.t_tank) * (self.t_hot + self.t_tank) / (2 * N_layers + 2)) + \
            (Cg * P * (self.t_hot ** 0.52 - self.t_tank ** 0.52) / N_layers)
    self.Q_tot = (2 * np.pi * self.d_0/2 * (self.length - self.d_0) + 4 * np.pi * (self.d_0/2) ** 2) * q_tot

    self.t_insulation = N_layers / N_density * 0.01
    self.r3 = self.r2 + self.t_insulation
    self.outer_vol_insulation = m.pi * self.r3**2 * (self.length-self.d_0) + 4/3 * m.pi * (self.r3)**3
    self.vol_insulation = self.outer_vol_insulation - self.outer_vol_inner_wall # the volume of the insulation that will be used for computing the mass
    self.mass_insulation = self.vol_insulation * self.material_insulation.density

  def outer_wall(self):
    self.t_wall_outer = self.t_wall_inner      # we assume them equal for the moment
    self.r4 = self.r3 + self.t_wall_outer
    self.outer_vol_outer_wall = m.pi * self.r4**2 * (self.length-self.d_0) + 4/3* m.pi * self.r4**3
    self.vol_outer_wall = self.outer_vol_outer_wall - self.outer_vol_insulation
    self.mass_outer_wall = self.vol_outer_wall * self.material_outer.density

  def tank_design(self):
    self.inner_wall()
    self.insulation()
    self.outer_wall()
    self.mass_tank = self.mass_inner_wall + self.mass_insulation + self.mass_outer_wall


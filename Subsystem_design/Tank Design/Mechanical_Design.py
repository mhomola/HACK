"""
This files contains the classes and functions for determining the mechanical design of the tank (thickness of the wall).
"""
import numpy as np
import matplotlib.pyplot as plt
import math as m

#Hydrogen tank storage data
p_outside = 0.19028  # [bar] at 12100 [m] of altitude
p_tank = 3           # [bar]
t_outside = 216.65   # [k]
t_tank = 20          # [k]
d_0 = 4              # [m] outside diameter -- dummy
e_w = 0.8            # weld efficiency from Barron, 1985
s_a = 200            # [MPa] allowable stress dummy value for now!!

dp = np.abs(p_tank-p_outside)
dt = np.abs(t_tank-t_outside)

# def thickness(dp, d_0, s_a, e_w, K):
#   t_wall = dp * d_0 / (2*s_a*e_w + 0.8*dp)
#   t_caps = dp * d_0 * K / (2*s_a*e_w + 2*dp*(K-0.1))
#   return t_wall, t_caps

#choose tank dimensions depending on constraints
#import the constraints for the Main file
#we maximize the volume of central tanks
#we assume circular cylinders with two hemisphere heads at ends

class Tank():
  def __init__(self,constraints,dp, s_a, e_w):
    """

    :param constraints: object of class constraints
    :param dp: pressure difference [bar]
    :param s_a: [MPa] allowable stress
    :param e_w: weld efficiency
    :param K:
    """
    self.d_0 = min(constraints.width,constraints.height)
    self.length = constraints.length #including the hemisphere ends
    self.dp = dp
    self.s_a = s_a
    self.e_w = e_w
    self.K = 1/6 * (2 + 1)
    self.outer_volume = m.pi * dp**2/4 * (self.length-self.d_0) + 4/3* m.pi * (self.d_0/2)**3 #cylinder volume + sphere volume


  def thickness(self): #D_Verstraete_Thesis_2009
    t_wall = self.dp * d_0 / (2 * s_a * e_w + 0.8 * self.dp)
    t_caps = self.dp * d_0 * self.K / (2 * s_a * e_w + 2 * self.dp * (self.K - 0.1))
    self.t_wall = t_wall
    self.t_caps = t_caps

  def tankwall_volume(self):
    self.inner_volume = m.pi * (dp/2-self.t_wall)**2 * self.length + 4/3* m.pi * (self.d_0/2-self.t_caps)**3 #volume which will be occupied by
                                                                                                             #LH2 and insulation
    self.wall_volume = self.outer_volume - self.inner_volume #this will determine the mass of the tank






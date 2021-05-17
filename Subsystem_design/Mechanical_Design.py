"""
This files contains the classes and functions for determining the mechanical design of the tank (thickness of the wall).
"""
import numpy as np
import matplotlib.pyplot as plt

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


K = 1/6 * (2 + )

def thickness(dp, d_0, s_a, e_w, K):
  t_wall = dp * d_0 / (2*s_a*e_w + 0.8*dp)
  t_caps = dp * d_0 * K / (2*s_a*e_w + 2*dp*(K-0.1))
  return t_wall, t_caps

#choose tank dimensions depending on constraints




"""
This files contains the classes and functions for determining the mechanical design of the tank (thickness of the wall).
"""
import numpy as np
import matplotlib.pyplot as plt

#Hydrogen tank storage data
p_out = 0.19028  # [bar] at 12100 [m] of altitude
p_tank = 3         # [bar]

t_out = 216.65   #[k]
t_tank = 20      #[k]


dp = np.abs(p_tank-p_out)
dt = np.abs(t_tank-t_tank)


def wall_thickness():
    s_a =




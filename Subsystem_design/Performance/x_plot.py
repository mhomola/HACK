from Subsystem_design.common_constants import Constants
import numpy as np
import matplotlib.pyplot as plt

class X_diagram(Constants):
    def __init__(self):
        super().__init__()
        self.VhV = 1                                    # tail/wing speed ratio
        self.x_ac_cruise = 0.0                          # aerodynamic center of the aircraft excluding the tail at cruise [MAC]
        self.x_ac_land = 0.0                            # aerodynamic center of the aircraft excluding the tail at landing [MAC]                   
        self.SM = 0.0
        self.ShS = 0.0
        self.ShLhSc = 0.0
        self.lh_over_c =  self.ShLhSc/self.ShS

      
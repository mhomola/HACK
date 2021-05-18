import numpy as np
import matplotlib.pyplot as plt

central_H2_vol = 15 # [m^3] the volume of LH2 in the central tanks area
aft_H2_vol = 30 #[m^3] the volume of LH2 that must be stored behind the passenger cabin

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
#ACT1 = spacial_constraint(length=,width=,height=)
#ACT1 = spacial_constraint(length=,width=,height=)
#aft = spacial_constraint(length=,width=,height=)
from common_constants import Constants
import numpy as np


class AerodynamicCharacteristics(Constants):
    def __init__(self):
        super().__init__()
        self.M = 0.78  # Cruising Mach number
        self.cruise_alt = 11280  # Initial cruising altitude [m]
        self.visc = 1.458 * 10**(-5)  # Air viscosity N*s/m^2
        



# Try out the class

if __name__ == '__main__':
    Ae = AerodynamicCharacteristics()

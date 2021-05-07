
#Tool to estimate for each configuration the added weight due to the hydrogen tanks.

import numpy as np
import Constants


def hydrogen_tank_mass(State):

    if State == 'Liquid':

        Grav_eff_L = np.arange(0.25,0.38,0.1)                      #Gravimetric efficiency if we choose liquid H2




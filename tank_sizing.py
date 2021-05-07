
#Tool to estimate for each configuration the added weight due to the hydrogen tanks.

import numpy as np
import Constants
import matplotlib.pyplot as plt


def hydrogen_tank_mass(State):

    if State == 'Liquid':

        Grav_eff_L = np.arange(0.25,0.38,0.1)                      #Gravimetric efficiency if we choose liquid H2







# Visulize some stuff about gas H2 tank given H2 volume as input
def gas_H_tanks(H2_vol):
    
    mass = []
    volume = []
    massH = []
    volumeH = []
    p_ran = np.arange(350,705,5)
    
    for p in range(350,705,5):
        H_mass_frac = 7 - 1/350 * p
        #d = 4 +17/350 * p
        d = -3.552714e-15 + 0.07845238*p - 0.00004285714*p**2 + 1.190476e-8*p**3
        H2_mass = H2_vol*0.001*d # code variable density
        tank_mass = (H2_mass*100)/H_mass_frac
        mass.append(tank_mass)
        massH.append(H2_mass)
        
    for p in range(350,705,5):
        H_vol_frac = 92 - 11/350 * p
        #d = 4 +17/350 * p
        d = -3.552714e-15 + 0.07845238*p - 0.00004285714*p**2 + 1.190476e-8*p**3
        tank_vol = (H2_vol*100)/H_vol_frac
        volume.append(tank_vol)
        volumeH.append(H2_vol*H_vol_frac)
         
    plt.figure()
    
    plt.subplot(411)
    plt.plot(p_ran, mass, color='tab:blue', marker='x',label='total mass')
    plt.xlabel("Pressure [bar]")
    plt.ylabel("Tank mass [kg]")
    plt.legend()
    plt.grid()
    
    plt.subplot(412)
    plt.plot(p_ran, massH, color='tab:red', marker='o',label='H2 mass')
    plt.xlabel("Pressure [bar]")
    plt.ylabel("Hydrogen mass [kg]")
    plt.legend()
    plt.grid()
    
    plt.subplot(413)
    plt.plot(p_ran, volume, color='tab:orange', marker='x',label='total volume')
    plt.xlabel("Pressure [bar]")
    plt.ylabel("Tank volume [l]")
    plt.legend()
    plt.grid()
    
    plt.subplot(414)
    plt.plot(p_ran, volumeH, color='tab:green', marker='o',label='H2 volume')
    plt.xlabel("Pressure [bar]")
    plt.ylabel("Hydrogen volume [l]")
    plt.legend()
    plt.grid()
    
    plt.show()
    
    return None

#Tool to estimate for each configuration the added weight due to the hydrogen tanks.

import numpy as np
from fuel_constants import *
from fuel_required import fuel_mass_calc
import matplotlib.pyplot as plt
import Constants
import matplotlib.pyplot as plt


def hydrogen_tank_mass(State):
    #This function takes the state of the hydrogen: 'liquid' or 'gas'; and computes the mass of the tank's mass
    #for different fuel ratios (ratio between kerosene and H2) and different gravimetric efficiencies
    if State == 'liquid':

        Grav_eff = np.arange(0.25,0.38,0.05)                            #Gravimetric efficiency if we choose liquid H2
        H2_mass = fuel_mass_calc(State=State,                           #Array of H2 masses for all
                                 d_k=k_d, d_LH2=LH2_d, d_GH2_g=GH2_d_g)[1]  #feasible volumes of H2
        System_mass = np.ones((len(H2_mass),len(Grav_eff)))
        Tank_mass = np.ones((len(H2_mass),len(Grav_eff)))

        for i in range(len(H2_mass)):
            for j in range(len(Grav_eff)):

                System_mass[i,j] = H2_mass[i]/Grav_eff[j]
                Tank_mass[i,j] = System_mass[i,j] - H2_mass[i]

    if State == 'gas':
        System_mass = 0
        H2_mass =  0
        Tank_mass = 0

    return System_mass,H2_mass, Tank_mass

def plotting_sys_mass(State):

    e_ratios = np.arange(0.1,0.3,0.001)
    Sys_mass,mass_H2,Mass_tank = hydrogen_tank_mass(State)

    fig1, ax1 = plt.subplots(3, 3,sharex = True,sharey = True)
    ax1[0,0].plot(e_ratios, Sys_mass[:,0])
    ax1[0,0].set_ylabel('Total Mass of System [kg]')
    ax1[0,0].set_title('Gravimetric coefficient = 0.25')

    ax1[0,1].plot(e_ratios, Sys_mass[:,1])
    ax1[0,1].set_title('Gravimetric coefficient = 0.30')

    ax1[0,2].plot(e_ratios, Sys_mass[:,2])
    ax1[0,2].set_title('Gravimetric coefficient = 0.35')

    ax1[1, 0].plot(e_ratios, mass_H2)
    ax1[1, 0].set_ylabel('H2 Mass [kg]')

    ax1[1, 1].plot(e_ratios, mass_H2)

    ax1[1, 2].plot(e_ratios, mass_H2)

    ax1[2, 0].plot(e_ratios, Mass_tank[:,0])
    ax1[2, 0].set_ylabel('Tank mass [kg]')
    ax1[2, 0].set_xlabel('Fuel Energy Ratio [-]')

    ax1[2, 1].plot(e_ratios, Mass_tank[:, 1])
    ax1[2, 1].set_xlabel('Fuel Energy Ratio [-]')

    ax1[2, 2].plot(e_ratios, Mass_tank[:, 2])
    ax1[2, 2].set_xlabel('Fuel Energy Ratio [-]')
    plt.show()

plotting_sys_mass(State='liquid')

def liquid_H_tanks(H2_vol):
    #input: volume available to store H2
    #output: tank mass, system mass and H2 mass
    Grav_eff = 0.25                                         # Gravimetric efficiency if we choose liquid H2
    H2_mass = H2_vol * 0.001 * LH2_d                        # Mass of Liquid hydrogen
    Tank_mass = H2_mass/Grav_eff                            # Tank's mass
    Tot_mass = H2_mass + Tank_mass                          # Total system mass

    return H2_mass,Tank_mass,Tot_mass







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

gas_H_tanks(41743.13)
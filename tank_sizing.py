
#Tool to estimate for each configuration the added weight due to the hydrogen tanks.

import numpy as np
from fuel_constants import *
from fuel_required import fuel_mass_calc
import matplotlib.pyplot as plt

def hydrogen_tank_mass(State):

    if State == 'liquid':

        Grav_eff = np.arange(0.25,0.38,0.05)                             #Gravimetric efficiency if we choose liquid H2
        H2_mass = fuel_mass_calc(State=State,                           #Array of H2 masses for all
                                 d_k=k_d, d_LH2=LH2_d, d_GH2=GH2_d)[1]  #feasible volumes of H2
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



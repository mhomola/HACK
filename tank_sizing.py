
#Tool to estimate for each configuration the added weight due to the hydrogen tanks.

import numpy as np
from fuel_constants import *
from Subsystem_design.fuel_required import fuel_mass_calc, fuel_volume_calc
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

    e_ratios = np.arange(0.11,0.5,0.001)
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
    Grav_eff = 0.5                                         # Gravimetric efficiency if we choose liquid H2
    H2_mass = H2_vol * 0.001 * LH2_d                        # Mass of Liquid hydrogen
    Tank_mass = H2_mass/Grav_eff                            # Tank's mass
    Tot_mass = H2_mass + Tank_mass                          # Total system mass

    return H2_mass,Tank_mass,Tot_mass

def liquid_H_tanks_2(Vol_avl):
    # input: volume available to store H2
    # output: tank mass, system mass and H2 mass
    Grav_eff = 0.5                                        # Gravimetric efficiency if we choose liquid H2
    rho_tank = 2266.248                                    # ASSUMED density of tank (based on density of
                                                           # insulation + wall, weighted average)
    H2_vol = (-Grav_eff * rho_tank / LH2_d * Vol_avl) / (-Grav_eff * rho_tank / LH2_d + Grav_eff - 1)
    H2_mass = H2_vol * LH2_d * 0.001
    Tank_mass = H2_mass/Grav_eff
    Tot_mass = H2_mass + Tank_mass

    return H2_mass,Tank_mass,Tot_mass

def compare_des_liquid():

    #Cargo + Optional Tanks
    H2_mass_1, Tank_mass_1, Tot_mass_1 = liquid_H_tanks(11640)
    #Raising aisle
    H2_mass_2, Tank_mass_2, Tot_mass_2 = 0,0,0
    #A321
    H2_mass_3, Tank_mass_3, Tot_mass_3 = liquid_H_tanks(37261.17)
    #Flat Bottom
    H2_mass_4, Tank_mass_4, Tot_mass_4 = liquid_H_tanks_2(4875)#4875 3250
    #Wing podded
    H2_mass_5, Tank_mass_5, Tot_mass_5 = liquid_H_tanks(37261.17)
    #Beluga
    H2_mass_6, Tank_mass_6, Tot_mass_6 = liquid_H_tanks(37261.17)

    H2Mass_array = np.array([H2_mass_1,H2_mass_2,H2_mass_3,H2_mass_4,H2_mass_5,H2_mass_6])
    TMass_array = np.array([Tank_mass_1, Tank_mass_2, Tank_mass_3, Tank_mass_4, Tank_mass_5, Tank_mass_6])
    TotMass_array = np.array([Tot_mass_1, Tot_mass_2, Tot_mass_3, Tot_mass_4, Tot_mass_5, Tot_mass_6])
    conf = np.arange(1,7)
    x = np.arange(len(conf))                                     # Label locations
    width = 0.25                                                    # Width of the bars

    plt.bar(x - width, H2Mass_array, color='tab:blue', width=0.25)
    plt.bar(x , TMass_array, color='tab:orange', width=0.25)
    plt.bar(x + width, TotMass_array, color='tab:red', width=0.25)
    plt.ylabel('Mass [kg]')
    plt.title('Design options Added weight')
    plt.xticks(x, ('Cargo', 'Raisle', 'A321', 'FlBttm', 'Wpodded','Beluga'))
    plt.legend(labels=['Hydrogen mass', 'Tank mass','Total mass'])
    plt.show()
compare_des_liquid()

# Visulize some stuff about gas H2 tank given H2 volume as input
def gas_H_tanks():

    mass = []
    volume = []
    massH = []
    volumeH = []
    p_ran = np.arange(350,705,5)

    for p in p_ran:
        H_vol_frac = 92 - 11 / 350 * p
        H_mass_frac = 7 - 1/350 * p
        d = -3.552714e-15 + 0.07845238*p - 0.00004285714*p**2 + 1.190476e-8*p**3
        H2_vol = fuel_volume_calc(d_LH2=LH2_d, d_GH2= GH2_d, d_GH2_g= d, d_k=k_d, Ed_H2=H2_ed, Ed_k=k_ed, tot_vol_k=fuel_capacity_a320neo,
                             e_ratio=0.3,state='gas')[1]
        H2_mass = H2_vol*0.001*d
        tank_mass = (H2_mass*100)/H_mass_frac
        tank_vol = (H2_vol * 100) / H_vol_frac
        volH = (H_vol_frac / 100) * tank_vol
        mass.append(tank_mass)
        massH.append(H2_mass)
        volume.append(tank_vol)
        volumeH.append(volH)


    plt.figure()

    plt.subplot(411)
    plt.plot(p_ran, mass, color='tab:blue', marker='x',label='total mass')
    plt.xlabel("Pressure [bar]")
    plt.ylabel("Total mass [kg]")
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
    plt.ylabel("Total volume [l]")
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

# Visualize stuff for the gas tanks given the max available volume

def gas_H_tanks_VOL(tot_Vol):

    mass = []
    volume = []
    massH = []
    volumeH = []
    p_ran = np.arange(350,705,5)

    for p in p_ran:
       H_vol_frac = 92 - 11/350 * p
       H_vol = H_vol_frac/100 * tot_Vol
       d = -3.552714e-15 + 0.07845238*p - 0.00004285714*p**2 + 1.190476e-8*p**3
       H_mass = H_vol*0.001*d

       H_mass_frac = 7 - 1/350 * p
       tank_mass = (H_mass*100)/H_mass_frac

       mass.append(tank_mass)
       volume.append(tot_Vol)
       massH.append(H_mass)
       volumeH.append(H_vol)

    plt.figure()

    plt.subplot(411)
    plt.plot(p_ran, mass, color='tab:blue', marker='x',label='total mass')
    plt.xlabel("Pressure [bar]")
    plt.ylabel("Total mass [kg]")
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
    plt.ylabel("Total volume [l]")
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

#Option 1: cargo + optional tanks
#gas_H_tanks_VOL(43220)
#Option 2: raising aisle
#gas_H_tanks_VOL(0)
#Option 3: A321
#gas_H_tanks_VOL(263400)
#Option 4: flat bottom (upper bound)
#gas_H_tanks_VOL(4875)
#Option 5: wing podded
gas_H_tanks()
#Option 6: beluga
gas_H_tanks()

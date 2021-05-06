#File to plot kerosene amount vs hydrogen for different fuel ratio mixtures in the the engine

import matplotlib.pyplot as plt
#import matplotlib.axes as ax
import numpy as np
import fuel_constants

def kerosene_vol():

    return

def hydrogen_vol():

    return

def plotting():#K_mass,K_vol,H2_mass,H2_vol):

    E_ratio = np.arange(10,60,10).astype(str)
    k_mass = np.array([1,20,5,4,6])
    h2_mass = np.array([50,9,8,4,5])
    k_vol = np.array([1, 20, 5, 4, 6])
    h2_vol = np.array([50, 9, 8, 4, 5])
    x = np.arange(len(E_ratio))             #Label locations
    width = 0.25                            #Width of the bars

    fig, ax = plt.subplots()

    plt.subplot(121)
    k_bar_m = plt.bar(x - width/2,k_mass,width, label = 'Kerosene Mass')
    h2_bar_m = plt.bar(x + width/2,h2_mass,width, label = 'Hydrogen Mass')
    plt.ylabel('Fuel mass [kg]')
    plt.xlabel('Energy ratio [%]')
    ax.set_xticks(x)
    ax.set_xticklabels(E_ratio)

    plt.legend('best')

    plt.subplot(122)
    k_bar_vol = plt.bar(x - width / 2, k_vol, width, label='Kerosene Volume')
    h2_bar_vol = plt.bar(x + width / 2, h2_vol, width, label='Hydrogen Volume')
    plt.ylabel('Fuel volume [$m^3$]')
    plt.xlabel('Energy ratio [%]')
    ax.set_xticks(x)
    ax.set_xticklabels(E_ratio)
    plt.legend('best')

    plt.show()

def range():


import numpy as np
import matplotlib.pyplot as plt
from fuel_constants import *

def fuel_volume_calc(d_LH2, d_GH2,d_GH2_g, d_k, Ed_H2, Ed_k, tot_vol_k, e_ratio, state,y = 0.95):
    """

    :param d_LH2: density of liquid hydrogen [kg/m^3]
    :param d_LH2: density of gaseous hydrogen [kg/m^3]
    :param d_k: density of kerosene [kg/m^3]
    :param Ed_H2: energy density of hydrogen [kWh/kg]
    :param Ed_k: energy density of kerosene [kWh/kg]
    :param tot_vol_k: volume of kerosene for original A320neo [l]
    :param f_ratio: fuel ratio in terms of energy E_H2/E_k
    :return: V_K [l], the volume of kerosene needed; V_H2 [l], the volume of hydrogen needed
    """

    E_tot = Ed_k * d_k / 1000 * tot_vol_k

    r = 1 / (1 + e_ratio)

    E_k = E_tot * r
    E_H2 = e_ratio * E_k

    V_k = E_k / Ed_k / d_k * 1000

    if state == 'gas':

        V_H2 = (E_H2/Ed_H2)/d_GH2_g *1000

    if state == 'liquid':

        V_LH2 = (E_H2/Ed_H2)*(1/(d_LH2+d_GH2*(1-y)/y))
        V_H2 = V_LH2*(1+(1-y)/y)*1000

    return V_k, V_H2


def fuel_mass_calc(State, d_k, d_LH2,d_GH2_g):

    #Function to compute mass of kerosene and H2, both for the case in which the hydrogen is stored
    #in liquid state, as well as when it is stored in gas (compressed) state.


    if State == 'liquid':
        Vk, VH2 = fuel_volume_calc(d_LH2=LH2_d, d_GH2=GH2_d,d_GH2_g=GH2_d_g, d_k=k_d, # Volume of kerosene and Liquid H2[l]
                                         Ed_H2=H2_ed, Ed_k=k_ed,
                                         tot_vol_k=fuel_capacity_a320neo,
                                         e_ratio=e_ratios, state=State)
        V_tot = Vk + VH2                                                              #Total volume [l]

        m_k, m_H2 = Vk * d_k * 0.001, VH2 * d_LH2 * 0.001                             #Mass of Kerosene and Liquid H2 [kg]
        m_tot = m_k + m_H2                                                            #Total mass

    if State == 'gas':
        Vk, VH2 = fuel_volume_calc(d_LH2=LH2_d, d_GH2=GH2_d,d_GH2_g=GH2_d_g, d_k=k_d,  # Volume of kerosene and Liquid H2[l]
                                         Ed_H2=H2_ed, Ed_k=k_ed,
                                         tot_vol_k=fuel_capacity_a320neo,
                                         e_ratio=e_ratios, state=State)
        V_tot = Vk + VH2                                               #Total volume [l]

        m_k, m_H2 = Vk * d_k * 0.001, VH2 * d_GH2_g * 0.001              #Mass of Kerosene and Gas H2 [kg]
        m_tot = m_k + m_H2                                               #Total mass

    return m_k, m_H2, m_tot, Vk, VH2,V_tot

def plotting_vol_mass():


    m_k_l, m_H2_l,m_tot_l,V_k_l, V_H2_l, V_tot_l = fuel_mass_calc(State='liquid',d_k=k_d,d_LH2=LH2_d,d_GH2_g=GH2_d_g)
    m_k_g, m_H2_g, m_tot_g, V_k_g, V_H2_g, V_tot_g = fuel_mass_calc(State='gas',d_k=k_d,d_LH2=LH2_d,d_GH2_g=GH2_d_g)

    fig1, ax1 = plt.subplots(3,2, sharex = True)
    ax1[0,0].plot(e_ratios, V_H2_l)
    ax1[0,0].set_ylabel('Liquid H2 Volume [l]')

    ax1[1,0].plot(e_ratios, V_k_l)
    ax1[1,0].set_ylabel('Kerosene Volume [l]')

    ax1[2,0].plot(e_ratios, V_tot_l)
    ax1[2,0].set_xlabel('Fuel Energy Ratio [-]')
    ax1[2,0].set_ylabel('Total Fuel Volume [l]')

    ax1[0, 1].plot(e_ratios, V_H2_g)
    ax1[0, 1].set_ylabel('Gas H2 Volume [l]')

    ax1[1, 1].plot(e_ratios, V_k_g)
    ax1[1, 1].set_ylabel('Kerosene Volume [l]')

    ax1[2, 1].plot(e_ratios, V_tot_g)
    ax1[2, 1].set_xlabel('Fuel Energy Ratio [-]')
    ax1[2, 1].set_ylabel('Total Fuel Volume [l]')

    plt.show()
    fig2, ax2 = plt.subplots(3, 2, sharex=True)
    ax2[0, 0].plot(e_ratios, m_H2_l)
    ax2[0, 0].set_ylabel('Liquid H2 Mass [kg]')

    ax2[1, 0].plot(e_ratios, m_k_l)
    ax2[1, 0].set_ylabel('Kerosene Mass [kg]')

    ax2[2, 0].plot(e_ratios, m_tot_l)
    ax2[2, 0].set_xlabel('Fuel Energy Ratio [-]')
    ax2[2, 0].set_ylabel('Total Fuel Mass [kg]')

    ax2[0, 1].plot(e_ratios, m_H2_g)
    ax2[0, 1].set_ylabel('Gas H2 Mass [kg]')

    ax2[1, 1].plot(e_ratios, m_k_g)
    ax2[1, 1].set_ylabel('Kerosene Mass [kg]')

    ax2[2, 1].plot(e_ratios, m_tot_g)
    ax2[2, 1].set_xlabel('Fuel Energy Ratio [-]')
    ax2[2, 1].set_ylabel('Total Fuel Mass [kg]')

    plt.show()

def VH2_calc():
    e_ratios = np.arange(0.1,0.5,0.001)
    return fuel_volume_calc(d_LH2=LH2_d, d_GH2= GH2_d, d_GH2_g= GH2_d_g,d_k=k_d, Ed_H2=H2_ed, Ed_k=k_ed
                            , tot_vol_k=fuel_capacity_a320neo,e_ratio=e_ratios,state='liquid')[1]

e_ratios = np.arange(0.1,0.5,0.001) #ratio of energy stored in H2 to kerosene (NOT GLOBAL!)

plotting_vol_mass()

# Compute volumes

E_ratio = 0.5 # [-] #ratio of energy stored in H2 to kerosene (NOT GLOBAL!)

V_k, V_H2 = fuel_volume_calc(d_LH2=LH2_d, d_GH2= GH2_d, d_GH2_g= GH2_d_g, d_k=k_d, Ed_H2=H2_ed, Ed_k=k_ed, tot_vol_k=fuel_capacity_a320neo,
                             e_ratio=E_ratio,state='liquid')

print('Volume of H2 is ', V_H2, ' l')
print('Volume of kerosene is ', V_k, ' l')



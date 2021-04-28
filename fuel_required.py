import numpy as np
from fuel_constants import *

def fuel_volume_calc(d_H2, d_k, Ed_H2, Ed_k, tot_vol_k, f_ratio):
    """

    :param d_H2: density of hydrogen [kg/m^3]
    :param d_k: density of kerosene [kg/m^3]
    :param Ed_H2: energy density of hydrogen [kWh/kg]
    :param Ed_k: energy density of kerosene [kWh/kg]
    :param tot_vol_k: volume of kerosene for original A320neo [l]
    :param f_ratio: fuel ratio in terms of energy E_H2/E_k
    :return: V_K [l], the volume of kerosene needed; V_H2 [l], the volume of hydrogen needed
    """

    E_tot = Ed_k * d_k / 1000 * tot_vol_k
    E_k = (1 - f_ratio) * E_tot
    E_H2 = f_ratio * E_tot

    V_k = E_k / Ed_k / d_k * 1000
    V_H2 = E_H2 / Ed_H2 / d_H2 * 1000

    return V_k, V_H2

# Compute volumes

Fuel_ratio = 0.2 # [-]
fuel_capacity_a320neo = 20000 # [l]
V_k, V_H2 = fuel_volume_calc(d_H2=LH2_d, d_k=k_d, Ed_H2=H2_ed, Ed_k=k_ed, tot_vol_k=fuel_capacity_a320neo,
                             f_ratio=Fuel_ratio)

print('Volume of LH2 is ', V_H2, ' l')
print('Volume of kerosene is ', V_k, ' l')

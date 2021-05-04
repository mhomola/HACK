import numpy as np
import matplotlib.pyplot as plt
from fuel_constants import *

def fuel_volume_calc(d_LH2, d_GH2, d_k, Ed_H2, Ed_k, tot_vol_k, e_ratio, y = 0.95):
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

    r = 1/(1+e_ratio)

    E_k = E_tot*r
    E_H2 = e_ratio * E_k

    V_k = E_k / Ed_k / d_k * 1000

    V_LH2 = (E_H2/Ed_H2)*(1/(d_LH2+d_GH2*(1-y)/y))
    V_H2 = V_LH2*(1+(1-y)/y)*1000

    return V_k, V_H2

# Compute volumes

E_ratio = 0.15 # [-]
fuel_capacity_a320neo = 26730 # [l]
V_k, V_H2 = fuel_volume_calc(d_LH2=LH2_d, d_GH2= GH2_d, d_k=k_d, Ed_H2=H2_ed, Ed_k=k_ed, tot_vol_k=fuel_capacity_a320neo,
                             e_ratio=E_ratio)

print('Volume of H2 is ', V_H2, ' l')
print('Volume of kerosene is ', V_k, ' l')

e_ratios = np.arange(0.1,0.3,0.001)
V_k, V_H2 = fuel_volume_calc(d_LH2=LH2_d, d_GH2= GH2_d, d_k=k_d, Ed_H2=H2_ed, Ed_k=k_ed, tot_vol_k=fuel_capacity_a320neo,
                             e_ratio=e_ratios)
V_tot = V_k+V_H2

f1 = plt.figure()
ax1 = f1.add_subplot(311)
ax1.plot(e_ratios,V_H2)
ax1.set_ylabel('H2 Volume [l]')
ax2 = f1.add_subplot(312)
ax2.plot(e_ratios, V_k)
ax2.set_ylabel('Kerosene Volume [l]')
ax3 = f1.add_subplot(313)
ax3.plot(e_ratios, V_tot)
ax3.set_xlabel('Fuel Energy Ratio [-]')
ax3.set_ylabel('Total Fuel Volume [l]')
plt.show()
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:53:39 2021

@author: sarar
"""

"""
    SOME NOTES:

    mr = mass ratio
    mf = mass flow
    _cc = during the combustion, so only the reactants
    _mix = after the cooling air is added

"""

import math as m
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

''' DEFINITION OF THE FUNCTION '''
def cp_regression(data, T):
    T_before = data[ np.where(data[:,0] < T)[0][-1] ][0]
    T_after = data[ np.where(data[:,0] > T)[0][0] ][0]
    cp_before = data[ np.where(data[:,0] < T)[0][-1] ][1]
    cp_after = data[ np.where(data[:,0] > T)[0][0] ][1]
    
    slope = (cp_after - cp_before) / (T_after - T_before)
    intersection = cp_after - slope*T_after
    
    cp_T = slope*T + intersection   
    
    index_after = np.where(data[:,0] > T)[0][0]
    
    return cp_T, index_after  

def cp_temperature(data, T):
    for i in range(len(data)):
        if data[i][0] == T:
            cp = data[i][0]
            index_after = i+1
    
    return cp, index_after

def cp_between(data, i0, cp_array, T_array, T, T_max):
    # input T is only to initialise loop
    for i in range(i0, len(data)):
        if T < T_max:
            cp_array = np.append(cp_array, data[i][1])
            T_array = np.append(T_array, data[i][0])
            T = data[i][0] 
        else:
            break
    
    return cp_array, T_array

def cp_first_last(data, T, cp_array, T_array):
    
    if not(T in data[:][0]):
        cp, index_after = cp_regression(data, T)
    else:
        cp, index_after = cp_temperature(data, T)
    
    cp_array = np.append(cp_array, cp)
    T_array = np.append(T_array, T)
       
    return cp_array, T_array, index_after

def h(h0, cp_array, T_array):
    cp_integral = np.array([])
    for i in range(len(cp_array)):
        cp_integral = np.append(cp_integral, cp_array[i]*T_array[i])
    
    h = h0 + np.sum(cp_integral) 
    
    return h


''' BEGINNING OF CODE '''
#print('hello')
N2_cp_data = np.array(np.genfromtxt('N2_cp.dat'))
#print(N2_cp_data)
h2_cp_data = np.array(np.genfromtxt('C:\\Users\\sarar\\OneDrive\\Ambiente de Trabalho\\Folders\\DELFT\\3rd year\\DSE\\h2_cp.dat'))

Tair = 800 # [K] # Temperature of air when injected (both in PZ and SZ)
Tcc = 2000 # [K] # Temperature of fuel after the combustion, i.e. maximum temperature
T0 = 298 # [K] --> 15 degrees C, where does it come from? At sea level? # To be used in the formulae

# Mass flows
# Before combustion
mf_h2_cc = 1
mf_ker_cc = 1
mf_air_cc = 1
mf_cc = mf_h2_cc + mf_ker_cc + mf_air_cc
# After cooling air is injected
mf_cool = 1
mf_end = mf_cc + mf_cool

# Mass ratios
# Before combustion
mr_h2_cc = mf_h2_cc / mf_cc
mr_ker_cc = mf_ker_cc / mf_cc
mr_air_cc = mf_air_cc / mf_cc
# After cooling air is injected
mr_h2_mix = mf_h2_cc / mf_end
mr_ker_mix = mf_ker_cc / mf_end
mr_air_mix = (mf_air_cc + mf_cool) / mf_end
# When we mix reactants to air to cool engine
mr_cc = mf_cc / mf_end
mr_cool = mf_cool/ mf_end


# Find h_cc: find h_ker, h_h2, h_air [J/kg]
h_decane = -249.9*1000 + (Tcc-T0)*235 # [J/mol] --> these values are at 25C, is this okay?
h_benzene = 8*1000 + (Tcc-T0)*154 # [J/mol]
# h_hexane =

h_ker_cc = 0.74/(0.74+0.15)*h_decane + 0.15/(0.74+0.15)*h_benzene # [J/mol] --> Cp not found for C8H19
                                        # can I use this formula? Or maybe not dividing by the sum
# Find h_h2_cc [J/kg]
cp_h2_cc, T_h2_cc = np.array([]), np.array([])
# Find initial data to retrieve
cp_h2_cc, T_h2_cc, index_after = cp_first_last(h2_cp_data, T0, cp_h2_cc, T_h2_cc)
# Find all other cp's and temperatures
cp_h2_cc, T_h2_cc = cp_between(h2_cp_data, index_after, cp_h2_cc, T_h2_cc, T0, Tcc)   
# Find final data
if T_h2_cc[-1] != Tcc:
    cp_h2_cc, T_h2_cc, index_after = cp_first_last(h2_cp_data, Tcc, cp_h2_cc, T_h2_cc)

print("\nHydrogen\nTempreature [K]\n", T_h2_cc)
print("Specific heat [kJ/(kg K)]\n",cp_h2_cc)
    
h_h2_cc = h(0, cp_h2_cc, T_h2_cc)*1000 # [J/kg]

# Find h_air_cc [J/kg]
cp_air_cc, T_air_cc = np.array([]), np.array([])
# Find initial data to retrieve
cp_air_cc, T_air_cc, index_after = cp_first_last(N2_cp_data, T0, cp_air_cc, T_air_cc)
# Find all other cp's and temperatures
cp_air_cc, T_air_cc = cp_between(N2_cp_data, index_after, cp_air_cc, T_air_cc, T0, Tcc)   
# Find final data
if T_h2_cc[-1] != Tcc:
    cp_air_cc, T_air_cc, index_after = cp_first_last(N2_cp_data, Tcc, cp_air_cc, T_air_cc)

print("\nAir\nTempreature [K]\n", T_air_cc)
print("Specific heat [kJ/(kg K)]\n",cp_air_cc)

h_air_cc = h(0, cp_air_cc, T_air_cc)*1000 # [J/kg]
        
# Total enthalpy in combustion
h_cc = h_ker_cc*mr_ker_cc + h_h2_cc*mr_h2_cc + h_air_cc*mr_air_cc



# Find h_cool [J/kg]
cp_cool, T_cool = np.array([]), np.array([])
# Find initial data to retrieve
cp_cool, T_cool, index_after = cp_first_last(N2_cp_data, T0, cp_cool, T_cool)
# Find all other cp's and temperatures
cp_cool, T_cool = cp_between(N2_cp_data, index_after, cp_cool, T_cool, T0, Tair)   
# Find final data
if T_cool[-1] != Tair:
    cp_cool, T_cool, index_after = cp_first_last(N2_cp_data, Tair, cp_cool, T_cool)

print("\nAir to cool\nTempreature [K]\n", T_cool)
print("Specific heat [kJ/(kg K)]\n",cp_cool)

h_cool = h(0, cp_cool, T_cool) # [kJ/kg]


# Find h_mix
h_mix = mr_cc*h_cc + mr_cool*h_cool

# Find delta_h_mix (only contribution comes from kerosene, but now mr_ker changes because of the cooling air added)


# Find cp_mix
# Assume a temperature at the end of the cc: T_assumed
# Find cp_mix for that temperature: C_p_mix(T_assumed) = SUM( c_p_k(T_assumed)*mr_k )
# Find the temperature at the end of the combustion chamber
T_end = T0 + (h_mix - delta_h_mix) / Cp_mix
# Check if T_end ~ T_assumed, if not iterate T_assumed

# Find the temperature in case we case delta Cp is sufficiently small
T_linear = T_cc*mr_cc + T_cool*mr_cool

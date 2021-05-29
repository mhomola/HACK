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


def h(h0, data, T0, Tmax):    
    cp_array, T_array = np.array([]), np.array([])
    # Find initial data to retrieve
    cp_array, T_array, index_after = cp_first_last(data, T0, cp_array, T_array)
    # Find all other cp's and temperatures
    cp_array, T_array = cp_between(data, index_after, cp_array, T_array, T0, Tmax)   
    # Find final data
    if T_array[-1] != Tmax and data[-1][1] >= Tmax:
        cp_array, T_array, index_after = cp_first_last(data, Tmax, cp_array, T_array)     
    
    cp_integral = np.array([])
    for i in range(len(cp_array)):
        cp_integral = np.append(cp_integral, cp_array[i]*T_array[i])
    
    h = h0 + np.sum(cp_integral)
    
    return h


''' BEGINNING OF CODE '''
N2_cp_data = np.array(np.genfromtxt('N2_cp.dat')) # T[K]; cp[kJ/(kg*K)]
h2_cp_data = np.array(np.genfromtxt('h2_cp.dat')) # T[K]; cp[kJ/(kg*K)]
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C124185&Units=SI&Mask=7
C10H22_cp_data = np.array(np.genfromtxt('C10H22_cp.dat'))  # T[K]; cp[J/(mol*K)]
h0_C10H22 = -249.7 # [kJ/mol]
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C95636&Units=SI&Mask=7
C9H12_cp_data = np.array(np.genfromtxt('C9H12_cp.dat'))  # T[K]; cp[J/(mol*K)]
h0_C9H12 = -13.9 # [kJ/mol]


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
h_decane = h(h0_C10H22*1000, C10H22_cp_data, T0, Tcc) # [J/mol] C10H22
h_benzene = h(h0_C9H12*1000, C9H12_cp_data, T0, Tcc) # [J/mol] C9H12
# h_hexane: couldn't find values

h_ker_cc = 0.74/(0.74+0.15)*h_decane + 0.15/(0.74+0.15)*h_benzene # [J/mol] --> Cp not found for C8H19
                                        # can I use this formula? Or maybe not dividing by the sum
# Find h_h2_cc [J/kg]
h_h2_cc = h(0, h2_cp_data, T0, Tcc)*1000 # [J/kg]

# Find h_air_cc [J/kg]
h_air_cc = h(0, h2_cp_data, T0, Tcc)*1000 # [J/kg]
 
# Total enthalpy in combustion
h_cc = h_ker_cc*mr_ker_cc + h_h2_cc*mr_h2_cc + h_air_cc*mr_air_cc

# Find h_cool [J/kg]
h_cool = h(0, h2_cp_data, T0, Tair)*1000 # [J/kg]

# Find h_mix
h_mix = mr_cc*h_cc + mr_cool*h_cool

# Find delta_h_mix (only contribution comes from kerosene, but now mr_ker changes because of the cooling air added)
r_h0_C9H12 = 1
r_C10H22 = 1
delta_h_mix = mr_ker_mix*r_h0_C9H12*h0_C10H22+mr_ker_mix*r_C10H22*h0_C9H12

# Find cp_mix
# Assume a temperature at the end of the cc: T_assumed
# Find cp_mix for that temperature: C_p_mix(T_assumed) = SUM( c_p_k(T_assumed)*mr_k )
cp_mix = mr_air_mix*cp_air_mix + mr_h2_mix*cp_h2_mix + mr_ker_mix*cp_ker_mix
# Find the temperature at the end of the combustion chamber
T_end = T0 + (h_mix - delta_h_mix) / cp_mix
# Check if T_end ~ T_assumed, if not iterate T_assumed

# Find the temperature in case delta Cp is sufficiently small
T_linear = T_cc*mr_cc + T_cool*mr_cool

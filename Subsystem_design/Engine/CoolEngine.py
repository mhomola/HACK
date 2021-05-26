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

print('hello')
N2_cp_data = np.array(np.genfromtxt('N2_cp.dat'))
print(N2_cp_data)
h2_cp_data = np.array(np.genfromtxt('C:\\Users\\sarar\\OneDrive\\Ambiente de Trabalho\\Folders\\DELFT\\3rd year\\DSE\\h2_cp.dat'))

Tair = 800 # [K] # Temperature of air when injected (both in PZ and SZ)
Tcc = 2000 # [K] # Temperature of fuel after the combustion, i.e. maximum temperature
T0 = 298 # [K] --> 15 degrees C, where does it come from? At sea level? # To be used in the formulae

# Mass flows
# Before combustion
mf_h2_cc = 0
mf_ker_cc = 0
mf_air_cc = 0
mf_cc = mf_h2_cc + mf_ker_cc + mf_air_cc
# After cooling air is injected
mf_cool = 0
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


# Find h_cc: find h_ker, h_h2, h_air [J/kg]
h_decane = -249.9*1000 + (Tcc-T0)*235 # [J/mol] --> these values are at 25C, is this okay?
h_benzene = 8*1000 + (Tcc-T0)*154 # [J/mol]
# h_hexane =

h_ker = 0.74/(0.74+0.15)*h_decane + 0.15/(0.74+0.15)*h_benzene # [J/mol] --> Cp not found for C8H19
                                        # can I use this formula? Or maybe not dividing by the sum

h_h2 = (Tcc-T0)

def cp_regression(data, T):
    T_before = data[ np.where(data[:,0] < T)[0][-1] ][0]
    T_after = data[ np.where(data[:,0] > T)[0][0] ][0]
    cp_before = data[ np.where(data[:,0] < T)[0][-1] ][1]
    cp_after = N2_cp_data[ np.where(data[:,0] > T)[0][0] ][1]
    
    slope = (cp_after - cp_before) / (T_after - T_before)
    intersection = cp_after - slope*T_after
    
    cp_T = slope*T + intersection   
    
    index_after = np.where(data[:,0] > T)[0][0]
    
    return cp_T, index_after  

def cp_temperature(data, T):
    for i in range(len(data)):
        if data[i][0] == T:
            cp = N2_cp_data[i][0]
            index_after = i+1
    
    return cp, index_after

def check_temp(T_assumed):
    if T_assumed - 10 <= T_end <= T_assumed + 10:
    print("T_assumed = T_end", T_assumed, T_end)

cp_air_cc = np.array([])
T_air_cc = np.array([])

# Find initial data to retrieve
if not(T0 in N2_cp_data[:][0]):
    cp_T0, index_after = cp_regression(N2_cp_data, T0)
else:
    cp_T0, index_after = cp_temperature(N2_cp_data, T0)

cp_air_cc = np.append(cp_air_cc, cp_T0)
T_air_cc = np.append(T_air_cc, T0)


# Find all other cp's and temperatures
for i in range(len(N2_cp_data)):
    while N2_cp_data[i][0] < Tcc:
        
        
    
        
        
    
# make it a sum of (T*cp)i between Tcc and T0. If the
                #initial or final temperatures are not in the file, interpolate between before and after

h_air = (Tcc-T0) # approximate to air to only N2
        # make it a sum of (T*cp)i between Tcc and T0. If the initial or final
        # temperatures are not in teh file, interpolate between before and after

h_cc = h_ker*mr_ker_cc + h_h2*mr_h2_cc + h_air*mr_air_cc

# Find h_cool [J/kg]
h_cool = (T_air-T0)

#for i in range(len(N2_cp_data))

# Find h_mix
h_mix = mr_cc*h_cc + mr_cool*h_cool

# Find delta_h_mix (only contribution comes from kerosene, but now mr_ker changes because of the cooling air added)


# Find cp_mix
# Assume a temperature at the end of the cc: T_assumed
T_assumed = 2000 #K
# Find cp_mix for that temperature: C_p_mix(T_assumed) = SUM( c_p_k(T_assumed)*mr_k )
cp_mix = cp_regression(h2_cp_data,T_assumed) * mr_h2_mix + cp_regression(N2_cp_data,T_assumed) * mr_air_mix + cp_ker * \
         (T_assumed - T_0) * mr_ker_mix     #cp_ker does not have a list with values yet. will fix that later

# Find the temperature at the end of the combustion chamber
T_end = T0 + (h_mix - delta_h_mix) / cp_mix
# Check if T_end ~ T_assumed, if not iterate T_assumed


#blabla this is a check


# Find the temperature in case we case delta Cp is sufficiently small
T_linear = T_cc*mr_cc + T_cool*mr_cool

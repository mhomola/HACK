# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 09:39:04 2021

@author: sarar
"""

# TODO
# Find estimates for all constants
# Find altitude and speed of each phase
# Determine stoichiometric ratio --> find equivalence ratio
# Fuel flow of each fuel

# DONE
# Determine LHV_fuel according to mass ratio of fuels
# Determine atmospheric conditions


from Subsystem_design.common_constants import Constants

import math as m
import numpy as np
import pandas as pd
# from atmosphere import atmos_prop



        

    
''' CONSTANTS '''
# D_fan = 4 # [m]
# A_fan = m.pi * D_fan**2 / 4 # [m2]

# LHV_h2 = 120 # [MJ/kg]
# LHV_ker = 43 # [MJ/kg]

#aircraft = "HACK"
constants = Constants()
print(constants.R)

aircraft = input('Would you like to do the analysis for the engine of A320neo or A320-HACK? Answer neo or HACK: ')

# A320neo
if aircraft == "neo":
    neo_data = engine_data_neo()
    eta_inlet = neo_data.eta_inlet
    PR_fan = neo_data.PR_fan
    eta_fan = neo_data.eta_fan
    BR = neo_data.BR
    eta_LPC = neo_data.eta_LPC # or "booster"
    eta_HPC = neo_data.eta_HPC
    eta_LPT = neo_data.eta_LPT
    eta_HPT = neo_data.eta_HPT
    eta_mech = neo_dataeta_mech
    eta_cc = neo_data.eta_cc
    PR_LPC = neo_data.PR_LPC # or "booster"
    PR_HPC = neo_data.PR_HPC
    eta_nozzle = neo_data.eta_nozzle
    PR_cc = neo_data.PR_cc
    T04 = neo_data.T04 # [K] Temperature at the end of the combustion chamber
    LHV_f = neo_data.LHV_f

# A320-HACK
if aircraft == "HACK":
    hack_data = engine_data_hack()
    eta_inlet = hack_data.eta_inlet
    PR_fan = hack_data.PR_fan
    eta_fan = hack_data.eta_fan
    BR = hack_data.BR
    eta_LPC = hack_data.eta_LPC # or "booster"
    eta_HPC = hack_data.eta_HPC
    eta_LPT = hack_data.eta_LPT
    eta_HPT = hack_data.eta_HPT
    eta_mech = hack_dataeta_mech
    eta_cc = hack_data.eta_cc
    PR_LPC = hack_data.PR_LPC # or "booster"
    PR_HPC = hack_data.PR_HPC
    eta_nozzle = hack_data.eta_nozzle
    PR_cc = hack_data.PR_cc
    T04 = hack_data.T04 # [K] Temperature at the end of the combustion chamber
    LHV_f = hack_data.LHV_f
    

    




''' VARIABLES THAT VARY WITH FLIGHT PHASE '''
phases = np.array(['idle', 'taxi out', 'takeoff', 'climb', 'cruise', 'approach', 'taxi in'])

# Speed and mass of air
M0 = np.array([ 0.2, 0.2, 0.5, 0.5, 0.78, 0.5, 0.2]) # [-] Mach number

h = np.array([ 10, 10, 50, 3000, 11280, 3000, 10]) # [m] altitude
# # calculate density and temperature at that altitude
# T0, p0, rho0 = np.array([]), np.array([]), np.array([])
# for altitude in h:
#     T, p, rho = atmos_prop(altitude)
#     T0 = np.append(T0, T)
#     p0 = np.append(p0, p)
#     rho0 = np.append(rho0, rho)
    
    
# rho0 = np.array([ 1, 1, 1, 1, 1, 1, 1, 1]) # [kg/m3] ambient density
# T0 = np.array([ 1, 1, 1, 1, 1, 1, 1, 1]) # [K] ambient temperature
# p0 = np.array([ 1, 1, 1, 1, 1, 1, 1, 1]) # [Pa] ambient pressure
#a = np.array([ 1, 1, 1, 1, 1, 1, 1, 1]) # [m/s] speed of air

a0 = np.sqrt( 1.4 * 8.314 * T0 ) # [m/s] | 8.314 = R [J/mol-K], 1.4 = gamma
v0 = M0*a0 # [m/s]

mf_air_init = rho0 * A_fan * v0





''' FORMULAE

Stations:
    0 - ambient
    13 - bypassed air
    2 - entrance of the fan
    21 - entrance of LPC
    25 - entrance of HPC
    3 - entrance of cc
    4 - exit of cc
    45 - exit of HPT
    5 - exit of LPT
    18 - exit of bypassed air
    8 - exit of nozzle from core
'''

i = np.where( phases == 'cruise' )[0][0]

# Total temperature and pressure at inlet
T00 = T0[i] * ( 1 + (k_air-1)/2 * M0[i]**2 )
p00 = p0[i] * ( 1 + (k_air-1)/2 * M0[i]**2 ) ** (k_air / (k_air-1) )
print('\nInlet:\nT0 = ', T0[i], '[K]; p0 = ', p0[i], ',[Pa]; T00 = ', T00, '[K]; p00 = ', p00, '[Pa]')

# Entrance of the fan
T02 = T00
p02 = p0[i] * (1 + eta_inlet * (k_air-1)/2 * M0[i]**2 ) ** (k_air / (k_air-1) )
print('\nEntrance of fan:\nT02 = ', T02, '[K]; p02 = ', p02, '[Pa]')

# Exit of the fan - Entrance of LPC
T021 = T02 + ( T02/eta_fan ) * ( PR_fan ** ( (k_air-1)/k_air ) - 1 )
p021 = p02 * PR_fan
print('\nEntrance of LPC:\nT021 = ', T021, '[K]; p021 = ', p021, '[Pa]')

# Hot and cold mass flow of air
mf_hot = mf_air_init[i] / (BR+1)
mf_cold = mf_air_init[i] * BR / (BR+1)
print('\nMass flow of air:\nTotal = ', mf_air_init[i], '[kg/s]; Core = ', mf_hot, '[kg/s]; Bypassed = ', mf_cold,'[kg/s]')

# Exit of LPC - Entrance of HPC
T025 = T021 + ( T021/eta_LPC ) * ( PR_LPC ** ( (k_air-1)/k_air ) - 1 )
p025 = p021 * PR_LPC
print('\nEntrance of HPC:\nT025 = ', T025, '[K]; p025 = ', p025, '[Pa]')

# Exit of HPC - Entrance of cc
T03 = T025 + ( T025/eta_HPC ) * ( PR_HPC ** ( (k_air-1)/k_air ) - 1 )
p03 = p025 * PR_HPC
print('\nEntrance of CC:\nT03 = ', T03, '[K]; p03 = ', p03, '[Pa]')

# Exit of cc - Entrance of HPT
mf_fuel = (mf_hot * ratio_air_cc * cp_gas * (T04-T03)) / (LHV_f[i]*10**6 * eta_cc)
mf_airfuel = mf_hot + mf_fuel # at the end of the cc
print('\nMass flow CC:\nFuel = ', mf_fuel, '[kg/s]; air CC = ', mf_hot*ratio_air_cc, '[kg/s]; Total end of CC = ', mf_airfuel,'[kg/s]')

# T04 = 1500 [K], is given
p04 = p03 * PR_cc
print('\nEntrance of HPT:\nT04 = ', T04, '[K]; p04 = ', p04, '[Pa]')

# Power to drive fan, LPC, HPC, HPT, LPT [W]
W_fan = mf_air_init[i] * cp_air * (T021-T00)
W_LPC = mf_hot * cp_air * (T025-T021)
W_HPC = mf_hot * cp_air * (T03-T025)
W_HPT = W_HPC / eta_mech
W_LPT = (W_fan + W_LPC) / eta_mech

# Exit of HPT - Entrance of LPT
T045 = T04 - W_HPT / ( mf_airfuel * cp_gas )
p045 = p04 * ( 1 - ( 1 - T045/T04 ) / eta_HPT ) ** ( k_gas / (k_gas-1) )
print('\nEntrance of LPT:\nT045 = ', T045, '[K]; p045 = ', p045, '[Pa]')

# Exit of LPT - Entrance of nozzle
T05 = T045 - W_LPT / (mf_airfuel * cp_gas)
p05 = p045 * ( 1 - ( 1 - T05/T045 ) / eta_LPT ) ** ( k_gas / (k_gas-1) )
print('\nEntrance of nozzle:\nT05 = ', T05, '[K]; p05 = ', p05, '[Pa]')

# Is the nozzle chocked?
PR_cr_nozzle = 1 / ( ( 1 - (k_gas-1)/(eta_nozzle*(k_gas+1)) ) ** (k_gas / (k_gas-1)) )
if p05/p0[i] > PR_cr_nozzle:
    print('The nozzle is chocked')
else:
    print('The nozzle is NOT chocked')
    
# Is the fan chocked?
PR_cr_fan = 1 / ( ( 1 - (k_air-1)/(eta_nozzle*(k_air+1)) ) ** (k_air / (k_air-1)) )
if p021/p0[i] > PR_cr_fan:
    print('The fan is chocked')
else:
    print('The fan is NOT chocked')

# Exit of the nozzle
TR_cr_nozzle = ( k_gas + 1 ) / 2
T8 = T05 / TR_cr_nozzle
p8 = p05 / PR_cr_nozzle
v8 = m.sqrt( k_gas * R * T8 )

# Exit of bypassed air
TR_cr_bypassed = ( k_air + 1 ) / 2
T18 = T021 / TR_cr_bypassed
p18 = p021 / PR_fan
v18 = m.sqrt( k_air * R * T18 )

# Find thrust of fan, core, and total
A18 = (mf_cold * R * T18) / (p18 * v18)
T_fan = mf_cold * (v18 - v0[i]) + A18 * (p18 - p0[i]) # [N]
A8 = (mf_airfuel * R * T8) / (p8 * v8)
T_core = mf_airfuel * (v8 - v0[i]) + A8 * (p8 - p0[i]) # [N]
T_total = T_fan + T_core # [N]
print('\nProvided Trhust:\nFan = ', T_fan, '[N]; Core = ', T_core, '[N]; Total = ', T_total, '[N]')
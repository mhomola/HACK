from Subsystem_design.Environment.Climate_impact_assessment import Climate_assess
from Subsystem_design.Engine.CoolEngine import Engine_Cool
from Subsystem_design.Engine.EngineCycle import Engine_Cycle
from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.EnergySplit import Energy_Split
from Subsystem_design.Engine.thrust_req import thrust_req
import numpy as np
import matlab.engine


'''Run analysis for both neo and HACK'''
#Play with the Energy split between each phase to achieve requirements and with the mixture between air and fuel.

#1. Run cycle analysis for the selected split between hydrogen and kerosene. Should design for harmonic range of neo
# which is 3200 km
    #Checks:
        #1. Enoungh Thrust is provided at each phase. If not change (tbd)
        #2. Range requirement is achieved. If not change (tbd).
        #3. Emissions requirement is achieved. [note: for this must run A320neo vs A320HACK]
            # 3 points on LTO;
            # 3 points cruise;

#2. Run code to cool the engine and check if error is smaller than 5% some iteration needs to be done.
    #Checks:
        #1. Temperature is reduced enough?? (Ask Sara)


#3. Iterate until it converges: Compare old equivalence ratio with new equivalence ratio. 1% difference.

climate = Climate_assess(t = 2135)
eng = matlab.engine.start_matlab()
cycle = Engine_Cycle()
const = Constants()
cool = Engine_Cool()

aircraft = ['neo', 'hack']
phases = ['taxi_out', 'take_off', 'climb', 'cruise', 'approach', 'taxi_in']

#First run it for neo, store the ATR value to later compare with HACK
'-------------------------NEO----------------------------------'
print("\n= = = = Analysis for A320", aircraft[0], "= = = =")

T_neo = []
TSFC_neo = []

for b in phases:
    print("\n", b)
    '''Run cycle analysis, to get TSFC, T_tot'''
    cycle.cycle_analysis(aircraft[0],b)
    print('Entrance of CC: T03 = ', round(cycle.T03, 3), '[K]; p03 = ', round(cycle.p03, 3), '[Pa]; OPR = ', round(cycle.OPR, 3))
    print('Mass flow CC: Fuel = ', round(cycle.mf_fuel, 3), '[kg/s]; air CC = ', round(cycle.mf_hot, 3),
          '[kg/s]; Total end of CC = ', round(cycle.mf_airfuel, 3), '[kg/s]')
    print('LHV fuel = ', round(cycle.LHV_f, 3), 'm air to cool / m air core', round(cycle.mr_SZair_simpl, 4),
          round(cycle.mr_SZair_simpl1, 4))
    print('TPZ = ', round(cycle.TPZ, 3))

    # Get TPZ from Ivan's code. Inputs are ( gas, P, T, phi )
    (TPZ, MF, MF_names) = eng.reactor1('neo', float(cycle.p03), float(cycle.T03), float(cycle.equivalence_ratio))

    cool.SZ_air(cycle.mf_hot, cycle.mf_h2, cycle.mf_ker, cycle.T03, TPZ, cycle.T04, cycle.mr_SZair_simpl1,
                cycle.LHV_f)





'''--------------------------------HACK--------------------------------------------'''
print("\n= = = = Analysis for A320", aircraft[1], "= = = =")

Range_requirement = False
Thrust_requirement = False
Emissions_requirement = False
Flamability_requirement = False


for b in phases:
    print("\n", b)
    '''Run cycle analysis, to get TSFC, T_tot'''
    cycle.cycle_analysis(aircraft[1],b)
    print('Entrance of CC: T03 = ', round(cycle.T03, 3), '[K]; p03 = ', round(cycle.p03, 3), '[Pa]; OPR = ', round(cycle.OPR, 3))
    print('Mass flow CC: Fuel = ', round(cycle.mf_fuel, 3), '[kg/s]; air CC = ', round(cycle.mf_hot, 3),
          '[kg/s]; Total end of CC = ', round(cycle.mf_airfuel, 3), '[kg/s]')
    print('LHV fuel = ', round(cycle.LHV_f, 3), 'm air to cool / m air core', round(cycle.mr_SZair_simpl, 4),
          round(cycle.mr_SZair_simpl1, 4))
    print('TPZ = ', round(cycle.TPZ, 3))
    #Compare Thrust from cycle analysis to Thrust required
    if Thrust_cycle_analysis< Thrust_required:
        print('Thrust available for '+b+ ' is too little, change...')
    else:
        print('Thrust available for '+b+ 'is enough: Thrust_cycle_analysis.)
        Thrust_requirement = True

if range_analysis< 3200:

    print('Harmonic Range is:',range_analysis)
    print('Requirement of range not met.')
else:
    Range_requirement = True






from Subsystem_design.Environment.Climate_impact_assessment import Climate_assess
from Subsystem_design.Engine.CoolEngine import Engine_Cool
from Subsystem_design.Engine.EngineCycle import Engine_Cycle
from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.EnergySplit import Energy_Split
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

#2. Run code to cool the engine and check if error is smaller than 5% some iteration needs to be done.
    #Checks:
        #1. Temperature is reduced enough?? (Ask Sara)

#3. Iterate until it converges

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
for b in phases:
    print("\n", b)
    cycle.cycle_analysis(aircraft[0],b)
    # Get TPZ from Ivan's code. Inputs are ( gas, P, T, phi )
    (TPZ, MF, MF_names) = eng.reactor1('neo', float(cycle.p03), float(cycle.T03), float(cycle.equivalence_ratio))
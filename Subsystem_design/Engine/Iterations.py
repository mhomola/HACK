from Subsystem_design.Environment.Climate_impact_assessment import Climate_assess
from Subsystem_design.Engine.CoolEngine import Engine_Cool, get_TPZ
from Subsystem_design.Engine.EngineCycle import Engine_Cycle
from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.DataEngine import DataFrame
from Subsystem_design.Engine.EnergySplit import Energy_Split
from Subsystem_design.Engine.Thrust_Required import thrust_req
from Subsystem_design.aerodynamic_subsys import cd0clean, wingar
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
ec = Engine_Cycle()
const = Constants()
cool = Engine_Cool()
T_required = thrust_req(cd0clean,wingar)

aircraft = ['neo', 'hack']
phases = np.array(['taxi_out', 'take_off', 'climb', 'cruise', 'approach', 'taxi_in'])
printing = False                                                            #Change to True for printing stuff

#First run it for neo, store the ATR value to later compare with HACK
'-------------------------NEO----------------------------------'
print("\n= = = = Analysis for A320", aircraft[0], "= = = =")

'CH4, CO, CO2, H2O, NO, NO2, H2'
Emissions_array_neo = np.zeros(7)

for b in phases:
    print("\n", b)
    '''Run cycle analysis, to get TSFC, T_tot'''
    ec.cycle_analysis(aircraft[0],b)

    if printing== True:
        print('\nInlet: T0 = ', round(ec.T0, 3), '[K]; p0 = ', round(ec.p0, 3), '[Pa]; v0 = ', round(ec.v0, 3), '[m/s]')
        print('T00 = ', round(ec.T00, 3), '[K]; p00 = ', round(ec.p00, 3), '[Pa]')
        print('Entrance of fan: T02 = ', round(ec.T02, 3), '[K]; p02 = ', round(ec.p02, 3), '[Pa]')
        print('Entrance of LPC: T021 = ', round(ec.T021, 3), '[K]; p021 = ', round(ec.p021, 3), '[Pa]')
        print('Mass flow of air: Total = ', round(ec.mf_air_init, 3), '[kg/s]; Core = ', round(ec.mf_hot, 3),
              '[kg/s]; Bypassed = ', round(ec.mf_cold, 3), '[kg/s]')
        print('Entrance of HPC: T025 = ', round(ec.T025, 3), '[K]; p025 = ', round(ec.p025, 3), '[Pa]')
        print('Entrance of CC: T03 = ', round(ec.T03, 3), '[K]; p03 = ', round(ec.p03, 3), '[Pa]; OPR = ',
              round(ec.OPR, 3))
        print('Mass flow CC: Fuel = ', round(ec.mf_fuel, 3), '[kg/s]; air CC = ', round(ec.mf_hot, 3),
              '[kg/s]; Total end of CC = ', round(ec.mf_airfuel, 3), '[kg/s]')
        print('LHV fuel = ', round(ec.LHV_f, 3), 'm air to cool / m air core', round(ec.mr_SZair_simpl, 4),
              round(ec.mr_SZair_simpl1, 4))
        print('Initial estimate TPZ = ', round(ec.TPZ, 3))
        print('Power: Fan = ', round(ec.W_fan, 3), '[W]; LPC = ', round(ec.W_LPC, 3), '[W]; HPC = ', round(ec.W_HPC, 3),
              '[W]')
        print('LPT = ', round(ec.W_LPT, 3), '[W]; HPT = ', round(ec.W_HPT, 3), '[W]')
        print('Entrance of HPT: T04 = ', round(ec.T04, 3), '[K]; p04 = ', round(ec.p04, 3), '[Pa]')
        print('Entrance of LPT: T045 = ', round(ec.T045, 3), '[K]; p045 = ', round(ec.p045, 3), '[Pa]')
        print('Entrance of nozzle: T05 = ', round(ec.T05, 3), '[K]; p05 = ', round(ec.p05, 3), '[Pa]')
        print('Exit of nozzle: T07 = ', round(ec.T07, 3), '[K]; p07 = ', round(ec.p07, 3), '[Pa]; PR_cr_noz = ',
              ec.PR_cr_noz_core)
        print('Exit of nozzle: T8 = ', round(ec.T8, 3), '[K]; p8 = ', round(ec.p8, 3), '[Pa]; v8 = ', round(ec.v8, 3),
              '[m/s]')
        print('Exit of fan: T016 = ', round(ec.T016, 3), '[K]; p016 = ', round(ec.p016, 3), '[Pa]; PR_cr_fan = ',
              ec.PR_cr_fan)
        print('Exit of fan: T18 = ', round(ec.T18, 3), '[K]; p18 = ', round(ec.p18, 3), '[Pa]; v18 = ',
              round(ec.v18, 3), '[m/s]')
        print('Provided Thrust: Fan = ', round(ec.T_fan, 3), '[N]; Core = ', round(ec.T_core, 3), '[N]; Total = ',
              round(ec.T_total, 3), '[N]')
        print('Thrust SFC = ', round(ec.TSFC, 5), '[g/kN/s]; Equivalence ratio = ', round(ec.equivalence_ratio, 4))

    '''Getting moles/second of H2 and kerosene'''
    n_h2,n_ker,n_O2,n_N2 = ec.n_h2,ec.n_ker,ec.n_O2,ec.n_N2
    print(ec.mf_h2,ec.mf_ker)
    print(n_h2,n_ker,n_O2,n_N2)

    cool.SZ_air(aircraft[0], b, ec.TPZ)
    eqr_old = cool.eqr
    print('Initial TPZ [K]:', ec.TPZ, ' Initial mr_cool', cool.mr_SZair, ' Initial eqr', cool.eqr)

    ''' LOOP FOR CONVERGENCE OF EQUIVALENCE RATIO '''
    eqr_old = cool.eqr.copy()
    err = 1

    while err > 0.02:  # error larger than 2%
        print(err)

        TPZ, Emissions = get_TPZ(aircraft[0], b, ec.p03, ec.T03, cool.eqr)
        cool.SZ_air(aircraft[0], b, TPZ)

        err = abs(cool.eqr - eqr_old) / cool.eqr
        eqr_old = cool.eqr.copy()

        print('Updated TPZ:', TPZ, ' Updated MR:', cool.mr_SZair, 'Updated eqr:', cool.eqr)

    '''ASSESSING THE CLIMATE IMPACT OF NEO'''
    time = T_required.durations[phases==b]
    print('Time',time)
    U_ker = climate.number_aircraft_kerosene()
    ATR = climate.ATR(h=ec.h, e_CO2=e_CO2, e_H2O=e_H2O, e_NOx=e_NOx, e_soot=e_soot, e_sulfate=e_SO4, U=U_ker, plot=True)
    print('The average temperature response, A_100, for the LTO of the HACK is:', ATR, '[K]')







'''--------------------------------HACK--------------------------------------------'''
print("\n= = = = Analysis for A320", aircraft[1], "= = = =")

T_HACK = []
TSFC_HACK = []
Range_requirement = False
Thrust_requirement = False
Emissions_requirement = False
Flamability_requirement = False

'CH4, CO, CO2, H2O, NO, NO2, H2'
Emissions_array_HACK = np.zeros(7)

while Range_requirement != True:
    for b in phases:
        print("\n", b)
        '''Run cycle analysis, to get TSFC, T_tot'''
        ec.cycle_analysis(aircraft[1],b)
        if printing== True:
            print('\nInlet: T0 = ', round(ec.T0, 3), '[K]; p0 = ', round(ec.p0, 3), '[Pa]; v0 = ', round(ec.v0, 3), '[m/s]')
            print('T00 = ', round(ec.T00, 3), '[K]; p00 = ', round(ec.p00, 3), '[Pa]')
            print('Entrance of fan: T02 = ', round(ec.T02, 3), '[K]; p02 = ', round(ec.p02, 3), '[Pa]')
            print('Entrance of LPC: T021 = ', round(ec.T021, 3), '[K]; p021 = ', round(ec.p021, 3), '[Pa]')
            print('Mass flow of air: Total = ', round(ec.mf_air_init, 3), '[kg/s]; Core = ', round(ec.mf_hot, 3),
                  '[kg/s]; Bypassed = ', round(ec.mf_cold, 3), '[kg/s]')
            print('Entrance of HPC: T025 = ', round(ec.T025, 3), '[K]; p025 = ', round(ec.p025, 3), '[Pa]')
            print('Entrance of CC: T03 = ', round(ec.T03, 3), '[K]; p03 = ', round(ec.p03, 3), '[Pa]; OPR = ',
                  round(ec.OPR, 3))
            print('Mass flow CC: Fuel = ', round(ec.mf_fuel, 3), '[kg/s]; air CC = ', round(ec.mf_hot, 3),
                  '[kg/s]; Total end of CC = ', round(ec.mf_airfuel, 3), '[kg/s]')
            print('LHV fuel = ', round(ec.LHV_f, 3), 'm air to cool / m air core', round(ec.mr_SZair_simpl, 4),
                  round(ec.mr_SZair_simpl1, 4))
            print('Initial estimate TPZ = ', round(ec.TPZ, 3))
            print('Power: Fan = ', round(ec.W_fan, 3), '[W]; LPC = ', round(ec.W_LPC, 3), '[W]; HPC = ', round(ec.W_HPC, 3),
                  '[W]')
            print('LPT = ', round(ec.W_LPT, 3), '[W]; HPT = ', round(ec.W_HPT, 3), '[W]')
            print('Entrance of HPT: T04 = ', round(ec.T04, 3), '[K]; p04 = ', round(ec.p04, 3), '[Pa]')
            print('Entrance of LPT: T045 = ', round(ec.T045, 3), '[K]; p045 = ', round(ec.p045, 3), '[Pa]')
            print('Entrance of nozzle: T05 = ', round(ec.T05, 3), '[K]; p05 = ', round(ec.p05, 3), '[Pa]')
            print('Exit of nozzle: T07 = ', round(ec.T07, 3), '[K]; p07 = ', round(ec.p07, 3), '[Pa]; PR_cr_noz = ',
                  ec.PR_cr_noz_core)
            print('Exit of nozzle: T8 = ', round(ec.T8, 3), '[K]; p8 = ', round(ec.p8, 3), '[Pa]; v8 = ', round(ec.v8, 3),
                  '[m/s]')
            print('Exit of fan: T016 = ', round(ec.T016, 3), '[K]; p016 = ', round(ec.p016, 3), '[Pa]; PR_cr_fan = ',
                  ec.PR_cr_fan)
            print('Exit of fan: T18 = ', round(ec.T18, 3), '[K]; p18 = ', round(ec.p18, 3), '[Pa]; v18 = ',
                  round(ec.v18, 3), '[m/s]')
            print('Provided Thrust: Fan = ', round(ec.T_fan, 3), '[N]; Core = ', round(ec.T_core, 3), '[N]; Total = ',
                  round(ec.T_total, 3), '[N]')
            print('Thrust SFC = ', round(ec.TSFC, 5), '[g/kN/s]; Equivalence ratio = ', round(ec.equivalence_ratio, 4))


        '''Get TPZ from Ivan's code. Inputs are ( gas, P, T, phi )'''
        eqr_old = ec.equivalence_ratio
        TPZ, Emissions = get_TPZ(aircraft[0], b, ec.p03, ec.T03, ec.equivalence_ratio)
        print('1st TPZ from Matlab:', TPZ)
        cool.SZ_air(aircraft[0], b, TPZ)
        print('1st MR from engine cycle:', ec.mr_SZair_simpl1, 'MR with this new TPZ:', cool.mr_SZair)

        ''' LOOP FOR CONVERGENCE OF EQUIVALENCE RATIO '''
        eqr_new = cool.eqr
        err = abs(eqr_new - eqr_old) / eqr_new
        eqr_old = eqr_new  # to start while loop
        print('Error', err)
        while err > 0.02:  # error larger than 2%
            print(err)
            TPZ, Emissions = get_TPZ(aircraft[0], b, ec.p03, ec.T03, cool.eqr)
            cool.SZ_air(aircraft[0], b, TPZ)
            err = abs(cool.eqr - eqr_old) / cool.eqr
            eqr_old = cool.eqr
            print('Updated TPZ:', TPZ, ' Updated MR:', cool.mr_SZair, 'Updated eqr:', cool.eqr)


        #Compare Thrust from cycle analysis to Thrust required
        if ec.T_total< Thrust_required:
            print('Thrust available for '+b+ ' is too little, change...')
        else:
            print('Thrust available for '+b+ 'is enough:',ec.T_total)
            Thrust_requirement = True

        T_HACK.append(ec.T_total)
        TSFC_HACK.append(ec.TSFC)

    if range_analysis< 3200:

        print('Harmonic Range is:',range_analysis)
        print('Requirement of range not met.')
    else:
        Range_requirement = True






from Subsystem_design.Environment.Climate_impact_assessment import Climate_assess
from Subsystem_design.Engine.CoolEngine import Engine_Cool, get_TPZ
from Subsystem_design.Engine.EngineCycle import Engine_Cycle
from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.DataEngine import DataFrame
from Subsystem_design.Engine.EnergySplit import Energy_Split
from Subsystem_design.Engine.Thrust_Required import thrust_req
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
from Subsystem_design.aerodynamic_subsys import cd0clean, wingar
from Subsystem_design.Performance.Weight_estimation_NEW import performance
from Subsystem_design.Performance.Weight_estimation_NEW import Compute_weight
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
T = thrust_req(cd0clean,wingar)
T.drag()
ae = AerodynamicCharacteristics()
ae.aero_functions(AoA_cruise=2)

aircraft = ['neo', 'hack']
phases = np.array(['taxi_out', 'take_off', 'climb', 'cruise', 'approach', 'taxi_in'])
#phases = np.array(['taxi_out', 'take_off1', 'take_off2', 'climb1', 'climb2', 'cruise1', 'cruise2', 'approach1',
          #'approach2', 'taxi_in', 'idle'])
#Time_phases = 60* np.array([8.,0.5,1.,1.,18,217,16.5,0.5,1.,8.])
#Time_CA_analysis = 60*np.array([3.75,8.,8.5,9.5,10.5,28.5,245.5,262,262.5,263.5,264.5])
Time_phases = 60* np.array([7.5,1.,2.,217,1.5,8.])
printing = False                #Change to True for printing stuff

#Emission indexes of soot and UHC
'taxi_out,take_off,climb,cruise,approach,taxi_in'
EI_soot = np.array([2.0 * 10**-4,2.0 * 10**-4,2.0 * 10**-4,2.0 * 10**-4,2.0 * 10**-4,2.0 * 10**-4])
EI_UHC = np.array([0.28*10**-3,0.02*10**-3,0.02*10**-3,0.04*10**-3,0.04*10**-3,0.28*10**-3])
#First run it for neo, store the ATR value to later compare with HACK
'-------------------------NEO----------------------------------'
print("\n= = = = Analysis for A320", aircraft[0], "= = = =")

'CH4, CO, CO2, H2O, NO, NO2, H2, soot, UHC'
Emissions_array_neo = np.zeros((len(phases),9))



for b in phases:
    print("\n", b)
    '''Run cycle analysis, to get TSFC, T_tot'''
    ec.cycle_analysis(aircraft[0],b,flag=False, alph=0)

    if printing== True:
        print('\nInlet: T0 = ', round(ec.T0, 3), '[K]; p0 = ', round(ec.p0, 3), '[Pa]; v0 = ', round(ec.v0, 3),
              '[m/s]')
        print('T00 = ', round(ec.T00, 3), '[K]; p00 = ', round(ec.p00, 3), '[Pa]')
        print('Entrance of fan: T02 = ', round(ec.T02, 3), '[K]; p02 = ', round(ec.p02,3), '[Pa]')
        print('Entrance of LPC: T021 = ', round(ec.T021, 3), '[K]; p021 = ', round(ec.p021,3), '[Pa]')
        print('Mass flow of air: Total = ', round(ec.mf_air_init, 3), '[kg/s]; Core = ', round(ec.mf_hot, 3),
        '[kg/s]; Bypassed = ', round(ec.mf_cold, 3), '[kg/s]')
        print('Entrance of HPC: T025 = ', round(ec.T025, 3), '[K]; p025 = ', round(ec.p025, 3), '[Pa]')
        print('Entrance of CC: T03 = ', round(ec.T03, 3), '[K]; p03 = ', round(ec.p03, 3), '[Pa]; OPR = ',
              round(ec.OPR, 3))
        print('Mass flow CC: Fuel = ', round(ec.mf_fuel, 3), '[kg/s]; mf_h2 = ', round(ec.mf_h2,3),
        '[kg/s]; mf_ker = ', round(ec.mf_ker, 3), '[kg/s]')
        print('Moles reacting per second of kerosene:', ec.n_ker, 'of H2:', ec.n_h2, 'of O2:', ec.n_O2, 'of N2:',
              ec.n_N2)
        print('LHV fuel = ', round(ec.LHV_f, 3))
        print('Initial estimate TPZ = ', round(ec.TPZ, 3))
        print('Power: Fan = ', round(ec.W_fan, 3), '[W]; LPC = ', round(ec.W_LPC, 3), '[W]; HPC = ',
              round(ec.W_HPC, 3), '[W]')
        print('LPT = ', round(ec.W_LPT, 3), '[W]; HPT = ', round(ec.W_HPT, 3), '[W]')
        print('Entrance of HPT: T04 = ', round(ec.T04, 3), '[K]; p04 = ', round(ec.p04, 3), '[Pa]')
        print('Entrance of LPT: T045 = ', round(ec.T045, 3), '[K]; p045 = ', round(ec.p045, 3), '[Pa]')
        print('Entrance of nozzle: T05 = ', round(ec.T05, 3), '[K]; p05 = ', round(ec.p05, 3), '[Pa]')
        print('Exit of nozzle: T07 = ', round(ec.T07, 3), '[K]; p07 = ', round(ec.p07, 3), '[Pa]; PR_cr_noz = ',
              round(ec.PR_cr_noz_core, 3))
        print('Exit of nozzle: T8 = ', round(ec.T8, 3), '[K]; p8 = ', round(ec.p8, 3), '[Pa]; v8 = ',
              round(ec.v8, 3), '[m/s]')
        print('Exit of fan: T016 = ', round(ec.T016, 3), '[K]; p016 = ', round(ec.p016, 3), '[Pa]; PR_cr_fan = ',
              round(ec.PR_cr_fan, 3))
        print('Exit of fan: T18 = ', round(ec.T18, 3), '[K]; p18 = ', round(ec.p18, 3), '[Pa]; v18 = ',
              round(ec.v18,3), '[m/s]')
        print('Provided Thrust: Fan = ', round(ec.T_fan, 3), '[N]; Core = ', round(ec.T_core, 3), '[N]; Total = ',
              round(ec.T_total, 3), '[N]')
        print('Thrust SFC = ', round(ec.TSFC_m, 5), '[g/kN/s];', round(ec.TSFC_e, 5), '[MJ/kN/s]')
        print('Eqr PZ:', round(ec.eqr_PZ, 3), 'Eqr Overall:', round(ec.eqr_overall, 3))

    '''Getting moles/second of H2 and kerosene'''
    n_h2,n_ker,n_O2,n_N2 = ec.n_h2,ec.n_ker,ec.n_O2,ec.n_N2                     # Number of moles/sec for a
                                                                                # stoichiometric reaction
    TPZ, Emissions = get_TPZ(aircraft[0], b, ec.p03, ec.T03, ec.eqr_PZ, n_h2, n_ker, n_O2, n_N2)
    print('TPZ is:', TPZ)

    '''ASSESSING THE EMISSIONS OF NEO'''
    #Getting time of phase
    time = Time_phases[phases==b]

    #Getting total emissions in that phase
    mf_reactants = ec.mf_fuel + ec.mf_hot
    Emissions = np.array(Emissions)* mf_reactants * time
    E_soot = ec.mf_ker * EI_soot[phases==b] * time
    E_UHC = ec.mf_ker * EI_UHC[phases==b] * time

    Emissions = np.append(Emissions,[E_UHC,E_soot])
    Emissions_array_neo[phases == b] = Emissions


    print(np.array(['CH4', 'CO', 'CO2', 'H2O', 'NO', 'NO2', 'H2','soot','UHC']))
    print('Emissions array:',Emissions_array_neo)



'''ASSESSING THE CLIMATE IMPACT OF NEO'''
#For the LTO cycle [taxi-out,take-off,clim-out(until 900m),descend(after 900m),landing,taxi-in]
eCO2 = np.sum(Emissions_array_neo[:,2])-Emissions_array_neo[3,2]
eCO = np.sum(Emissions_array_neo[:,1])-Emissions_array_neo[3,1]
eH2O = np.sum(Emissions_array_neo[:,3])-Emissions_array_neo[3,3]
eNOx = np.sum(Emissions_array_neo[:,4])-Emissions_array_neo[3,4] + np.sum(Emissions_array_neo[:,5])-Emissions_array_neo[3,5]
esoot = np.sum(Emissions_array_neo[:,7])-Emissions_array_neo[3,7]
eHUC = np.sum(Emissions_array_neo[:,8])-Emissions_array_neo[3,8]
eSO4 = 0
TE_LTO_neo= np.array([eCO,eCO2,eH2O,eNOx,esoot,eHUC])
print('[eCO,eCO2,eH2O,eNOx,esoot,eHUC]')
print(TE_LTO_neo)

U_ker = climate.number_aircraft_H2()
h = 450
ATR_LTO_neo = climate.ATR(h=h, e_CO2=eCO2, e_H2O=eH2O, e_NOx=eNOx, e_soot=esoot, e_sulfate=eSO4, U=U_ker, plot=True)
print('The average temperature response, A_100, for the LTO of the neo is:', ATR_LTO_neo, '[K]')

#For the cruise phase
eCO2 = Emissions_array_neo[3,2]
eCO = Emissions_array_neo[3,1]
eH2O = Emissions_array_neo[3,3]
eNOx = Emissions_array_neo[3,4] + Emissions_array_neo[3,5]
esoot = Emissions_array_neo[3,7]
eHUC = Emissions_array_neo[3,8]
eSO4 = 0
TE_cruise_neo = np.array([eCO,eCO2,eH2O,eNOx,esoot,eHUC])
print(TE_cruise_neo)
h = 11600
ATR_cruise_neo = climate.ATR(h=h, e_CO2=eCO2, e_H2O=eH2O, e_NOx=eNOx, e_soot=esoot, e_sulfate=eSO4, U=U_ker, plot=True)
print('The average temperature response, A_100, for the cruise of the neo is:', ATR_cruise_neo, '[K]')

'''Storing Neo emissions'''

file = open('Emissions_neo', 'w')
file.write(str(Emissions_array_neo) + 'ATR during LTO' + str(ATR_LTO_neo)+ 'ATR during cruise' + str(ATR_cruise_neo))
file.close()


'''--------------------------------HACK--------------------------------------------'''
print("\n= = = = Analysis for A320", aircraft[1], "= = = =")

T_HACK = []
TSFC_HACK = []
data_list = []

#Thrust required at each phase of the flight
T_req = np.array( [29570/2, 156900/2, 93800/2, 41066/2,67158/2, 30125.4/2])

Range_requirement = False
Thrust_requirement = False
Emissions_requirement_LTO = False
Emissions_requirement_cruise = False
ATR_requirement = False
Flamability_requirement = False

'CH4, CO, CO2, H2O, NO, NO2, H2,soot,UHC'
Emissions_array_HACK = np.zeros((len(phases),9))



while Emissions_requirement_LTO != True and Emissions_requirement_cruise!= True and ATR_requirement!= True:
    while Range_requirement != True:
        for b in phases:
            print("\n", b)
            '''Run cycle analysis, to get TSFC, T_tot'''
            ec.cycle_analysis(aircraft[1],b,flag=False, alph=1)
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
            n_h2, n_ker, n_O2, n_N2 = ec.n_h2, ec.n_ker, ec.n_O2, ec.n_N2                   # Number of moles/sec for a
                                                                                            # stoichiometric reaction
            TPZ, Emissions = get_TPZ(aircraft[1], b, ec.p03, ec.T03, ec.eqr_PZ, n_h2, n_ker, n_O2, n_N2)
            print('TPZ is:',TPZ)

            '''ASSESSING THE EMISSIONS OF NEO'''
            # Getting time of phase
            time = Time_phases[phases == b]

            # Getting total emissions in that phase
            mf_reactants = ec.mf_fuel + ec.mf_hot
            Emissions = np.array(Emissions)* mf_reactants * time
            E_soot = ec.mf_ker * EI_soot[phases == b] * time
            E_UHC = ec.mf_ker * EI_UHC[phases == b] * time

            Emissions = np.append(Emissions, [E_UHC, E_soot])
            Emissions_array_HACK[phases == b] = Emissions

            print(np.array(['CH4', 'CO', 'CO2', 'H2O', 'NO', 'NO2', 'H2']))
            print('Emissions array:', Emissions_array_HACK)

            #Compare Thrust from cycle analysis to Thrust required
            if ec.T_total< T_req[phases==b]:
                print('Thrust available for '+b+ ' is too little, change energy split')
            else:
                print('Thrust available for '+b+ 'is enough:',ec.T_total)
                Thrust_requirement = True

            T_HACK.append(ec.T_total)
            TSFC_HACK.append(ec.TSFC_m)

        '''ASSESSING RANGE OF HACK'''
        # m_5 = T_required.m[T_required.t_array<=T_required.cruise_end_time][-1]
        # m_4 = T_required.m[T_required.t_array<=T_required.cruise_start_time][-1]
        perf = performance()
        perf.mission_profile(T.durations,aircraft[1])
        m_5 = perf.fs_h2_cruise[-1] + perf.fs_k_cruise[-1] + Compute_weight().weight_break_down_HACK().MZFW_HACK
        m_4 = perf.fs_h2_climb[-1] + perf.fs_k_climb[-1] + Compute_weight().weight_break_down_HACK().MZFW_HACK
        Range_HACK = performance().Range(ae.L_D_ratio_HACK,m_5/m_4,TSFC_HACK[phases=='cruise']*10**-6)

        if Range_HACK< 3200*10**3:

            print('Harmonic Range is:',Range_HACK)
            print('Requirement of range not met. Change energy split')
        else:
            print('Range requirement is met')
            Range_requirement = True

    '''ASSESSING THE CLIMATE IMPACT OF HACK'''
    #For the LTO cycle [taxi-out,take-off,clim-out(until 900m),descend(after 900m),landing,taxi-in]
    eCO2 = np.sum(Emissions_array_HACK[:,2])-Emissions_array_HACK[3,2]
    eCO = np.sum(Emissions_array_HACK[:,1])-Emissions_array_HACK[3,1]
    eH2O = np.sum(Emissions_array_HACK[:,3])-Emissions_array_HACK[3,3]
    eNOx = np.sum(Emissions_array_HACK[:,4])-Emissions_array_HACK[3,4] + np.sum(Emissions_array_HACK[:,5])-Emissions_array_HACK[3,5]
    esoot = np.sum(Emissions_array_HACK[:, 7]) - Emissions_array_HACK[3, 7]
    eHUC = np.sum(Emissions_array_HACK[:, 8]) - Emissions_array_HACK[3, 8]
    eSO4 = 0
    TE_LTO_HACK = np.array([eCO, eCO2, eH2O, eNOx, esoot, eHUC])
    print('[eCO,eCO2,eH2O,eNOx,esoot,eHUC]')
    print(TE_LTO_HACK)
    U_H2 = climate.number_aircraft_H2()
    h = 450
    ATR_LTO_HACK = climate.ATR(h=h, e_CO2=eCO2, e_H2O=eH2O, e_NOx=eNOx, e_soot=esoot, e_sulfate=eSO4, U=U_H2, plot=True)
    print('The average temperature response, A_100, for the LTO of the HACK is:', ATR_LTO_HACK, '[K]')

    #For the cruise phase
    eCO2 = Emissions_array_HACK[3,2]
    eCO = Emissions_array_HACK[3,1]
    eH2O = Emissions_array_HACK[3,3]
    eNOx = Emissions_array_HACK[3,4] + Emissions_array_HACK[3,5]
    esoot = Emissions_array_HACK[3, 7]
    eHUC = Emissions_array_HACK[3, 8]
    eSO4 = 0
    TE_cruise_HACK = np.array([eCO, eCO2, eH2O, eNOx, esoot, eHUC])
    print(TE_cruise_HACK)
    h = 11600
    ATR_cruise_HACK = climate.ATR(h=h, e_CO2=eCO2, e_H2O=eH2O, e_NOx=eNOx, e_soot=esoot, e_sulfate=eSO4, U=U_H2, plot=True)
    print('The average temperature response, A_100, for the cruise of the HACK is:', ATR_cruise_HACK, '[K]')

    '''Storing hack emissions'''

    file = open('Emissions_HACK', 'w')
    file.write(str(Emissions_array_HACK) + 'ATR during LTO' + str(ATR_LTO_HACK) + 'ATR during cruise' + str(ATR_cruise_HACK))
    file.close()

    'Checking cruise emissions'
    if (TE_cruise_HACK[0]<0.3 *TE_cruise_neo[0]) and  (TE_cruise_HACK[3]<0.3 *TE_cruise_neo[3]) and (TE_cruise_HACK[4]<0.3 *TE_cruise_neo[4]) and (TE_cruise_HACK[5]<0.3 *TE_cruise_neo[5]):
        print('Emissions requirement met during cruise')
        Emissions_requirement_cruise = True
    else:
        print('Emissions requirement during cruise is not met. Change Energy split')

    if (TE_LTO_HACK[0]<0.5 *TE_LTO_neo[0]) and  (TE_LTO_HACK[3]<0.5 *TE_LTO_neo[3]) and (TE_LTO_HACK[4]<0.5 *TE_LTO_neo[4]) and (TE_LTO_HACK[5]<0.5 *TE_LTO_neo[5]):
        print('Emissions requirement met during cruise')
        Emissions_requirement_LTO = True
    else:
        print('Emissions requirement during LTO is not met. Change Energy split')

    if (ATR_cruise_HACK<0.3*ATR_cruise_neo) and (ATR_LTO_HACK<0.3*ATR_LTO_neo):
        ATR_requirement=True
        print('Climate_change requirement is met.')
    else:
        print('Climate_change requirement is not met. Change Energy split')


#--------------OLD CODE--------------------------------NEO-------------------------------
# #cool.SZ_air(aircraft[0], b, ec.TPZ)
#     #eqr_old = cool.eqr
#     #print('Initial TPZ [K]:', ec.TPZ, ' Initial mr_cool', cool.mr_SZair, ' Initial eqr', cool.eqr)
#
#     ''' LOOP FOR CONVERGENCE OF EQUIVALENCE RATIO '''
#
#     #eqr_old = 0.5350926260993377#0.7  # initial value
#     TPZ, MF = get_TPZ(aircraft[0], b, ec.p03, ec.T03, 0.7, n_h2, n_ker, n_O2, n_N2)
#     cool.SZ_air(aircraft[0], b, TPZ)
#     print('Initial TPZ [K]:', round(TPZ, 3), 'Initial mr_cool:', round(cool.mr_SZair, 3))
#     print(' Initial eqr:', 0.7, ' Updated eqr:', round(cool.eqr, 3))
#     err = 1
#
#     while err > 0.02:  # error larger than 2%
#
#         TPZ, Emissions = get_TPZ(aircraft[0], b, ec.p03, ec.T03, cool.eqr,n_h2,n_ker,n_O2,n_N2)
#
#         cool.SZ_air(aircraft[0], b, TPZ)
#         err = abs(cool.eqr - eqr_old) / cool.eqr
#         eqr_old = cool.eqr.copy()
#
#         print('Error at each iteration:', err * 100, '[%]')
#         print('Updated TPZ:', TPZ, ' Updated MR:', cool.mr_SZair, 'Updated eqr:', cool.eqr)
#
#
#     data_list.append([1 - cool.mr_SZair])
#
#     print('\nFINAL\nMass ratio of air injected on DZ:', round(cool.mr_SZair, 3))
#     print('TPZ = ', round(TPZ, 3))
#
#     np.savetxt('mr_cc_neo.dat', np.array(data_list))

#--------------OLD CODE--------------------------------HACK-------------------------------
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

Time_phases = 60* np.array([7.5,1.,2.,217,1.5,8.])
printing = False               #Change to True for printing stuff

#Emission indexes of soot and UHC
'taxi_out,take_off,climb,cruise,approach,taxi_in'
EI_soot = np.array([4.0* 10**-5,4.0* 10**-5,4.0* 10**-5,4.0* 10**-5,4.0* 10**-5,4.0* 10**-5])
EI_UHC = np.array([0.28*10**-3,0.02*10**-3,0.02*10**-3,0.04*10**-3,0.04*10**-3,0.28*10**-3])
#First run it for neo, store the ATR value to later compare with HACK

'-------------------------NEO----------------------------------'
# print("\n= = = = Analysis for A320", aircraft[0], "= = = =")
#
# 'CH4, CO, CO2, H2O, NO, NO2, H2, soot, UHC'
# Emissions_array_neo = np.zeros((len(phases),9))
#
#
#
# for b in phases:
#     print("\n", b)
#     '''Run cycle analysis, to get TSFC, T_tot'''
#     ec.cycle_analysis(aircraft[0],b,flag=False, alph=0)
#
#     if printing== True:
#         print('\nInlet: T0 = ', round(ec.T0, 3), '[K]; p0 = ', round(ec.p0, 3), '[Pa]; v0 = ', round(ec.v0, 3),
#               '[m/s]')
#         print('T00 = ', round(ec.T00, 3), '[K]; p00 = ', round(ec.p00, 3), '[Pa]')
#         print('Entrance of fan: T02 = ', round(ec.T02, 3), '[K]; p02 = ', round(ec.p02,3), '[Pa]')
#         print('Entrance of LPC: T021 = ', round(ec.T021, 3), '[K]; p021 = ', round(ec.p021,3), '[Pa]')
#         print('Mass flow of air: Total = ', round(ec.mf_air_init, 3), '[kg/s]; Core = ', round(ec.mf_hot, 3),
#         '[kg/s]; Bypassed = ', round(ec.mf_cold, 3), '[kg/s]')
#         print('Entrance of HPC: T025 = ', round(ec.T025, 3), '[K]; p025 = ', round(ec.p025, 3), '[Pa]')
#         print('Entrance of CC: T03 = ', round(ec.T03, 3), '[K]; p03 = ', round(ec.p03, 3), '[Pa]; OPR = ',
#               round(ec.OPR, 3))
#         print('Mass flow CC: Fuel = ', round(ec.mf_fuel, 3), '[kg/s]; mf_h2 = ', round(ec.mf_h2,3),
#         '[kg/s]; mf_ker = ', round(ec.mf_ker, 3), '[kg/s]')
#         print('Moles reacting per second of kerosene:', ec.n_ker, 'of H2:', ec.n_h2, 'of O2:', ec.n_O2, 'of N2:',
#               ec.n_N2)
#         print('LHV fuel = ', round(ec.LHV_f, 3))
#         print('Initial estimate TPZ = ', round(ec.TPZ, 3))
#         print('Power: Fan = ', round(ec.W_fan, 3), '[W]; LPC = ', round(ec.W_LPC, 3), '[W]; HPC = ',
#               round(ec.W_HPC, 3), '[W]')
#         print('LPT = ', round(ec.W_LPT, 3), '[W]; HPT = ', round(ec.W_HPT, 3), '[W]')
#         print('Entrance of HPT: T04 = ', round(ec.T04, 3), '[K]; p04 = ', round(ec.p04, 3), '[Pa]')
#         print('Entrance of LPT: T045 = ', round(ec.T045, 3), '[K]; p045 = ', round(ec.p045, 3), '[Pa]')
#         print('Entrance of nozzle: T05 = ', round(ec.T05, 3), '[K]; p05 = ', round(ec.p05, 3), '[Pa]')
#         print('Exit of nozzle: T07 = ', round(ec.T07, 3), '[K]; p07 = ', round(ec.p07, 3), '[Pa]; PR_cr_noz = ',
#               round(ec.PR_cr_noz_core, 3))
#         print('Exit of nozzle: T8 = ', round(ec.T8, 3), '[K]; p8 = ', round(ec.p8, 3), '[Pa]; v8 = ',
#               round(ec.v8, 3), '[m/s]')
#         print('Exit of fan: T016 = ', round(ec.T016, 3), '[K]; p016 = ', round(ec.p016, 3), '[Pa]; PR_cr_fan = ',
#               round(ec.PR_cr_fan, 3))
#         print('Exit of fan: T18 = ', round(ec.T18, 3), '[K]; p18 = ', round(ec.p18, 3), '[Pa]; v18 = ',
#               round(ec.v18,3), '[m/s]')
#         print('Provided Thrust: Fan = ', round(ec.T_fan, 3), '[N]; Core = ', round(ec.T_core, 3), '[N]; Total = ',
#               round(ec.T_total, 3), '[N]')
#         print('Thrust SFC = ', round(ec.TSFC_m, 5), '[g/kN/s];', round(ec.TSFC_e, 5), '[MJ/kN/s]')
#         print('Eqr PZ:', round(ec.eqr_PZ, 3), 'Eqr Overall:', round(ec.eqr_overall, 3))
#
#     '''Getting moles/second of H2 and kerosene'''
#     n_h2,n_ker,n_O2,n_N2 = ec.n_h2,ec.n_ker,ec.n_O2,ec.n_N2                     # Number of moles/sec for a
#                                                                                 # stoichiometric reaction
#     TPZ, Emissions = get_TPZ(aircraft[0], b, ec.p03, ec.T03, ec.eqr_PZ, n_h2, n_ker, n_O2, n_N2)
#     print('TPZ is:', TPZ)
#
#     '''ASSESSING THE EMISSIONS OF NEO'''
#     #Getting time of phase
#     time = Time_phases[phases==b]
#
#     #Getting total emissions in that phase
#     mf_reactants = ec.mf_fuel + ec.mf_hot
#     Emissions = np.array(Emissions)* mf_reactants * time
#     E_soot = ec.mf_ker * EI_soot[phases==b] * time
#     E_UHC = ec.mf_ker * EI_UHC[phases==b] * time
#
#     Emissions = np.append(Emissions,[E_UHC,E_soot])
#     Emissions_array_neo[phases == b] = Emissions
#
#
#     print(np.array(['CH4', 'CO', 'CO2', 'H2O', 'NO', 'NO2', 'H2','soot','UHC']))
#     print('Emissions array:',Emissions_array_neo)
#
#
# # Emissions_array_neo = np.array([[3.43528163e-13, 1.09462356e+00, 8.10842320e+02, 3.40139179e+02,4.54630880e+01, 3.22032032e-01, 1.62490711e-02, 2.76365322e-02,1.97403801e-02],
# #                                 [5.19236348e-15, 5.54736858e-01, 3.20939593e+02, 1.34590891e+02,1.86476927e+01, 1.25431250e-01, 8.04912409e-03, 7.82287029e-04,7.82287029e-03],
# #                                 [2.55190928e-13, 2.33653572e+01, 5.54140085e+02, 2.41622124e+02,2.55039721e+01, 4.19759572e-02, 2.89806817e-01, 1.45242691e-03,1.45242691e-02],
# #                                 [5.34918943e-10, 6.24238542e+01, 2.79690164e+04, 1.17315741e+04,1.37682797e+03, 5.49879546e+00, 9.34129264e-01, 1.36566020e-01,6.82830098e-01],
# #                                 [1.62824490e-12, 3.72257359e-01, 3.10388188e+02, 1.30246914e+02,1.56034080e+01, 9.25452711e-02, 5.72144294e-03, 1.51071509e-03,7.55357543e-03],
# #                                 [7.20756730e-13, 3.14348188e+00, 2.27802709e+03,9.55566892e+02,1.28671607e+02, 9.11869372e-01, 4.65405453e-02, 7.76490364e-02,5.54635975e-02]])
#
# '''ASSESSING THE CLIMATE IMPACT OF NEO'''
# #For the LTO cycle [taxi-out,take-off,clim-out(until 900m),descend(after 900m),landing,taxi-in]
# eCO2 = np.sum(Emissions_array_neo[:,2])-Emissions_array_neo[3,2]
# eCO = np.sum(Emissions_array_neo[:,1])-Emissions_array_neo[3,1]
# eH2O = np.sum(Emissions_array_neo[:,3])-Emissions_array_neo[3,3]
# eNOx = np.sum(Emissions_array_neo[:,4])-Emissions_array_neo[3,4] + np.sum(Emissions_array_neo[:,5])-Emissions_array_neo[3,5]
# esoot = np.sum(Emissions_array_neo[:,7])-Emissions_array_neo[3,7]
# eHUC = np.sum(Emissions_array_neo[:,8])-Emissions_array_neo[3,8]
# eSO4 = 0
# TE_LTO_neo= np.array([eCO,eCO2,eH2O,eNOx,esoot,eHUC])
# print('[eCO,eCO2,eH2O,eNOx,esoot,eHUC]')
# print(TE_LTO_neo)
#
# U_ker = climate.number_aircraft_H2()
# h = 450
# ATR_LTO_neo = climate.ATR(h=h, e_CO2=eCO2, e_H2O=eH2O, e_NOx=eNOx, e_soot=esoot, e_sulfate=eSO4, U=U_ker, plot=True)
# print('The average temperature response, A_100, for the LTO of the neo is:', ATR_LTO_neo, '[K]')
#
# #For the cruise phase
# eCO2 = Emissions_array_neo[3,2]
# eCO = Emissions_array_neo[3,1]
# eH2O = Emissions_array_neo[3,3]
# eNOx = Emissions_array_neo[3,4] + Emissions_array_neo[3,5]
# esoot = Emissions_array_neo[3,7]
# eHUC = Emissions_array_neo[3,8]
# eSO4 = 0
# TE_cruise_neo = np.array([eCO,eCO2,eH2O,eNOx,esoot,eHUC])
# print(TE_cruise_neo)
# h = 11600
# ATR_cruise_neo = climate.ATR(h=h, e_CO2=eCO2, e_H2O=eH2O, e_NOx=eNOx, e_soot=esoot, e_sulfate=eSO4, U=U_ker, plot=True)
# print('The average temperature response, A_100, for the cruise of the neo is:', ATR_cruise_neo, '[K]')
#
# TE_neo = TE_LTO_neo + TE_cruise_neo
# ATR_HACK = climate.ATR(h=5800, e_CO2=TE_neo[1], e_H2O=TE_neo[2], e_NOx=TE_neo[3], e_soot=TE_neo[4], e_sulfate=eSO4, U=U_ker, plot=True)
# print('The average temperature response, A_100, for the  HACK is:', ATR_HACK, '[K]')

# '''Storing Neo emissions'''
#
# file = open('Emissions_neo', 'w')
# file.write(str(Emissions_array_neo) + 'ATR during LTO' + str(ATR_LTO_neo)+ 'ATR during cruise' + str(ATR_cruise_neo))
# file.close()

# TE_LTO_neo= np.array([2.85304569e+01, 4.27433727e+03, 1.80216600e+03, 2.35383622e+02,1.09030998e-01, 1.05104692e-01])
# TE_cruise_neo = np.array([6.24238542e+01, 2.79690164e+04, 1.17315741e+04, 1.38232676e+03,1.36566020e-01, 6.82830098e-01])
# ATR_LTO_neo = 0.003257
# ATR_cruise_neo = 0.1680

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
        weight = Compute_weight()
        weight.weight_break_down_HACK()
        perf.mission_profile(T.durations,aircraft[1])
        m_5 = perf.fs_h2_cruise[-1] + perf.fs_k_cruise[-1] + weight.MZFW_HACK
        m_4 = perf.fs_h2_climb[-1] + perf.fs_k_climb[-1] + weight.MZFW_HACK
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

    # file = open('Emissions_HACK', 'w')
    # file.write(str(Emissions_array_HACK) + 'ATR during LTO' + str(ATR_LTO_HACK) + 'ATR during cruise' + str(ATR_cruise_HACK))
    # file.close()

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
# perf = performance()
# weight = Compute_weight()
# weight.weight_break_down_HACK()
# perf.mission_profile(T.durations,aircraft[1])
# m_5 = perf.fs_h2_cruise[-1] + perf.fs_k_cruise[-1] + weight.MZFW_HACK
# m_4 = perf.fs_h2_climb[-1] + perf.fs_k_climb[-1] + weight.MZFW_HACK
# Range_HACK = performance().Range(ae.L_D_ratio_HACK,m_5/m_4,7.48422*10**-6)
#
# if Range_HACK< 3200*10**3:
#
#     print('Harmonic Range is:',Range_HACK)
#     print('Requirement of range not met. Change energy split')
# else:
#     print('Range requirement is met')
#     Range_requirement = True

# Emissions_array_HACK = np.array([[1.09470261e-29, -1.07880063e-28, -2.08289027e-25,  3.71830040e+02,5.22690132e-02,  1.61523861e-03,  1.61690079e-04,  0.00000000e+00,0.00000000e+00],
#                                  [1.48047946e-15,  2.47318095e-03,  7.14873092e+01,  1.23084305e+02,3.51796955e-01,  7.36592851e-03,  2.06200518e-04,  3.35684411e-04,3.35684411e-03],
#                                  [5.63429496e-14,  3.95458917e-02,  1.21509547e+02,  2.08927530e+02,3.63203861e+00,  3.52819995e-02,  2.74284849e-03,  5.71697266e-04,5.71697266e-03],
#                                  [4.36250306e+00,  5.75904202e-01,  1.93117131e-02,  1.92193626e+01,3.89798408e-17,  1.45220486e-14,  8.75946702e+02,  5.65983158e-02,2.82991579e-01],
#                                  [4.59563163e-16,  1.67435003e-03,  7.10990312e+01,  1.22444015e+02,1.09032174e-02,  1.97628952e-04,  1.49401515e-04 , 6.67737257e-04,3.33868629e-03],
#                                  [2.96528131e-28, -8.74485923e-27, -5.77667200e-24,  1.04544534e+03,1.72552027e-01,  5.30907128e-03,  4.74618146e-04,  0.00000000e+00,0.00000000e+00]])
# Emissions_array_HACK = np.array([[0.0,0.0,0.0, 3.71830090e+02,1.71724585e-01,  5.30720233e-03,  1.61615557e-04,  0.00000000e+00,0.00000000e+00],
#                         [.46552547e-15,  2.45755906e-03,  7.14873338e+01,  1.23084786e+02,9.18585857e-01,  1.92507151e-02,  2.05040423e-04,  3.35684411e-04,3.35684411e-03],
#                         [5.17405799e-14,  3.77747186e-02,  1.21512342e+02,  2.08946913e+02,1.07831124e+01,  1.05224900e-01,  2.63334891e-03,  5.71697266e-04,5.71697266e-03],
#                         [1.61746677e+01,  2.92725400e+01,  7.98918511e-01,  1.91311316e+02,8.71145061e-16,  6.91415269e-13,  8.62472332e+02,  5.65983158e-02,2.82991579e-01],
#                         [4.58333835e-16,  1.67192138e-03,  7.10990350e+01,  1.22444100e+02,1.26993623e-01,  2.30254391e-03,  1.49208681e-04, 6.67737257e-04,3.33868629e-03],
#                         [0.0, 0.0, 0.0,  1.04544551e+03,5.44846914e-01,  1.67656632e-02,  4.74353881e-04, 0.00000000e+00,.00000000e+00]])
# Emissions_array_HACK = np.array([[0.0,0.0,0.0, 3.71830090e+02,1.15435702e+00,  3.57021390e-02,  1.60920563e-04,  0.00000000e+00,0.00000000e+00],
#                         [ 1.37999240e-15,  2.36780384e-03,  7.14873338e+01,  1.23084786e+02,4.23726084e+00,  8.92696844e-02,   1.98360280e-04,  3.35684411e-05,3.35684411e-03],
#                         [4.40489443e-14,  3.77747186e-02,  1.21512342e+02,  2.08946913e+02,1.07831124e+01,  1.05224900e-01,  2.63334891e-03,  5.71697266e-05,5.71697266e-03],
#                         [1.33152992e-14,  1.23914231e-01,  6.02674212e+03,  1.03794728e+04,2.33774233e+01,  3.48750321e-01,  1.15407759e-02,  5.65983158e-03,2.82991579e-01],
#                         [4.58333835e-16,  1.67192138e-03,  7.10990350e+01,  1.22444100e+02,1.26993623e-01,  2.30254391e-03,  1.49208681e-04, 6.67737257e-04,3.33868629e-03],
#                         [0.0, 0.0, 0.0,  1.04544551e+03,5.44846914e-01,  1.67656632e-02,  4.74353881e-04, 0.00000000e+00,.00000000e+00]])

# '''ASSESSING THE CLIMATE IMPACT OF HACK'''
# #For the LTO cycle [taxi-out,take-off,clim-out(until 900m),descend(after 900m),landing,taxi-in]
# eCO2 = np.sum(Emissions_array_HACK[:,2])-Emissions_array_HACK[3,2]
# eCO = np.sum(Emissions_array_HACK[:,1])-Emissions_array_HACK[3,1]
# eH2O = np.sum(Emissions_array_HACK[:,3])-Emissions_array_HACK[3,3]
# eNOx = np.sum(Emissions_array_HACK[:,4])-Emissions_array_HACK[3,4] + np.sum(Emissions_array_HACK[:,5])-Emissions_array_HACK[3,5]
# esoot = np.sum(Emissions_array_HACK[:, 7]) - Emissions_array_HACK[3, 7]
# eHUC = np.sum(Emissions_array_HACK[:, 8]) - Emissions_array_HACK[3, 8]
# eSO4 = 0
# TE_LTO_HACK = np.array([eCO, eCO2, eH2O, eNOx, esoot, eHUC])
# print('[eCO,eCO2,eH2O,eNOx,esoot,eHUC]')
# print(TE_LTO_HACK)
# U_H2 = climate.number_aircraft_H2()
# h = 450
# ATR_LTO_HACK = climate.ATR(h=h, e_CO2=eCO2, e_H2O=eH2O, e_NOx=eNOx, e_soot=esoot, e_sulfate=eSO4, U=U_H2, plot=True)
# print('The average temperature response, A_100, for the LTO of the HACK is:', ATR_LTO_HACK, '[K]')
#
# #For the cruise phase
# eCO2 = Emissions_array_HACK[3,2]
# eCO = Emissions_array_HACK[3,1]
# eH2O = Emissions_array_HACK[3,3]
# eNOx = Emissions_array_HACK[3,4] + Emissions_array_HACK[3,5]
# esoot = Emissions_array_HACK[3, 7]
# eHUC = Emissions_array_HACK[3, 8]
# eSO4 = 0
# TE_cruise_HACK = np.array([eCO, eCO2, eH2O, eNOx, esoot, eHUC])
# print(TE_cruise_HACK)
# h = 11600
# ATR_cruise_HACK = climate.ATR(h=h, e_CO2=eCO2, e_H2O=eH2O, e_NOx=eNOx, e_soot=esoot, e_sulfate=eSO4, U=U_H2, plot=True)
# print('The average temperature response, A_100, for the cruise of the HACK is:', ATR_cruise_HACK, '[K]')
#
# TE_HACK = TE_LTO_HACK + TE_cruise_HACK
# ATR_HACK = climate.ATR(h=5800, e_CO2=TE_HACK[1], e_H2O=TE_HACK[2], e_NOx=TE_HACK[3], e_soot=TE_HACK[4], e_sulfate=eSO4, U=U_H2, plot=True)
# print('The average temperature response, A_100, for the  HACK is:', ATR_HACK, '[K]')
# '''Storing hack emissions'''
#
# file = open('Emissions_HACK', 'w')
# file.write(str(Emissions_array_HACK) + 'ATR during LTO' + str(ATR_LTO_HACK) + 'ATR during cruise' + str(ATR_cruise_HACK))
# file.close()
#
# 'Checking cruise emissions'
# if (TE_cruise_HACK[0]<0.3 *TE_cruise_neo[0]) and  (TE_cruise_HACK[3]<0.3 *TE_cruise_neo[3]) and (TE_cruise_HACK[4]<0.3 *TE_cruise_neo[4]) and (TE_cruise_HACK[5]<0.3 *TE_cruise_neo[5]):
#     print('Emissions requirement met during cruise')
#     Emissions_requirement_cruise = True
# else:
#     print('Emissions requirement during cruise is not met. Change Energy split')
#
# if (TE_LTO_HACK[0]<0.5 *TE_LTO_neo[0]) and  (TE_LTO_HACK[3]<0.5 *TE_LTO_neo[3]) and (TE_LTO_HACK[4]<0.5 *TE_LTO_neo[4]) and (TE_LTO_HACK[5]<0.5 *TE_LTO_neo[5]):
#     print('Emissions requirement met during LTO')
#     Emissions_requirement_LTO = True
# else:
#     print('Emissions requirement during LTO is not met. Change Energy split')
#
# if (ATR_cruise_HACK<0.3*ATR_cruise_neo) and (ATR_LTO_HACK<0.3*ATR_LTO_neo):
#     ATR_requirement=True
#     print('Climate_change requirement is met.')
# else:
#     print('Climate_change requirement is not met. Change Energy split')

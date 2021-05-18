

LHV_k = 43.2
LHV_H2 = 119.96
og_SFC = 14.4

Max_Engine_Thrust =  155.6 #kN
#m_flow_air = 37.7  #kg/s using a hard coded 0.27 equivalence ratio CFM567
m_flow_air = 31.75 #using a hard coded 0.3 value for equivalence ratio (e.g from the A330)

T_frac_Taxi = 0.124
T_frac_TO = 0.772
T_frac_Climb = 0.77
T_frac_Cruise = 0.29
T_frac_Descent = 0.28
T_frac_Arr_Roll_Taxi = 0.068

FF_Taxi = 0.1775
FF_T0 = 1.875
FF_Climb = 1.9
FF_Cruise = 0.5
FF_Descent = 0.45
FF_Arrival_Roll_Taxi = 0.145

def sfc(H2_frac,Ker_frac):
    SFC_Ker = og_SFC*Ker_frac
    SFC_H2 = og_SFC*LHV_k*H2_frac/LHV_H2
    SFC_total = SFC_H2 + SFC_Ker


    thrust = Max_Engine_Thrust * T_frac_Cruise #_____ HERE COMPLETE WITH WHAT PHASE YOU WANT e.g. T_frac_Taxi
    m_flow_fuel = (SFC_total * thrust)/1000

    stoichiometric_ratio = H2_frac * 34.3 + Ker_frac * 15.66

    equivalence_ratio = (m_flow_fuel * stoichiometric_ratio)/m_flow_air

    #FAR = m_flow_fuel/m_flow_air

    return (SFC_total, m_flow_fuel,m_flow_air,equivalence_ratio)#equivalence_ratio)



print(sfc(0,1))







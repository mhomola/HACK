import numpy as np
import matplotlib.pyplot as plt
import fuel_constants as fc

LHV_k = 43.2
LHV_H2 = 119.96
og_SFC = 14.4

Max_Engine_Thrust =  155.6 #kN
m_flow_air = 31.75 #using a hard coded 0.3 value for equivalence ratio (e.g from the A330)

#Fuel fraction
FF_Taxi = 0.1775
FF_T0 = 1.875
FF_Climb = 1.9
FF_Cruise = 0.5
FF_Descent = 0.45
FF_Arrival_Roll_Taxi = 0.145

#Thrust
T_frac_Taxi = 0.124
T_frac_TO = 0.772
T_frac_Climb = 0.77
T_frac_Cruise = 0.29
T_frac_Descent = 0.28
T_frac_Arr_Roll_Taxi = 0.068

t_idle = 1229
t_taxi = 26*60 #[s]
t_taxi_out = 17.94*60
t_taxi_in = 8.06*60
t_climb = 25*60
t_descend = 25*60

V_ker = fc.fuel_capacity_a320neo
m_ker_all = fc.fuel_capacity_a320neo*fc.k_d/1000

E_tot = m_ker_all*LHV_k
E_h2 = (1/3)*E_tot

m_h2 = (E_h2/LHV_H2)

sfc_ker = 14.4*0.9 #Altitude has to be added!!!
sfc_h2 = sfc_ker*fc.k_ed/fc.H2_ed

m_ker = (2/3)*m_ker_all*0.98 #considering APU

def stoch_rat(H2_frac):
    return H2_frac * 34.3 + (1-H2_frac) * 15.66

def sfc(H2_frac):
    return sfc_h2*H2_frac+(1-H2_frac)*sfc_ker

def eq_r(m_flow_fuel,m_flow_air,stoichiometric_ratio):
    equivalence_ratio = (m_flow_fuel * stoichiometric_ratio) / m_flow_air
    return equivalence_ratio

def fuel_flow(T_frac, SFC_total):
    thrust = Max_Engine_Thrust * T_frac #_____ HERE COMPLETE WITH WHAT PHASE YOU WANT e.g. T_frac_Taxi
    return (SFC_total * thrust)/1000

#take-off H2 volume
ff_idle = 0.033131044
h2m_idle = 2*ff_idle*t_idle

ff_taxi_o = fuel_flow(T_frac_Taxi,sfc_h2)
h2m_taxi_o = 2*ff_taxi_o*t_taxi_out

ff_taxi_i = fuel_flow(T_frac_Arr_Roll_Taxi,sfc_h2)
h2m_taxi_i = 2*ff_taxi_i*t_taxi_in

h2_m_rem = m_h2-(h2m_idle+h2m_taxi_o+h2m_taxi_i)

r = h2_m_rem/(h2_m_rem+m_ker)

#masses of kerosene
mk_init = m_ker
ff_climb_k = 2*fuel_flow(T_frac_Climb,sfc_ker*(1-r))
ff_cruise_k = 2*fuel_flow(T_frac_Cruise,sfc_ker*(1-r))
ff_desc_k = 2*fuel_flow(T_frac_Descent,sfc_ker*(1-r))



mk_climb = mk_init - ff_climb_k*t_climb
mk_cruise = mk_init - ff_climb_k*t_climb - ff_desc_k*t_descend #mk_climb - 2*ff_cruise_k*t_cruise
mk_desc = mk_cruise - ff_desc_k*t_descend
# ---

ff_climb_h2 = r*ff_climb_k #fuel_flow(T_frac_Climb,sfc_h2*r)
ff_cruise_h2 = r*ff_cruise_k #fuel_flow(T_frac_Cruise,sfc_h2*r)
ff_desc_h2 = r*ff_desc_k #fuel_flow(T_frac_Descent,sfc_h2*r)
#ff_climb = fuel_flow(T_frac_Climb,s)

h2m_climb = 2*ff_climb_h2*t_climb
h2m_desc = 2*ff_desc_h2*t_descend

h2_m_rem2 = h2_m_rem - (h2m_climb+h2m_desc)

t_cruise = h2_m_rem2/ff_cruise_h2

#representative times in flight
t0 = 0
t1 = t_idle
t2 = t1 + t_taxi_out
t3 = t2 + t_climb
t4 = t3 + t_cruise
t5 = t4 + t_descend
t6 = t5 + t_taxi_in

#masses of hydrogen in the flight
mh2_idle = m_h2 - h2m_idle
mh2_taxi_o = mh2_idle - h2m_taxi_o
mh2_climb = mh2_taxi_o - h2m_climb
mh2_cruise = mh2_climb - h2_m_rem2
mh2_desc = mh2_cruise - h2m_desc
mh2_taxi_i = mh2_desc - h2m_taxi_i

print(1-r)

t = [t0,t1,t2,t3,t4,t5,t6]
fuel_h2 = [m_h2,mh2_idle, mh2_taxi_o, mh2_climb, mh2_cruise, mh2_desc, mh2_taxi_i]
fuel_k = [m_ker, mk_init, mk_init, mk_climb, mk_cruise, mk_desc, mk_desc]

plt.plot(t,fuel_h2,label = 'H2')
plt.plot(t,fuel_k,label = 'Kerosene')
plt.legend()
plt.show()
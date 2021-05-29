import numpy as np
import matplotlib.pyplot as plt
import fuel_constants as fc

LHV_k = 43.2
LHV_H2 = 119.96
og_SFC = 14.4

Max_Engine_Thrust = 155.6 #kN
m_flow_air = 31.75 #using a hard coded 0.3 value for equivalence ratio (e.g from the A330)

#Fuel fraction
FF_Taxi = 0.1775
FF_T0 = 1.875
FF_Climb = 1.9
FF_Cruise = 0.5
FF_Descent = 0.45
FF_Arrival_Roll_Taxi = 0.145
F_cruise = 18.8

#Thrust
T_frac_Taxi = 0.124
T_frac_TO = 0.772
T_frac_Climb = 0.55
T_frac_Cruise = 0.29
T_frac_Descent = 0.28
T_frac_Arr_Roll_Taxi = 0.068

t_idle = 1229
t_taxi = 26*60 #[s]
t_taxi_out = 17.94*60
t_taxi_in = 8.06*60
t_climb = 23*60
t_descend = 25*60

V_ker = fc.fuel_capacity_a320neo
m_ker_all = fc.fuel_capacity_a320neo*fc.k_d/1000

E_tot = m_ker_all*LHV_k
E_h2 = (1/3)*E_tot

m_h2 = (E_h2/LHV_H2)
print(m_h2)

sfc_ker = 14.4*0.9 #Altitude has to be added!!!
sfc_h2 = sfc_ker*fc.k_ed/fc.H2_ed


print(sfc_h2/(sfc_ker+sfc_h2))

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

print((h2m_idle+h2m_taxi_o+h2m_taxi_i))

h2_m_rem = m_h2-(h2m_idle+h2m_taxi_o+h2m_taxi_i)
r = h2_m_rem*LHV_H2/(h2_m_rem*LHV_H2+m_ker*LHV_k)
print(r)
#r = h2_m_rem/(h2_m_rem+m_ker)


ff_climb_h2 = fuel_flow(T_frac_Climb,sfc_h2*r)
ff_cruise_h2_OLD = fuel_flow(T_frac_Cruise,sfc_h2*r)
ff_cruise_h2 = (r * sfc_h2 * F_cruise)/1000
print(ff_cruise_h2_OLD)
print(ff_cruise_h2)
ff_desc_h2 = fuel_flow(T_frac_Descent,sfc_h2*r)
#ff_climb = fuel_flow(T_frac_Climb,s)

h2m_climb = 2*ff_climb_h2*t_climb
h2m_desc = 2*ff_desc_h2*t_descend

h2_m_rem2 = h2_m_rem - (h2m_climb+h2m_desc)

t_cruise = h2_m_rem2/(2*ff_cruise_h2)

#masses of kerosene
mk_init = m_ker
ff_climb_k = fuel_flow(T_frac_Climb,sfc_ker*(1-r))
ff_cruise_k_OLD = fuel_flow(T_frac_Cruise,sfc_ker*(1-r))
ff_cruise_k = ((1-r)*sfc_ker* F_cruise)/1000
ff_desc_k = fuel_flow(T_frac_Descent,sfc_ker*(1-r))

mk_climb = mk_init - 2*ff_climb_k*t_climb
mk_cruise = mk_climb - 2*ff_cruise_k*t_cruise
mk_desc = mk_cruise - 2*ff_desc_k*t_descend
# ---

#representative times in flight
t0 = 0
t1 = t_idle
t2 = t1 + t_taxi_out
t3 = t2 + t_climb
t4 = t3 + t_cruise
t5 = t4 + t_descend
t6 = t5 + t_taxi_in

x=t5-t2
print(x)

print(t2,'s - ',t5)

print(m_h2)
print(h2m_idle+h2m_taxi_i+h2m_taxi_o, 'kilogram')

#masses of hydrogen in the flight
mh2_idle = m_h2 - h2m_idle
mh2_taxi_o = mh2_idle - h2m_taxi_o
mh2_climb = mh2_taxi_o - h2m_climb
mh2_cruise = mh2_climb - h2_m_rem2
mh2_desc = mh2_cruise - h2m_desc
mh2_taxi_i = mh2_desc - h2m_taxi_i


t = [t0,t1,t2,t3,t4,t5,t6]
fuel_h2 = np.array([m_h2,mh2_idle, mh2_taxi_o, mh2_climb, mh2_cruise, mh2_desc, mh2_taxi_i])
fuel_k = np.array([m_ker, mk_init, mk_init, mk_climb, mk_cruise, mk_desc, mk_desc])
fuel_tot = fuel_h2+fuel_k

Eh = 119.96*np.array([0,h2m_idle, h2m_taxi_o, h2m_climb, h2_m_rem2, h2m_desc, h2m_taxi_i])
Ek = 43.2*np.array([0,0,0,2*ff_climb_k*t_climb,2*ff_cruise_k*t_cruise,2*ff_desc_k*t_descend,0])

Em_sum = np.array([])
Eh_sum = np.array([])

for i in range(len(Eh)):
    Em_sum = np.append(Em_sum,sum(Ek[:i+1]))
    Eh_sum = np.append(Eh_sum,sum(Eh[:i+1]))
    print(Eh_sum)

plt.plot(t,fuel_h2,label = 'H2')
plt.plot(t,fuel_k,label = 'Kerosene')
plt.plot(t,fuel_tot,label = 'Total')
plt.ylabel('Mass of fuel on-board [kg]')
plt.xlabel('Time [s]')
plt.legend()
plt.show()

E_tot = Em_sum + Eh_sum

plt.plot(t,Em_sum,label = 'H2')
plt.plot(t,Eh_sum,label = 'Kerosene')
plt.plot(t,E_tot,label = 'Total')
plt.ylabel('Energy consumed over mission [MJ]')
plt.xlabel('Time [s]')
plt.legend()
plt.show()

h2ff = [ff_idle,ff_idle, ff_taxi_o, ff_taxi_o, ff_climb_h2, ff_climb_h2, ff_cruise_h2, ff_cruise_h2, ff_desc_h2, ff_desc_h2, ff_taxi_i, ff_taxi_i]
t2 = [t0,t1,t1,t2,t2,t3,t3,t4,t4,t5,t5,t6]
kff = [0,0,0,0,ff_climb_k,ff_climb_k,ff_cruise_k,ff_cruise_k,ff_desc_k,ff_desc_k,0,0]
plt.plot(t2,h2ff,label = 'H2')
plt.plot(t2,kff,label = 'H2')
plt.legend()
plt.ylabel('Fuel flow [kg/s]')
plt.xlabel('Time [s]')
plt.show()



# Mission fuel fractions
#Hydrogen
fr_h2_idle = h2m_idle/m_h2
fr_h2_taxi_out = h2m_taxi_o/m_h2
fr_h2_taxi_in = h2m_taxi_i/m_h2
fr_h2_climb = h2m_climb/m_h2
fr_h2_cruise = h2_m_rem2/m_h2
fr_h2_desc = h2m_desc/m_h2

print('h2 idle ',h2m_idle)
print('h2 taxi out ',h2m_taxi_o)
print('h2 taxi in ',h2m_taxi_i)
print('h2 climb ',h2m_climb)
print('h2 cruise ',h2_m_rem2)
print('h2 descend ',h2m_desc)
print()

#Kerosene
fr_k_idle = 0
fr_k_taxi_o = 0
fr_k_taxi_i = 0
fr_k_climb = 2*ff_climb_k*t_climb/m_ker
fr_k_cruise = 2*ff_cruise_k*t_cruise/m_ker
fr_k_desc = 2*ff_desc_k*t_descend/m_ker

print('k idle ',0)
print('k taxi out ',0)
print('k taxi in ',0)
print('k climb ', 2*ff_climb_k*t_climb)
print('k cruise ',2*ff_cruise_k*t_cruise)
print('k descend ',2*ff_desc_k*t_descend)
print()

#total
m_total = m_h2+m_ker
fr_t_idle = (h2m_idle+0)/m_total
fr_t_taxi_out = (h2m_taxi_o+0)/m_total
fr_t_taxi_in = (h2m_taxi_i+0)/m_total
fr_t_climb = (h2m_climb+2*ff_climb_k*t_climb)/m_total
fr_t_cruise = (h2_m_rem2+2*ff_cruise_k*t_cruise)/m_total
fr_t_desc = (h2m_desc+2*ff_desc_k*t_descend)/m_total


print(fr_t_idle+fr_t_taxi_out+fr_t_taxi_in+fr_t_climb+fr_t_cruise+fr_t_desc)

print('Fuel fractions')
print('Hydrogen:')
print('Idle: ',fr_h2_idle)
print('Taxi take off: ',fr_h2_taxi_out)
print('Climb: ',fr_h2_climb)
print('Cruise: ',fr_h2_cruise)
print('Descend: ',fr_h2_desc)
print('Taxi landing: ',fr_h2_taxi_in)
print()
print('Kerosene:')
print('Idle: ',fr_k_idle)
print('Taxi take off: ',fr_k_taxi_o)
print('Climb: ',fr_k_climb)
print('Cruise: ',fr_k_cruise)
print('Descend: ',fr_k_desc)
print('Taxi landing: ',fr_k_taxi_i)
print()
print('Total:')
print('Idle: ',fr_t_idle)
print('Taxi take off: ',fr_t_taxi_out)
print('Climb: ',fr_t_climb)
print('Cruise: ',fr_t_cruise)
print('Descend: ',fr_t_desc)
print('Taxi landing: ',fr_t_taxi_in)
import numpy as np
import matplotlib.pyplot as plt

EOCV = 1.05    # [V]
b = 0.067   # [V/decade]
R = 0.11    # [Ohm*cm2]
iloss = 0.00123 # [A*cm2]
m = 0.005   # [V]
n = 2.24    # [cm2/A]
cp = 1004
gamma = 1.4

i = np.arange(0,1.6,0.01)

def Ecell(i):
    return EOCV-b*np.log10((i+iloss)/iloss)-R*i-m*(np.exp(n*i)-1)

def ISA_trop(h):
    Tb = 15 + 273.15
    pb = 101325
    dT = -0.0065

    T = Tb+dT*h
    p = pb*(1-0.0065*(h/Tb))**5.2561

    return T,p

def ISA_strat(h):

    g = 9.81
    R = 287.05

    T11 = ISA_trop(11000)[0]
    p11 = ISA_trop(11000)[1]

    p = p11*np.exp(-(g/(R*T11))*(h-11000))

    return T11,p

def ISA(h):

    if h <= 11000:
        return ISA_trop(h)
    return ISA_strat(h)

def ISA_p(T):
    p0 = 101325
    T0 = 15 + 273.15
    R = 287.05
    lam = 0.0065
    g = 9.81

    return p0*(T/T0)**(g/(lam*R))

def ISA_t(h):
    T0 = 15 + 273.15
    dT = -0.0065

    return T0 + dT * h

def P_com(h,p2):

    T1 = ISA_t(h)

    m = 0.02
    nc = 0.7
    gamma = 1.4
    p1 = ISA_p(T1)
    
    return cp*(T1/nc)*((p2/p1)**((gamma-1)/gamma)-1)*m


print(ISA_p(281.65))

V = np.array([Ecell(j) for j in i])
P = i*V
n_eff = 0.95
th_eff = n_eff*V*(1/1.25)*100

plt.plot(i,V,'k', label='Actual value')
plt.title('Fuel cell voltage as a function of current')
plt.plot(i,[1.23 for j in i],':',label='Max theoretical value')
plt.ylabel('Voltage [V]')
plt.xlabel('Current density [A/cm2]')
plt.legend()
plt.ylim([0,1.4])
plt.show()

plt.plot(i,P,'k')
plt.title('Fuel cell power as a function of current')
plt.ylabel('Power density [W/cm2]')
plt.xlabel('Current density [A/cm2]')
plt.show()

plt.plot(P,th_eff,'k')
plt.title('Fuel cell efficiency as a function of power')
plt.ylabel('Thermal efficiency [-]')
plt.xlabel('Power [W/cm2]')
plt.show()

h = np.arange(0,11001)
P_comp_req = P_com(h,200000)

plt.plot(h/1000,P_comp_req/1000,'k')
plt.title('Compressor power as a function of altitude, m = 1, nc = 0.7')
plt.ylabel('Compressor power [kW]')
plt.xlabel('Altitude [km]')
plt.show()


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

def P_com(h,p2,m):

    T1 = ISA_t(h)

    nc = 0.7
    gamma = 1.4
    g = 9.81
    R = 287.05

    if h <= 11000:
        p1 = ISA_p(T1)
    else:
        p11 = ISA_p(ISA_t(11000))
        p1 = p11*np.exp(-g/(R*ISA_t(11000))*(h-11000))
    
    return cp*(T1/nc)*((p2/p1)**((gamma-1)/gamma)-1)*m, p1

def inverse(Vsearch):

    i = np.arange(0, 1.6, 0.01)
    V = np.array([Ecell(j) for j in i])

    idx = 0
    for j in range(len(V)):
        if abs(V[j]-Vsearch)/Vsearch < 0.0015:
            idx = j
            break
    return i[idx]

print(ISA_p(281.65))

A = 0.5*100*100

i = np.arange(0,1.6,0.01)

V = np.array([Ecell(j) for j in i])
P = i*V
n_eff = 0.95
th_eff = n_eff*V*(1/1.25)*100

# Iact = i*A
# Vact = np.array([Ecell(j) for j in Iact])
# Pact = Iact*Vact

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
plt.ylabel('Thermal efficiency [%]')
plt.xlabel('Power [W/cm2]')
plt.show()


plt.plot(V,th_eff,'k')
plt.title('Fuel cell efficiency as a function of voltage')
plt.ylabel('Thermal efficiency [%]')
plt.xlabel('Voltage [V]')
plt.show()

plt.plot(i*A,i*A*V*160,'k')
plt.title('Actual fuel cell power as a function of current')
plt.ylabel('Power density [W/cm2]')
plt.xlabel('Current density [A/cm2]')
plt.show()

h = np.arange(0,13001)
# P_comp_req1 = P_com(h,200000,0.05)
# P_comp_req2 = P_com(h,200000,0.1)
# P_comp_req3 = P_com(h,200000,0.15)
# P_comp_req4 = P_com(h,200000,0.2)

P_comp_req1 = np.array([P_com(i,200000,0.05)[0] for i in h])
P_comp_req2 = np.array([P_com(i,200000,0.1)[0] for i in h])
P_comp_req3 = np.array([P_com(i,200000,0.15)[0] for i in h])
P_comp_req4 = np.array([P_com(i,200000,0.2)[0] for i in h])
p_check = np.array([P_com(i,200000,0.2)[1] for i in h])

plt.plot(h/1000,P_comp_req1/1000,'#6f83e3', label = 'm = 0.05')
plt.plot(h/1000,P_comp_req2/1000,'#3c4a8c', label = 'm = 0.1')
plt.plot(h/1000,P_comp_req3/1000,'#1d2759', label = 'm = 0.15')
plt.plot(h/1000,P_comp_req4/1000,'#060b24', label = 'm = 0.2')
plt.legend()
plt.title('Compressor power as a function of altitude, nc = 0.7')
plt.ylabel('Compressor power [kW]')
plt.xlabel('Altitude [km]')
plt.show()

plt.plot(h, p_check)
plt.show()


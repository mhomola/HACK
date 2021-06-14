import numpy as np
import matplotlib.pyplot as plt
import FCParameters as fc
import math as m

cp = fc.cp

cost = 85 #EUR/kWh
spec_m = 300 #Wh/kg
spec_V = 280 #Wh/l

def ISA_trop(h):
    Tb = 15 + 273.15
    pb = 101325
    dT = -0.0065

    T = Tb + dT * h
    p = pb * (1 - 0.0065 * (h / Tb)) ** 5.2561

    return T, p


def ISA_strat(h):
    g = 9.81
    R = 287.05

    T11 = ISA_trop(11000)[0]
    p11 = ISA_trop(11000)[1]

    p = p11 * np.exp(-(g / (R * T11)) * (h - 11000))

    return T11, p


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

    return p0 * (T / T0) ** (g / (lam * R))


def ISA_t(h):
    T0 = 15 + 273.15
    dT = -0.0065

    return T0 + dT * h


def P_com(h, p2, m):

    T1 = ISA_t(h)
    nc = 0.75
    gamma = 1.4
    g = 9.81
    R = 287.05

    # if h <= 11000:
    #     p1 = ISA_p(T1)
    if h <= 1800:
        p1 = ISA_p(T1)
    else:
        p1 = ISA_p(ISA_t(1800))

        # p11 = ISA_p(ISA_t(11000))
        # p1 = p11 * np.exp(-g / (R * ISA_t(11000)) * (h - 11000))

    return cp * (T1 / nc) * ((p2 / p1) ** ((gamma - 1) / gamma) - 1) * m, p1

def Req_power(mair):

    # Power in diferent mission stages
    scale = 180/150
    p_ground = 195.41*scale*1000
    p_take_off = 90.13*scale*1000
    p_climb = 223.5*scale*1000
    p_cruise = 275.6*scale*1000
    p_descend = 233.5*scale*1000
    # loiter = 231.44*scale*1000
    p_landing = 96.35*scale*1000


    # Time of different mission stages
    t_ground = 20*60
    t_take_off = 60
    t_climb = 20*60
    t_cruise = 27740
    t_descend = 18*60
    # t_loiter = 0
    t_landing = 3*60

    # Altitude of different mission stages
    a_ground = 0
    a_climb = np.arange(0,11000,1)
    a_cruise = 11000

    pFC = 101325

    # Power of different mission stages
    pc_ground = P_com(a_ground,pFC,mair)[0]
    pc_climb = np.average([P_com(a,pFC,mair)[0] for a in a_climb])
    pc_cruise = P_com(a_cruise,pFC,mair)[0]

    power = np.array([p_ground+pc_ground,p_take_off+pc_ground,p_climb+pc_climb,p_cruise+pc_cruise,p_descend+pc_climb,p_landing+pc_ground])
    time = np.array([t_ground,t_take_off,t_climb,t_cruise,t_descend,t_landing])

    average = sum(power*time)/sum(time)

    #print('Average power [W]: ', average)

    avg_arr = average*np.ones(len(power+1))
    difference = avg_arr - power

    time_all = np.concatenate(([0],time))
    power_all = np.concatenate(([0],difference*time))

    time_tot = [0]
    power_tot = [0]
    power_bat = []

    for i in range(len(time_all)-1):
        time_tot.append((time_tot[i]+time_all[i+1]/60))
        power_tot.append(power_tot[i]+power_all[i+1])

    power_tot_upd = []
    for i in range(len(time_all)):
        power_tot_upd.append(power_tot[i]-min(power_tot))
        power_bat.append(average-power_all[i])

    # #plt.plot(time_tot, power_tot)
    # plt.plot(time_tot, power_tot_upd)
    # #plt.plot(time_tot, avg_arr)
    # plt.ylabel('Energy stored in battery [kJ]')
    # plt.xlabel('Time [min]')
    # plt.show()

    #print('Maximum energy stored in battery: ', max(power_tot_upd), ' kJ')

    return average, max(power_tot_upd), time_tot, power_tot_upd

# Mass flow required
Mo2 = 31.998    #g/mol
Mh2 = 2.02      #g/mol
MN2 = 28        #g/mol
lam = 2         #stochiometric ratio
#Pelec = fc.SF*fc.Pele   #FC power output
Farad = 96485   #Faraday's constant
Vc = fc.Vc      #Cell voltage
po2 = 0.22
pN2 = 0.78
Mair = po2*Mo2+pN2*MN2
W_rat = Mo2*po2/Mair

def mfoAir(Pel):
    mfo2 = (Mo2*Pel*lam)/(4*Vc*Farad)
    mfoAir = mfo2/(1000*W_rat) # Required mass flow of air
    return mfoAir

def mfoH2(Pel):
    mfo2 = (Mh2*Pel)/(2*Vc*Farad*1000)
    return mfo2

def toWh(val):
    return val*0.000277778

#if __name__ == '__main__':

def final_bat_size():

    FC_power = 0
    FC_power_new, bat_E, time_tot1, power_tot1 = Req_power(0)

    mair = mfoAir(FC_power_new)
    mh2 = mfoH2(FC_power_new)

    print('Original estimate without compressor')
    print('Air mass flow: ',round(mair*1000,2),' kg/s')
    print('Hydrogen mass flow: ', round(mh2*1000,2), ' kg/s')
    print(FC_power_new / 1000, " kW")
    print(bat_E / 1000000, " MJ")
    bat_orig = toWh(bat_E)
    print('Battery cost: ', round(0.001*bat_orig*cost), " EUR")
    print('Battery mass: ', round(bat_orig/spec_m), " kg")
    print('Battery volume: ', round(bat_orig/spec_V), " l")
    print()

    while(abs(FC_power-FC_power_new)/FC_power_new >= 0.05):
        FC_power = FC_power_new
        mAir = mfoAir(FC_power)
        FC_power_new, bat_E, time_tot2, power_tot2 = Req_power(mAir)

    bat_E = bat_E/0.85
    mair = mfoAir(FC_power_new)
    mh2 = mfoH2(FC_power_new)
    bat_orig = toWh(bat_E)

    Vbcell = 3.7 # V
    Ahbcell = 2
    b_voltage = fc.FCV
    b_cells_s = m.ceil(b_voltage/Vbcell)
    b_Ah = bat_orig/b_voltage
    b_cells_par = m.ceil(b_Ah/Ahbcell)
    b_tot = b_cells_s*b_cells_par

    print('New estimate with compressor')
    print('Air mass flow: ', round(mair*1000,2), ' g/s')
    print('Hydrogen mass flow: ', round(mh2*1000,2), ' g/s')
    print(FC_power_new/1000, " kW")
    print(bat_E/1000000, " MJ")
    print()
    print('Battery cells in series: ',b_cells_s)
    print('Battery cells in parallel: ',b_cells_par)
    print('Battery cells - total: ',b_tot)
    print('Battery cost: ', round(0.001 * bat_orig * cost), " EUR")
    print('Battery mass: ', round(bat_orig / spec_m), " kg")
    print('Battery volume: ', round(bat_orig / spec_V), " l")
    print()

    plt.plot(time_tot1, power_tot1, label = 'No compressor')
    #plt.plot(time_tot2, power_tot2, 'k', label = 'Compressor')
    # plt.plot(time_tot, avg_arr)
    plt.legend()
    plt.ylabel('Energy stored in battery [kJ]')
    plt.xlabel('Time [min]')
    plt.show()

    return FC_power_new, round(0.001 * bat_orig * cost), round(bat_orig / spec_m), round(bat_orig / spec_V), mh2


def compressor_plotter():

    h = np.arange(0,11000,1)
    mair = 0.3
    pcom = [P_com(i,101325,mair)[0] for i in h]

    plt.plot(h,pcom,'k')
    plt.title('Compressor power, m_air = 0.3 kg/s')
    plt.ylabel('Power required [W]')
    plt.xlabel('Altitude [m]')
    plt.show()

def powers_req():

    scale = 18/15
    p_ground = 195.41 * scale * 1000
    p_take_off = 90.13 * scale * 1000
    p_climb = 223.5 * scale * 1000
    p_cruise = 275.6 * scale * 1000
    p_descend = 233.5 * scale * 1000
    # loiter = 231.44*scale*1000
    p_landing = 96.35 * scale * 1000

    print()
    print('Ground: ',round(p_ground*0.001),'W')
    print('Take-Off: ',round(p_take_off*0.001),'W')
    print('Climb: ',round(p_climb*0.001),'W')
    print('Cruise: ',round(p_cruise*0.001),'W')
    print('Descend: ', round(p_descend*0.001),'W')
    print('Landing: ',round(p_landing*0.001),'W')




from Subsystem_design.common_constants import Constants
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics


from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from cg_excursion import cg_range
from sklearn.linear_model import LinearRegression

import numpy as np
import math as m
import matplotlib.pyplot as plt

# class xplot(Constants):
#     def __init__(self):
#         super().__init__()
#         self.v_approach = 80    #[m/s]
#         ## -- from common constants ##

# to be imported
#cruise mach = self.M
#b = self.b
#c = self.mac
#xlemac = self.x_LEMAC
#S = self.S
#A = b/c
#sweep = self.sweep_025
#a = self.a_c
#v = self.V_cruise
#Snet = self.S - self.S_h
#taper = self.taper_out
#bh = self.b_h
#ch = tbd
#Ah = bh/ch
#sweep_h = self.sweep_LE_h
#deda = tbd
#lh = self.c_r_h
#sh = self.S_h
#Vratio = 0.85 #tbc
#SM = 0.05 # stability margin
#a0 = tbd
#Mapp = Vapp/a0*m.sqrt(Vratio)
#bf =
#lf =
#hf =


M = self.M  # Mach number at cruise
Ml = 0.207  # Mach number at landing assuming 137 kts at sea level
Aw = self.b**2/self.S  # wing aspect ratio
bh = self.b_h
ch = tbd
Ah = bh/ch # h tail aspect ratio
S = self.S  # wing surface area [m^2]
Snet = 66.96  # wing surface area minus part in the fuselage [m^2]
Swf = 19.51  # TE flap area [m^2]
n = 0.95  # wing efficiency (derived from slides)
sweep_25c = self.sweep_025 * pi / 180  # 25% chord sweep angle for the wing [rad]
sweep_25c_h = self.sweep_025_h * pi / 180  # 25% chord sweep angle for the h tail [rad]
Cr = self.c_root  # Root chord wing [m]
Cr_h = self.c_r_h  # Root chord h tail [m]
taper_ratio_w = 0.24  # Taper ratio wing a320neo from online
taper_ratio_h = self.taper_h  # Taper ratio h tail
b = self.b  # Wing span wing [m]
bf = 3.95  # width of fuselage or component of wing span inside fuselage [m]
hf = 4.14  # height of fuselage [m]
bh = 12.6  # Wing span h tail [m]
lfn = 11.88  # position of wing at the fuselage from nose [m]
lh = 12.35  # tail arm [m]
lf = 26.5  # fuselage length [m]
MAC = self.mac  # [m]
mtv = 5  # distance between downwash lines of the wing and h tail [m]
Vland = 125 * 0.514444  # landing speed [m/s]
MLW = 38566  # [kg]
MTOW = 42185  # [kg]
g = 9.8065  # [m/s^2]
rho = 1.225  # see level density [kg/m^3]
rho_cruise = self.rho_c

ln_in = 3.29
ln_out = 0
bn = 1.93
beta = sqrt(1 - M ** 2)
beta_low = sqrt(1 - Ml ** 2)

lemac_des = 11.63

""" DEFINE FUNCTIONS """
# Sweep angles


def sweep(xc, tail=False):
    if tail == False:
        sweep_LE = atan(tan(sweep_25c) + 0.5 * Cr * (1 - taper_ratio_w) / b)
        return atan(tan(sweep_LE) - 2 * xc * Cr * (1 - taper_ratio_w) / b)
    elif tail == True:
        sweep_LE = atan(tan(sweep_25c_h) + 0.5 * Cr_h * (1 - taper_ratio_h) / bh)
        return atan(tan(sweep_LE) - 2 * xc * Cr_h * (1 - taper_ratio_h) / bh)


# Calculating the wing lift curve slope
def CLaw(beta):
    return (2 * pi * Aw) / (2 + sqrt(4 + (1 + (tan(sweep(0.5, False)) / beta) ** 2) * (Aw * beta / n) ** 2))


# Calculating the horizontal tail lift curve slope
def CLah(beta):
    return (2 * pi * Ah) / (2 + sqrt(4 + (1 + (tan(sweep(0.5, True)) / beta) ** 2) * (Ah * beta / n) ** 2))


# Calculating the lift curve slope of the whole aircraft minus the tail
def CLaAh(beta):
    return CLaw(beta) * (1 + 2.15 * bf / b) * Snet / S + 0.5 * pi * bf * bf / S


# Aerodynamic center of wing + fuselage
def Xac_w(landing=False):
    if landing == False:
        return 0.26
    elif landing == True:
        return 0.28


def fus_cont_1(beta):
    return (-1.8 * bf * hf * lfn) / (CLaAh(beta) * S * MAC)


def Xac_n(beta):
    return 2 * (-4 * ln_in * bn ** 2) / (S * MAC * CLaAh(beta)) + 2 * (-4 * ln_out * bn ** 2) / (S * MAC * CLaAh(beta))


def Xac(beta, landing=False):
    return Xac_w(landing) + fus_cont_1(beta) + fus_cont_2 + Xac_n(beta)


# For plotting
def stability(Xcg):
    denominator = (CLah(beta) * (1 - e_grad) * lh * 1) / (CLaAh(beta) * MAC)
    return Xcg / denominator - (Xac_stab - 0.05) / denominator, Xcg / denominator - (Xac_stab) / denominator


def controllability(Xcg):
    denominator = (CLh * lh * 1) / (CLAh * MAC)
    return Xcg / denominator + ((Cmac / CLAh) - Xac_contr) / denominator


# Find design condition
def regression(cg_list, regr_list):
    regr = LinearRegression().fit(cg_list.reshape(-1, 1), regr_list)
    m = regr.coef_[0]
    b = regr.intercept_
    return m, b


def F_min_cg(lemac):
    if lemac >= regr1_min[-1]:
        min_cg = (lemac - lemac1_min_b) / lemac1_min_m
    elif lemac < regr1_min[-1]:
        min_cg = (lemac - lemac2_min_b) / lemac2_min_m

    return min_cg


def F_max_cg(lemac):
    if lemac >= regr1_max[-1]:
        max_cg = (lemac - lemac1_max_b) / lemac1_max_m
    elif lemac < regr1_max[-1]:
        max_cg = (lemac - lemac2_max_b) / lemac2_max_m

    return max_cg


""" FIND Cmac """
# Aerodynamic center of wing + fuselage
CG = S / b
fus_cont_2 = tan(sweep_25c) * (0.273 * bf * CG * (b - bf)) / ((1 + taper_ratio_w) * (b + 2.15 * bf) * MAC ** 2)
Xac_stab = Xac(beta, False)
Xac_contr = Xac(beta_low, True)

# Calculating the downwash gradient
r = 2 * lh / b
Kea = ((0.1124 + 0.1265 * sweep_25c + 0.1766 * sweep_25c ** 2) / (r ** 2)) + 0.1024 / r + 2
Kea_0 = 0.1124 / r ** 2 + 0.1024 / r + 2
e_grad = (Kea / Kea_0) * ((r / (r ** 2 + mtv ** 2)) * (0.4876 / sqrt(r ** 2 + 0.6319 + mtv ** 2)) +
                          (1 + (1 - sqrt((mtv ** 2) / (1 + mtv ** 2))) * (
                                      (r ** 2) / (r ** 2 + 0.7915 + 5.0734 * mtv ** 2)) ** 0.3113)) * (
                     CLaw(beta) / (Aw * pi))

# Flap contribution to pitching moment
cf_c = 0.3
deltac_cf = 0.81
c_prime_c = 1 + cf_c * deltac_cf
deltaCl_max = 1.6 * c_prime_c
deltaCL_max = 0.9 * deltaCl_max * Swf * cos(sweep(0.7, False)) / S

S_land = S * (1 + (Swf / S) * (c_prime_c - 1))
CL_max_land = (MLW * g) / (0.5 * rho * (S_land) * Vland ** 2)
u1 = 0.21
u2 = 1.05
u3 = 0.275
deltaCm_25c = u2 * (-u1 * deltaCl_max * c_prime_c - (CL_max_land + deltaCl_max * (1 - Swf / S)) * c_prime_c * (
            c_prime_c - 1) / 8) + (0.7 * Aw * u3 * deltaCl_max * tan(sweep_25c)) / (1 + 2 / Aw)
deltaCm_flap = deltaCm_25c - CL_max_land * (0.25 - Xac_contr)
print('deltaCm_flap = ', deltaCm_flap)

# wing contribution to pitching moment
Cm0 = -0.3
deltaCm_wing = Cm0 * ((Aw * (cos(sweep_25c)) ** 2) / (Aw + 2 * cos(sweep_25c)))
print('deltaCm_wing = ', deltaCm_wing)

# fuselage contribution to pitching moment
delta_aoa_0l = -15 * pi / 180  # [rad]
aoa_0Lf = -7.5 * pi / 180 + delta_aoa_0l * Swf / S * cos(sweep(0.7, tail=False))  # [rad]
CLaf = CLaw(beta_low) * S_land / S
CL0f = -aoa_0Lf * CLaf
deltaCm_fuselage = -1.8 * (1 - 2.5 * bf / lf) * (pi * bf * hf * lf * CL0f) / (4 * S * MAC * CLaAh(beta_low))
print('deltaCm_fuselage', deltaCm_fuselage)

# nacelle contribution to pitching moment
deltaCm_nacelle = -0.1
print('deltaCm_nacelle', deltaCm_nacelle)

# total delta cm
Cmac = deltaCm_wing + deltaCm_flap + deltaCm_fuselage + deltaCm_nacelle
print('Cmac = ', Cmac)

CL_cruise = MTOW * g / (0.5 * rho_cruise * S * (420 * .51444) ** 2)

""" other variables needed for controllability line """
CLh = -0.35 * (Ah ** (1 / 3))
CLAh = CL_max_land

""" NEEDED FOR PLOTTING """
Xcg = np.linspace(-0.1, 1, 1000)

# from the potato diagram
lemac_lf = cg_range()[0]
cg_min = cg_range()[1]
cg_max = cg_range()[2]

""" FIND OPTIMAL DESIGN POINT """
contr_m, contr_b = regression(Xcg, controllability(Xcg))
stab_m, stab_b = regression(Xcg, stability(Xcg)[0])

# from the potato diagram
lemac_lf = np.flip(cg_range()[0])
cg_min = np.flip(np.array(cg_range()[1]))
cg_max = np.flip(np.array(cg_range()[2]))

# CG MIN LINE: BREAK INTO TWO SEGMENTS
slope1_min = (lemac_lf[0] - lemac_lf[1]) / (cg_min[1] - cg_min[0])
slope2_min = (lemac_lf[-2] - lemac_lf[-1]) / (cg_min[-1] - cg_min[-2])
regr1_min, cg1_min = np.array([lemac_lf[0]]), np.array([cg_min[0]])
regr2_min, cg2_min = np.array([]), np.array([])

for i in range(len(lemac_lf) - 1):

    if round(((lemac_lf[i] - lemac_lf[i + 1]) / (cg_min[i + 1] - cg_min[i])), 1) < 0.5:
        regr1_min = np.append(regr1_min, lemac_lf[i + 1])
        cg1_min = np.append(cg1_min, cg_min[i + 1])

    elif round(((lemac_lf[i] - lemac_lf[i + 1]) / (cg_min[i + 1] - cg_min[i])), 1) > 0.5:
        regr2_min = np.append(regr2_min, lemac_lf[i])
        cg2_min = np.append(cg2_min, cg_min[i])

regr2_min = np.append(regr2_min, lemac_lf[i + 1])
cg2_min = np.append(cg2_min, cg_min[i + 1])

lemac1_min_m, lemac1_min_b = regression(cg1_min, regr1_min)
lemac2_min_m, lemac2_min_b = regression(cg2_min, regr2_min)

# CG MAX LINE: BREAK INTO TWO SEGMENTS
slope1_max = (lemac_lf[0] - lemac_lf[1]) / (cg_max[1] - cg_max[0])
slope2_max = (lemac_lf[-2] - lemac_lf[-1]) / (cg_max[-1] - cg_max[-2])
regr1_max, cg1_max = np.array([lemac_lf[0]]), np.array([cg_max[0]])
regr2_max, cg2_max = np.array([]), np.array([])

for i in range(len(lemac_lf) - 1):

    if round(((lemac_lf[i] - lemac_lf[i + 1]) / (cg_max[i + 1] - cg_max[i])), 1) > 0.5:
        regr1_max = np.append(regr1_max, lemac_lf[i + 1])
        cg1_max = np.append(cg1_max, cg_max[i + 1])

    elif round(((lemac_lf[i] - lemac_lf[i + 1]) / (cg_max[i + 1] - cg_max[i])), 1) < 0.5:
        regr2_max = np.append(regr2_max, lemac_lf[i])
        cg2_max = np.append(cg2_max, cg_max[i])

regr2_max = np.append(regr2_max, lemac_lf[i + 1])
cg2_max = np.append(cg2_max, cg_max[i + 1])

lemac1_max_m, lemac1_max_b = regression(cg1_max, regr1_max)
lemac2_max_m, lemac2_max_b = regression(cg2_max, regr2_max)

# GET DESIGN CONDITION
ShS_opt = np.array([])
min_cg_opt = np.array([])
max_cg_opt = np.array([])

x_intersect = (stab_b - contr_b) / (contr_m - stab_m)
min_cg, max_cg = F_min_cg(lemac_des / lf), F_max_cg(lemac_des / lf)

ShS_contr = contr_m * min_cg + contr_b
ShS_stab = stab_m * max_cg + stab_b
print('ShS_contr = ', ShS_contr, '; ShS_stab = ', ShS_stab)

ShS_des = max(ShS_contr, ShS_stab)

if ShS_des == ShS_contr:
    print('\nCritical design condition is controllability')

elif ShS_des == ShS_stab:
    print('\nCritical design condition is stability')

print('Design condition found:\nShS = ', ShS_des, '\nx_LEMAC/l_f = ', lemac_des / lf, '\ncg_min/MAC = ', min_cg,
      '\ncg_max/MAC = ', max_cg)

if __name__ == "__main__":
    # PLOT SCISSORS PLOT
    lemac = np.linspace(lemac_lf[0], lemac_lf[-1], 301)

    fig, ax1 = plt.subplots()

    ax1.set_xlabel('${x_{CoG}}/{MAC}$', fontsize=20)
    ax1.set_ylabel('${S_h}/{S}$', color='tab:red', fontsize=20)
    ax1.plot(Xcg, stab_m * Xcg + stab_b, color='tab:red', label='Stability line', marker='8', markevery=70)
    ax1.plot(Xcg, stability(Xcg)[1], linestyle='--', color='tab:red', label='Neutral line', marker='x', markevery=70)
    ax1.plot(Xcg, contr_m * Xcg + contr_b, color='blue', label='Controllability line', marker='*', markevery=70)
    if len(ShS_opt) != 0:
        ax1.axhline(y=ShS_des, linestyle='--', color='k')
    # ax1.set_ylim(( ax1.get_ylim()[0]-ShS_des, ax1.get_ylim()[1]-ShS_des ))
    ax1.tick_params(axis='y', labelcolor='tab:red')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel('${x_{LEMAC}}/{l_f}$', color='tab:green', fontsize=20)  # we already handled the x-label with ax1
    ax2.plot(cg_min, lemac, color='tab:green', label='Min. CoG', marker='v', markevery=40)
    ax2.plot(cg_max, lemac, color='orange', label='Max. CoG', marker='.', markevery=40)
    ax2.tick_params(axis='y', labelcolor='tab:green')

    l1 = ax1.get_ylim()
    l2 = ax2.get_ylim()
    f = lambda x: l2[0] + (x - l1[0]) / (l1[1] - l1[0]) * (l2[1] - l2[0])
    ticks = f(ax1.get_yticks())
    ax2.yaxis.set_major_locator(ticker.FixedLocator(ticks))

    shs = ax1.get_yticks()
    lem = ax2.get_yticks()
    lem_m, lem_b = regression(shs, lem)

    ax2.axhline(y=lemac_des / lf, linestyle='--', color='k', label='Design condition')
    # ax2.set_ylim(( ax2.get_ylim()[0]-lemac_des, ax2.get_ylim()[1]-lemac_des ))
    deltaLEMAC = lemac_des / lf - (lem_m * ShS_des + lem_b)
    ax2.set_ylim(ax2.get_ylim()[0] + deltaLEMAC, ax2.get_ylim()[1] + deltaLEMAC)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    ax1.legend(fontsize=14, loc="lower left")
    ax2.legend(fontsize=14, loc="upper right")
    plt.show()
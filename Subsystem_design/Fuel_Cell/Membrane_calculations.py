'''
membrane calculations from Fuel Cell book
'''

import numpy as np
import sympy as sym
from scipy.integrate import *

#constants
E_thermo = 1.0  #[V]            thermodynamic voltage
j = 0.5         #[A?cm^2]       operating current density
T = 90 + 273    #[K]            operating temperature
x_H2 = 0.9      #[-]            hydrogen mole fraction
x_O2 = 0.19     #[-]            oxygen mole fraction
x_H2O = 0.1     #[-]            cathode water mole fraction
pC = 3          #[atm]          cathode pressure
pA = 3          #[atm]          anode pressure
alpha = 0.5     #[-]            transfer coefficient
j_0 = 0.0001    #[A/cm^2]       exchange current density
tM = 183        #[micro m]      electrolyte thickness
tA = 350        #[micro m]      anode thickness
tC = 350        #[micro m]      cathode thickness
rho_dry = 0.0021 #[kg/cm^3]     dry density of Nafion 117
M_m = 1         #[kg/mol]       Nafion equivalent weight
n_drag_SAT = 2.5 #[-]           electro-osmothic drag coefficient saturated Nafion
R = 8.314       #[J/mol K]      Gas constant
F = 96485       #[C/mol]        Faraday constant
p_operating = 3 #[bar]

# critical properties of gases
M_H2 = 2.016    #[g/mol]        molecular weight
M_air = 28.964
M_O2 = 31.999
M_H2O = 18.015
M_N2 = 28.013
Tc_H2 = 33.3    #[K]            temperatures
Tc_air = 132.4
Tc_O2 = 154.4
Tc_H2O = 647.3
Tc_N2 = 126.2
pc_H2 = 12.8    #[atm]          pressures
pc_air = 37
pc_O2 = 49.7
pc_H2O = 217.5
pc_N2 = 33.5
#constants a and b for non-polar gases and a combination of polar and non-polar:
a_non_polar = 0.0002745
b_non_polar = 1.823
a_polar = 0.000364
b_polar = 2.334

#General parameters
p_SAT = 10**(-2.1794 + 0.02953 * T -9.1837 * 10**(-5) * T**2 + 1.4454 * 10**(-7) * T**3)    #[bar]
x_H2O_SAT = p_SAT/p_operating
a_w = x_H2O/x_H2O_SAT
if 0 < a_w <= 1:
    lambda = 14 * a_w
elif 1 < a_w <= 3:
    lambda  = 10 + 4 * a_w
else:
    print("a_w should be between 0 and 3")

D_lambda = e**(2416*(1/303 - 1/T))*(2.562 - 0.33*lambda + 0.0264 * lambda**2 - 0.000671 * lambda**3) * 10**(-6) #for lambda > 4 [cm^2/s]

D_H2_H2O_eff = (0.4**(1.5) * a_polar*(T/(Tc_H2*Tc_H2O))**b_polar*(pc_H2*pc_H2O)**(1/3)*(Tc_H2*Tc_H2O)**(5/12)*
                (1/M_H2*1/M_H2O)**(1/2))/p_operating
D_O2_H2O_eff = (0.4**(1.5) * a_polar*(T/(Tc_O2*Tc_H2O))**b_polar*(pc_O2*pc_H2O)**(1/3)*(Tc_O2*Tc_H2O)**(5/12)*
                (1/M_O2*1/M_H2O)**(1/2))/p_operating
D_O2_N2_eff = (0.4**(1.5) * a_non_polar*(T/(Tc_O2*Tc_N2))**b_non_polar*(pc_O2*pc_N2)**(1/3)*(Tc_O2*Tc_N2)**(5/12)*
                (1/M_O2*1/M_N2)**(1/2))/p_operating

#cathodic overpotential:
eta_cathodic = (R*T)/(4*alpha*F) * np.log(j/(j_0*pC*(x_O2-((tC*10**(-6) * j*10**(4) *R*T)/4*F*pC*101325*D_O2_N2_eff ))) #note that np.log = ln!

#ohmic overpotential:
sym.init_printing()
const_C,alpha_star = sym.symbols('C,alpha_star')
f = sym.Eq(14*(pA/p_SAT)*(x_H2O-tA*10**(-6)*((alpha_star*j/j_0*R*T)/(2*F*pC*101325*D_H2_H2O_eff*j_0))), 11* alpha_star/n_drag_SAT + C* np.exp(0))
g = sym.Eq(10 + 4 * (pc_air/p_SAT)*(x_H2O + tC*10**(-6) *((1+alpha_star)*j/j_0*R*T)/(2*F*pC*101325*D_O2_H2O_eff*j_0)),11* alpha_star/n_drag_SAT + C* np.exp((j*M_m*n_drag_SAT)/(22*F*rho_dry*D_lambda)))
sym.solve([f,g],(const_C,alpha_star))


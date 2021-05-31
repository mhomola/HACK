# -*- coding: utf-8 -*-
"""
Created on Sat May 29 23:08:22 2021

@author: sarar
"""

from math import *
import numpy as np

# https://ntrs.nasa.gov/api/citations/20000065620/downloads/20000065620.pdf

R = 8.31446261815324  # [J / (K * mol)]

def a_par(T):
    if T >= 300 and T <= 1000:
        a1 = 3.9508691
        a2 = 0.10207987
        a3 = 0.13124466 * 10**(-4)
        a4 = -0.76649284 * 10**(-7)
        a5 = 0.34503763 * 10**(-10)
        a6 = -0.52093574 * 10**(5)
        a7 = 0.21980951 * 10**(2)  
        
    elif T > 1000 and T <= 5000:
        a1 =  36.440206
        a2 = 0.054614801
        a3 = -0.16091151 * 10**(-4)
        a4 = 0.21478497 * 10**(-8)
        a5 = -0.10131180 * 10**(-12)
        a6 = -0.63890109 * 10**(5)
        a7 = -0.15798973 * 10**(3)
    
    return a1, a2, a3, a4, a5

T = np.arange(300,4000 ,50)
Cp = np.array([])
to_save_list = list()

for t in T:
    a1, a2, a3, a4, a5 = a_par(t)
    #Cp = np.append( Cp, (a1 + a2*t + a3*t*t + a4*t**3 + a5*t**4)*R )  # [J / (K*mol)]
    Cp = (a1 + a2*t + a3*t*t + a4*t**3 + a5*t**4)*R   # [J / (K*mol)]
    to_save_list.append([t,Cp])

print(to_save_list)
np.savetxt("C12H26_cp.dat", to_save_list)

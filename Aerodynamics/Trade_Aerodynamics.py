import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import math as m
""" Concept 1: Cargo Hold"""
# Aerodynamic impact is 0


""" Concept 3: Longer Fuselage """

""" Concept 4: Flat Bottom """

""" Concept 5: Wing Pods """

""" Concept 6: Beluga """




#######INSERT HERE DESIRED PARAMETERS

lf = 37.57 #fuselage length in [m]
df = 4.14  #fuselage maximum diameter in [m]
new_df = 5 #modified maximum diameter in [m]

cd_var_graph,curr_cd_graph,new_cd_graph = graph_method(lf, df, new_df) #cd variation calculated with the garphical drag

l_cockpit = 5.04 #[m]
l_cabin = 29.53 - l_cockpit # [m]
l_tail = lf - 29.53 #[m]


ct = 1.64 # [m] obtained from top view drawing
cr = 6.07 # [m] pbtained from top view drawing
sweep = 25 # [deg] wing sweep obtained from Elsevier data base

exposed_mgc = mgc(ct,cr,sweep,df)


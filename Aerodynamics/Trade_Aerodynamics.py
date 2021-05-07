import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import math as m


def fus_wet_surface(l_cockpit, l_cabin, l_tail, df):
    """
    :param l_cockpit: length of the cockpit in [m]
    :param l_cabin: length of the cabin in [m]
    :param l_tail: length of the tail in [m]
    :param df: diameter length in [m]
    :return: Wet surface of the fuselage in [m^2]
    """

    S_w_fus = m.pi * df/4 * (1/3/ l_cockpit **2 * ((4*l_cockpit**2+df**2/4)-df**3/8) - df + 4 * l_cabin + 2 * m.sqrt(l_tail**2+ df**2/4))

    return S_w_fus


""" Concept 1: Cargo Hold"""
# Aerodynamic impact is 0


""" Concept 3: Longer Fuselage """

""" Concept 4: Flat Bottom """

""" Concept 5: Wing Pods """
pods_number = 2
pods_cdo = 0.025
pods_interf = 1.3

pods_L1 = 0.5
pods_main_length = 5
pods_L3 = 0.5

S_w_pods = fus_wet_surface(pods_L1, pods_main_length, pods_L3)
pods_contrib_cdo = pods_cdo*pods_interf*S_w_pods #todo: divide by wet surface of the plane and by Cdo of A320

""" Concept 6: Beluga """


def graph_method(lf, df, new_df):
    """

    :param lf: fuselage length in [m]
    :param df: initial maximum fuselage diameter in [m]
    :param new_df: new maximum fuselage in [m]
    :return:
    """
    #We recreate the Toreenbeek figure 3-12 curve for Cd frontal

    Cd_points = [0.09, 0.043, 0.038, 0.07, 0.095] #points on y-axis
    ratio_fus = [1.5, 2, 2.3, 8, 10]               #ponts on x-axis

    cd_graph = interpolate.interp1d(ratio_fus, Cd_points, kind = "cubic") #graph interpolation
    ratios = np.arange(ratio_fus[0], ratio_fus[-1], 0.01) #interval of ratios

    #set the parameters of the aircraft

    curr_ratio = lf/df

    for i in ratios:
        if abs(curr_ratio-i) < 0.01:
            curr_cd = cd_graph(curr_ratio)

    new_ratio = lf/new_df

    for i in ratios:
        if abs(new_ratio-i) < 0.01:
            new_cd = cd_graph(new_ratio)

    plt.plot(ratios, cd_graph(ratios), label = "Front Cd")
    plt.xlabel("lf/df[-]")
    plt.ylabel("cd[-]")
    plt.xlim(0, 10)
    plt.ylim(0, 0.1)
    plt.plot(curr_ratio, curr_cd, "red", marker = "o", label = "Standard design")
    plt.plot(new_ratio, new_cd, "green", marker = "o", label = "Modified Design")
    plt.legend()
    plt.show()

    #Compute the drag variation

    cd_var = (new_cd-curr_cd)/curr_cd*100 # [%]
    return cd_var, curr_cd, new_cd

#######INSERT HERE DESIRED PARAMETERS

lf = 37.57 #fuselage length in [m]
df = 4.14  #fuselage maximum diameter in [m]
new_df = 5 #modified maximum diameter in [m]

cd_var_graph, curr_cd_graph, new_cd_graph = graph_method(lf, df, new_df) #cd variation calculated with the garphical drag

l_cockpit = 5.04 #[m]
l_cabin = 29.53 - l_cockpit # [m]
l_tail = lf - 29.53 #[m]

S_w_fus = fus_wet_surface(l_cockpit, l_cabin, l_tail, df) #wetted surface area of the standard configuration


def mgc(ct,cr,sweep,df):
    """
    This function calculates the mean geometrical chord of the exposed wing.
    :param ct: tip chord in [m]
    :param cr: root chord in [m]
    :param sweep: sweep angle in [deg]
    :param df: fuselage diameter in [m]
    :return: mean geometrical chord in [m]
    """

    delta_cr = df/2/m.tan(m.radians(90-sweep)) #difference between root chord and exposed root chord
    exposed_cr = cr + delta_cr
    exposed_taper = ct/exposed_cr
    #Formula for mean geometric chord obtained from ADSEEII course

    exposed_mgc = 2/3 * exposed_cr * (1 + exposed_taper + exposed_taper**2)/(1+exposed_taper)

    return exposed_mgc

ct = 1.64 # [m] obtained from top view drawing
cr = 6.07 # [m] pbtained from top view drawing
sweep = 25 # [deg] wing sweep obtained from Elsevier data base

exposed_mgc = mgc(ct,cr,sweep,df)


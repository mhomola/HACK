import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import math as m

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
        if abs(curr_ratio-i)<0.01:
            curr_cd = cd_graph(curr_ratio)

    new_ratio = lf/new_df

    for i in ratios:
        if abs(new_ratio-i)<0.01:
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
    return cd_var,curr_cd,new_cd

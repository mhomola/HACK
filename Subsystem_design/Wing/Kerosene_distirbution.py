#wing geometry tank
import numpy as np
from Subsystem_design.common_constants import Constants
from scipy import interpolate
from matplotlib import pyplot as plt

def width_y(y,cr,taper,l):
    """
    Calculate the kerosene tank width at the span location y
    :param y: span location in [m]
    :param cr: tank root chord [m]
    :param l: maximum location of tanks along half_span [m]
    :return:
    """
    return cr - cr * (1-taper)*(y/l)

def height_y(y,h1,h2,l):
    """
    Calculate the kerosene tank height at the span location y
    :param h1: height at root
    :param h2: height at tip
    :param l: maximum location of tanks [m]
    :return:
    """
    taper = h2/h1
    return h1 - h1 * (1 - taper) * (y/l)

def kerosene_calc():
    """

    :return: kerosene_weight_function [function] defined on 0,l absolute value in [N]
            spanwise [array] of spanwise locations
            l [int] half span value at which the wing tanks end
    """
    c = Constants()

    l = 12.75   #graphically obtained
    t_c = 0.10  #maxium one for that airofil is 0.15 so this is assumed

    h1 = t_c * c.c_root
    h2 = t_c * c.c_tip

    V = 0
    dy = 0.01

    cr = 3.4            #play with this until you obtain volume from Airport manual
    taper = 0.240       #wing taper ratio elsevier

    spanwise = np.arange(0, l, dy)

    rho = 810   #[kg/m^3] density of kerosene
    weight = [] #list for storing spanwise weight

    for y in spanwise:
        dV = height_y(y,h1,h2,l) * width_y(y,cr,taper,l) * dy #rectangular cross section
        V = V + dV                                            #calculate total volume for verificatin purposes
        weight.append(dV*rho*9.81)

    kerosene_weight_function = interpolate.interp1d(spanwise,weight)
    return kerosene_weight_function,spanwise,l

if __name__ == '__main__':
    kerosene_weight_function,spanwise,l = kerosene_calc()
    plt.plot(spanwise,kerosene_weight_function(spanwise))
    plt.xlabel("Half_spanwise location [m]")
    plt.ylabel("Weight[N]")
    plt.show()
    print("chord at 0.55",width_y(y=0.55*14, cr=7.0465, taper=0.240, l=14))





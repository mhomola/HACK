import numpy as np
from Subsystem_design.common_constants import Constants


def pylon():
    c = Constants()
    # Material: AL-2090-T81
    shear_strength = 320 * 10**6    # [Pa]
    tensile_strength = 520 * 10**6  # [Pa]
    density_al_2090 = 2590          # [kg/m^3]

    # Max Loads:
    density_h2 = 71.1
    W = c.g_0 * (c.pod_tank_mass + c.pod_V_H2 * density_h2)
    W_max = W * 4.5    # Maximum Vertical Load due to 4.5g
    F_max = W * 9      # Maximum Front Load due to 9g
    S_max = W * 1.5    # Maximum Side Load due to 1.5g

    # Dimensions:
    l_pylon = 1.5                # pylon length [m]
    h_pylon = c.pylon_height     # pylon height [m]
    w_pylon = 0.2                # dummy width [m]

    """
    t1 thickness of the vertical sides of the cross section due to weight
    t2 thickness of the cross section due to shear, front load
    t3 thickness of the cross section due to shear, side load 
    """

    # ----------------------- Vertical Load: ---------------------
    A1 = W_max/tensile_strength
    t1 = A1/(2*l_pylon)     # thickness [mm]

    # ------------------- Shear around X due to front load of 9g: ----------------------------
    bad = True
    t2 = 0.001
    while bad:

        Ixx = w_pylon * h_pylon ** 3 / 12 - ((w_pylon - 2 * t2) * (h_pylon - 2 * t2) ** 3) / 12
        Qx = (h_pylon / 2 - t2) ** 2 * t2 + w_pylon * t2 * (h_pylon / 2 - t2 / 2)
        tau_xx = F_max * Qx / (Ixx * t2)
        tau_z = W_max * 1.09 * h_pylon / 2 / Ixx

        if tau_xx+tau_z < shear_strength:
            bad = False

        t2 = t2 + 0.000001
    t2 = t2 - 0.000001    # [m]

    # --------- Shear around Y due to side load of 1.5g: ----------------
    realbad = True
    t3 = 0.0001
    while realbad:

        Iyy = h_pylon * w_pylon ** 3 / 12 - ((h_pylon - 2 * t3) * (w_pylon - 2 * t3) ** 3) / 12
        Qy = (w_pylon / 2 - t3) ** 2 * t3 + h_pylon * t3 * (w_pylon / 2 - t3 / 2)
        tau_yy = S_max * Qy / (Iyy * t3)

        if tau_yy < shear_strength:
            realbad = False

        t3 = t3 + 0.000001
    t3 = t3 - 0.000001  # [m]

    t_list = [t1, t2, t3]
    t_pylon = np.max(t_list) * 1.5      # all computed thicknesses
    mass_pylon = l_pylon * (h_pylon * w_pylon - (h_pylon - 2 * t_pylon) * (w_pylon - 2 * t_pylon)) * density_al_2090

    return t_pylon, l_pylon, h_pylon, w_pylon, mass_pylon


t_pylon, l_pylon, h_pylon, w_pylon, mass_pylon = pylon()

if __name__ == '__main__':

    print(" ----- PYLON ----- ")
    print("Thickness of the cross section: ", t_pylon, " [m]")
    print("Pylon length: ", l_pylon, " [m]")
    print("Pylon height: ", h_pylon, " [m]")
    print("Pylon width: ", w_pylon, " [m]")
    print("Pylon mass: ", mass_pylon, " [kg]")





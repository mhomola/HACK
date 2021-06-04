"""
This files contains the classes and functions for determining the insulation properties of the tank.
"""
import numpy as np
Q_req = np.load("Q_required.npy") #fist element is the Q_req for the central tank, the second is for the aft tank

T_h = 273.15 + 45
T_c = 19.75

r_mid = 1.234659/2
r_pod = 2/2
l_mid = 3.14
l_pod = 5.74

P = 1.3*10**-3 * 0.0075006168         # 1.3 mPa -> Torr

# For DAM with Glass-Tissue

Cs = 7.3 * 10**-8
Cr = 7.07 * 10**-10
Cg = 1.46 * 10**-4
epsl = 0.021            # Effective Emissivity

N_density = 22.7        # layers/cm


def heat_flux_front():

    N_layers = 60  # number of layers
    q_tot = (Cr*epsl*(T_h**4.67-T_c**4.67)/N_layers) + (Cs*N_density**2.63*(T_h-T_c)*(T_h+T_c)/(2*N_layers+2)) + \
            (Cg*P*(T_h**0.52-T_c**0.52)/N_layers)

    Q_tot = (2*np.pi*r_mid*(l_mid-r_mid*2) + 4*np.pi*r_mid**2)*q_tot
    t_insul = N_layers / N_density

    return q_tot, Q_tot, t_insul


def heat_flux_aft():

    N_layers = 20 # number of layers
    q_tot = (Cr*epsl*(T_h**4.67-T_c**4.67)/N_layers) + (Cs*N_density**2.63*(T_h-T_c)*(T_h+T_c)/(2*N_layers+2)) + \
            (Cg*P*(T_h**0.52-T_c**0.52)/N_layers)

    Q_tot = (2*np.pi*r_pod*(l_pod-r_pod*2) + 4*np.pi*r_pod**2)*q_tot
    t_insul = N_layers / N_density

    return q_tot, Q_tot, t_insul


q_tank, Q_tank, t_insul_front = heat_flux_front()
q_tank_aft, Q_tank_aft, t_insul_aft = heat_flux_aft()

if __name__ == '__main__':

    print("Front Tanks: ")
    print("Required heat flow: ", Q_req[0])
    print("Designed heat flow: ", Q_tank)
    print("Thickness of Insulaton: ", t_insul_front)
    print("...........")
    print("Aft Tank: ")
    print("Required heat flow: ", Q_req[1])
    print("Designed heat flow: ", Q_tank_aft)
    print("Thickness of Insulaton: ", t_insul_aft)


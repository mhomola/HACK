from Subsystem_design.Wing.loads import Loads_w
from Subsystem_design.common_constants import Constants
import scipy.integrate as spint
from Subsystem_design.Wing.inertia import Inertia_normal, Inertia_shear
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

class Deflections(Constants):

    def __init__(self,):
        super(Deflections, self).__init__()
        self.E = 135 * 10**9
        self.G = 5.2 * 10**9
        self.dz = 0.1

    def vertical_deflection(self, z):

        z_arr = np.arange(self.width_f/2, z + self.dz/2, self.dz)

        Ixx_arr = np.zeros(len(z_arr))
        Mx_arr = np.zeros(len(z_arr))
        int1 = np.zeros(len(z_arr))

        I = Inertia_normal(n_str=11)
        L = Loads_w()
        L.compute_loads()
        L.chord_integrals()
        L.Reaction_forces()
        L.Reaction_moments()

        for i, z_i in enumerate(z_arr):
            I.compute_inertia(x=z_i)
            Ixx_arr[i] = I.Ixx_normal
            Mx_arr[i] = L.M_x(x=z_i)
            integrand = Mx_arr[:i+1]/Ixx_arr[:i+1]
            int1[i] = spint.simps(y=integrand, x=z_arr[:i+1])

        v_slope = -1/self.E * int1[-1]
        v = -1/self.E * spint.simps(y=int1, x=z_arr)

        return v, v_slope

    def horizontal_deflection(self, z):

        z_arr = np.arange(self.width_f/2, z + self.dz/2, self.dz)
        Iyy_arr = np.zeros(len(z_arr))
        My_arr = np.zeros(len(z_arr))
        int1 = np.zeros(len(z_arr))

        I = Inertia_normal(n_str=11)
        L = Loads_w()
        L.compute_loads()
        L.chord_integrals()
        L.Reaction_forces()
        L.Reaction_moments()

        for i, z_i in enumerate(z_arr):
            I.compute_inertia(x=z_i)
            Iyy_arr[i] = I.Iyy_normal
            My_arr[i] = L.M_y(x=z_i)
            integrand = My_arr[:i+1]/Iyy_arr[:i+1]
            int1[i] = spint.simps(y=integrand, x=z_arr[:i+1])

        w_slope = 1/self.E * int1[-1]
        w = 1/self.E * spint.simps(y=int1, x=z_arr)

        return w, w_slope

    def twist(self, z):

        z_arr = np.arange(self.width_f/2, z + self.dz/2, self.dz)
        J_arr = np.zeros(len(z_arr))
        T_arr = np.zeros(len(z_arr))

        L = Loads_w()
        L.compute_loads()
        L.chord_integrals()
        L.Reaction_forces()
        L.Reaction_moments()
        L.Reaction_torque()

        for i, z_i in enumerate(z_arr):

            T_arr[i] = L.T_z(x=z_i)

        integrand = T_arr/J_arr
        theta = 1/self.G * spint.simps(y=integrand, x=z_arr)

        return theta

    def plot_deflections(self):

        z_values = np.arange(self.width_f / 2, self.b / 2 + self.dz, self.dz*2)
        w_slope_arr = np.zeros(len(z_values))
        w_arr = np.zeros(len(z_values))
        v_slope_arr = np.zeros(len(z_values))
        v_arr = np.zeros(len(z_values))

        for j, z_j in enumerate(z_values):

            w_j, w_slope_j = self.horizontal_deflection(z=z_j)
            w_arr[j] = w_j
            w_slope_arr[j] = w_slope_j
            v_j, v_slope_j = self.vertical_deflection(z=z_j)
            v_arr[j] = v_j
            v_slope_arr[j] = v_slope_j

        L = Loads_w()
        L.compute_loads()
        L.chord_integrals()
        L.Reaction_forces()
        L.Reaction_moments()
        L.Reaction_torque()
        L.plot_loads()
        plt.show()

        fig = plt.figure(constrained_layout=True, figsize=(15, 15))
        gs = GridSpec(3, 2, figure=fig)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, :])
        ax4 = fig.add_subplot(gs[2, :])

        ax1.plot(z_values, v_slope_arr)
        ax1.set_ylabel(r"$v'(z)$ [$m/m$]", size=15)
        ax1.set_xlabel(r'$z$ [$m$]', size=15)
        ax1.set_title(r'Slope of the wingbox in the vertical direction', size=20)

        ax2.plot(z_values, w_slope_arr)
        ax2.set_ylabel(r"$w'(z)$ [$m/m$]", size=15)
        ax2.set_xlabel(r'$z$ [$m$]', size=15)
        ax2.set_title(r'Slope of the wingbox in the horizontal direction', size=20)

        ax3.plot(z_values, v_arr)
        ax3.set_ylabel(r"$v(z)$ [$m$]", size=15)
        ax3.set_title(r'Deflection of the wingbox in the vertical direction', size=20)

        ax4.plot(z_values, w_arr)
        ax4.set_ylabel(r"$w(z)$ [$m$]", size=15)
        ax4.set_xlabel(r'$z$ [$m$]', size=15)
        ax4.set_title(r'Deflection of the wingbox in the horiozntal direction', size=20)

        for i in [ax1, ax2, ax3, ax4]:
            i.tick_params(axis='both', which='major', labelsize=12)

        plt.show()

if __name__ == '__main__':
    delfections = Deflections()
    delfections.plot_deflections()

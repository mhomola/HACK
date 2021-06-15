from Subsystem_design.common_constants import Constants
import numpy as np
import matplotlib.pyplot as plt
import math as m

class Inertia_normal_var(Constants):

    def __init__(self,n_str,t_str,t_sp,t_sk,h_str,w_str,w_sk_c):
        super().__init__()
        self.t_str = t_str   # Stringer thickness [m]                                          #VARIABLE
        self.t_sp = t_sp     # Spar thickness [m]                                              #VARIABLE
        self.t_sk = t_sk  # Skin thickness [m]                                              #FIXED
        self.h_str = h_str    # Stringer height [m]                                             #VARIABLE
        self.w_str = w_str     # Stringer width [m]
        self.n_str = n_str    # Number of stringers on top and bottom (n_str * 2 = total_n_str)
        self.w_sk_c = w_sk_c    # Width of the skin over the local chord length                   #FIXED

    def chord_inertia(self, x):
        return (self.c_tip - self.c_kink_out) / (0.5 * self.b_out) * (x - 0.5 * self.b_in) + self.c_kink_out

    def t_c_distribution(self, x):
        y_b_arr = np.array([0.1138, 0.2, 0.2576, 0.2787, 0.306, 0.3276, 0.37, 0.5, 0.6, 0.7, 0.9, 1])
        t_c_arr = np.array([0.163, 0.1267, 0.117, 0.11546, 0.115, 0.1152, 0.117, 0.1123, 0.110, 0.109, 0.108, 0.108])

        yb = x / self.b * 2

        if yb >= .9:
            self.tc = 0.108

        elif yb <= 0.1138:
            self.tc = 0.163

        else:
            loc = np.where(y_b_arr >= yb)[0][:2] - 1
            self.tc = np.diff(t_c_arr[loc]) / np.diff(y_b_arr[loc]) * (yb - y_b_arr[loc][0]) + t_c_arr[loc][0]

        self.h_sp_c = 0.83 * self.tc

    def wb_config(self, x):

        c = self.chord_inertia(x=x)
        self.t_c_distribution(x=x)

        sep_str = self.w_sk_c * c / (self.n_str + 1)
        self.x_loc_str = np.zeros(self.n_str)
        y_loc_str = np.ones(self.n_str) * self.h_sp_c * c / 2

        for i in range(self.n_str):
            self.x_loc_str[i] = self.w_sk_c * c / 2 - sep_str * (i + 1)

        self.x_loc_str = np.hstack((self.x_loc_str, self.x_loc_str))
        y_loc_str = np.hstack((y_loc_str, -y_loc_str))

        sp_left_y = np.linspace(-self.h_sp_c * c / 2, self.h_sp_c * c / 2, 3)
        sp_left_x = np.ones(len(sp_left_y)) * self.w_sk_c * c / 2
        sp_right_y = sp_left_y
        sp_right_x = - sp_left_x

        sk_top_x = np.linspace(-self.w_sk_c * c / 2, self.w_sk_c * c / 2, 5)
        sk_top_y = np.ones(len(sk_top_x)) * (-self.h_sp_c * c / 2)
        sk_bot_x = sk_top_x
        sk_bot_y = - sk_top_y

        # plt.plot(sp_left_x, sp_left_y, color="black")
        # plt.plot(sp_right_x, sp_right_y, color="black")
        # plt.plot(sk_top_x, sk_top_y, color="black")
        # plt.plot(sk_bot_x, sk_bot_y, color="black")
        # plt.scatter(x=self.x_loc_str, y=y_loc_str, color='red', marker='o')
        # plt.scatter(x=0, y=0, color='g', marker='+')
        # plt.axis('equal')
        # plt.show()


    def compute_inertia(self, x):

        c = self.chord_inertia(x=x)
        self.wb_config(x=x)

        # Inertia around x axis

        Ixx_sp = 1/12 * self.t_sp * (self.h_sp_c * c)**3
        Ixx_sk = 1/12 * self.t_sk**3 * (self.w_sk_c * c) + (self.t_sk * (self.w_sk_c * c)) * (self.h_sp_c/2 * c)**2
        Ixx_str = 1/12 * self.t_str**3 * self.w_str + 1/12 * self.t_str * self.h_str**3 \
                  + self.t_str * self.w_str * (self.h_sp_c/2 * c)**2 \
                  + self.t_str * self.h_str * (self.h_sp_c/2 * c - self.h_str/2)**2

        self.Ixx_normal = 2 * Ixx_sp + 2 * Ixx_sk + 2 * self.n_str * Ixx_str

        # Inertia around y axis
        Iyy_sp = 1/12 * self.t_sp**3 * (self.h_sp_c * c) + self.t_sp * (self.h_sp_c * c) * (self.w_sk_c * c / 2)**2
        Iyy_sk = 1/12 * self.t_sk * (self.w_sk_c * c)**3
        Iyy_str = 2 * self.n_str * (1/12 * self.t_str**3 * self.h_str + 1/12 * self.t_str * self.w_str**3) \
                  + self.t_str * self.w_str * self.x_loc_str**2 + self.t_str * self.h_str * self.x_loc_str**2

        self.Iyy_normal = 2 * Iyy_sp + 2 * Iyy_sk + np.sum(Iyy_str)

class Inertia_shear_var(Constants):

    def __init__(self,n_str,t_sp,t_sk,w_sk_c):
        super().__init__()
        self.t_sp = t_sp  # Spar thickness [m]                                              #VARIABLE
        self.t_sk = t_sk  # Skin thickness [m]                                              #FIXED
        self.w_sk_c = w_sk_c   # Width of the skin over the local chord length                   #FIXED

    def chord_inertia(self, x):
        return (self.c_tip - self.c_kink_out) / (0.5 * self.b_out) * (x - 0.5 * self.b_in) + self.c_kink_out

    def t_c_distribution(self, x):
        y_b_arr = np.array([0.1138, 0.2, 0.2576, 0.2787, 0.306, 0.3276, 0.37, 0.5, 0.6, 0.7, 0.9, 1])
        t_c_arr = np.array([0.163, 0.1267, 0.117, 0.11546, 0.115, 0.1152, 0.117, 0.1123, 0.110, 0.109, 0.108, 0.108])

        yb = x / self.b * 2

        if yb >= .9:
            self.tc = 0.108

        elif yb <= 0.1138:
            self.tc = 0.163

        else:
            loc = np.where(y_b_arr >= yb)[0][:2] - 1
            self.tc = np.diff(t_c_arr[loc]) / np.diff(y_b_arr[loc]) * (yb - y_b_arr[loc][0]) + t_c_arr[loc][0]

        self.h_sp_c = 0.83 * self.tc

    def compute_inertia(self, x):

        c = self.chord_inertia(x=x)
        self.t_c_distribution(x=x)

        # Inertia around x axis

        Ixx_sp = 1/12 * self.t_sp * (self.h_sp_c * c)**3
        Ixx_sk = 1/12 * self.t_sk**3 * (self.w_sk_c * c) + (self.t_sk * (self.w_sk_c * c)) * (self.h_sp_c/2 * c)**2

        self.Ixx_shear = 2 * Ixx_sp + 2 * Ixx_sk

        # Inertia around y axis
        Iyy_sp = 1/12 * self.t_sp**3 * (self.h_sp_c * c) + self.t_sp * (self.h_sp_c * c) * (self.w_sk_c * c / 2)**2
        Iyy_sk = 1/12 * self.t_sk * (self.w_sk_c * c)**3

        self.Iyy_shear = 2 * Iyy_sp + 2 * Iyy_sk

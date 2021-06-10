from Subsystem_design.common_constants import Constants
import numpy as np
import matplotlib.pyplot as plt

class Inertia(Constants):

    def __init__(self):
        super().__init__()
        self.t_str = 0.02 # Stringer thickness [m]
        self.t_sp = 0.05  # Spar thickness [m]
        self.t_sk = 9.1/1000  # Skin thickness [m]
        self.h_str = 0.03  # Stringer height [m]
        self.h_sp_c = 0.091  # Height of the spar over local chord length
        self.n_str = 13  # Number of stringers on top and bottom (n_str * 2 = total_n_str)

        self.w_sk_c = 0.43  # Width of the skin over the local chord length

    def chord_inertia(self, x):
        return (self.c_tip - self.c_kink_out) / (0.5 * self.b_out) * (x - 0.5 * self.b_in) + self.c_kink_out

    def wb_config(self, x):

        c = self.chord_inertia(x=x)

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

        plt.plot(sp_left_x, sp_left_y, color="black")
        plt.plot(sp_right_x, sp_right_y, color="black")
        plt.plot(sk_top_x, sk_top_y, color="black")
        plt.plot(sk_bot_x, sk_bot_y, color="black")
        plt.scatter(x=self.x_loc_str, y=y_loc_str, color='red', marker='o')
        plt.scatter(x=0, y=0, color='g', marker='+')
        plt.axis('equal')
        plt.show()


    def compute_inertia(self, x):

        c = self.chord_inertia(x=x)
        self.wb_config(x=x)

        # Inertia around x axis

        Ixx_sp = 1/12 * self.t_sp * (self.h_sp_c * c)**3
        Ixx_sk = 1/12 * self.t_sk**3 * (self.w_sk_c * c) + (self.t_sk * (self.w_sk_c * c)) * (self.h_sp_c/2 * c)**2
        Ixx_str = 1/12 * self.t_str**3 * self.h_str + 1/12 * self.t_str * self.h_str**3 \
                  + (self.t_str * self.h_str) * ((self.h_sp_c/2 * c)**2 + (self.h_sp_c/2 * c - self.h_str/2)**2)

        self.Ixx_no_str = 2 * Ixx_sp + 2 * Ixx_sk
        self.Ixx = 2 * Ixx_sp + 2 * Ixx_sk + 2 * self.n_str * Ixx_str

        # Inertia around y axis
        Iyy_sp = 1/12 * self.t_sp**3 * (self.h_sp_c * c) + self.t_sp * (self.h_sp_c * c) * (self.w_sk_c * c / 2)**2
        Iyy_sk = 1/12 * self.t_sk * (self.w_sk_c * c)**3
        Iyy_str = 2 * self.n_str * (1/12 * self.t_str**3 * self.h_str + 1/12 * self.t_str * self.h_str**3) \
                  + (2 * self.t_str * self.h_str) * self.x_loc_str**2

        self.Iyy_no_str = 2 * Iyy_sp + 2 * Iyy_sk
        self.Iyy = 2 * Iyy_sp + 2 * Iyy_sk + np.sum(Iyy_str)

if __name__ == '__main__':
    ia = Inertia()
    ia.wb_config(x=0)
    ia.wb_config(x=ia.b/2)
    ia.compute_inertia(x=0)
    print(ia.Iyy)
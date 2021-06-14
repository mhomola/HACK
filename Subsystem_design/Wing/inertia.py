from Subsystem_design.common_constants import Constants
import numpy as np
import matplotlib.pyplot as plt
import math as m

class Inertia_initial(Constants):

    def __init__(self,n_str):
        super().__init__()
        self.t_str = 0.0025   # Stringer thickness [m]                                          #VARIABLE
        self.t_sp = 0.01      # Spar thickness [m]                                              #VARIABLE
        self.t_sk = 9.1/1000  # Skin thickness [m]                                              #FIXED
        self.h_str = 0.03     # Stringer height [m]                                             #VARIABLE
        self.w_str = 0.03     # Stringer width [m]
        self.n_str = n_str    # Number of stringers on top and bottom (n_str * 2 = total_n_str)
        self.w_sk_c = 0.43    # Width of the skin over the local chord length                   #FIXED

    def chord_inertia(self, x):
        return (self.c_tip - self.c_kink_out) / (0.5 * self.b_out) * (x - 0.5 * self.b_in) + self.c_kink_out

    def t_c_distribution(self, x):
        y_b_arr = np.array([0.1138, 0.2, 0.2576, 0.2787, 0.306, 0.3276, 0.37, 0.5, 0.6, 0.7, 0.9, 1])
        t_c_arr = np.array([0.163, 0.1267, 0.117, 0.11546, 0.115, 0.1152, 0.117, 0.1123, 0.110, 0.109, 0.108, 0.108])


        yb = x / self.b * 2
        loc = np.where(y_b_arr >= yb)[0][:2] - 1
        self.tc = np.diff(t_c_arr[loc]) / np.diff(y_b_arr[loc]) * (yb - y_b_arr[loc][0]) + t_c_arr[loc][0]

        if yb >= .9:
            self.tc = 0.108

        self.h_sp_c = 0.83 * self.tc

        # plt.plot(y_b_arr, t_c_arr)
        # plt.scatter([yb], self.tc)
        # plt.show()

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
    i = Inertia_initial(n_str=12)
    y_arr = np.linspace(i.width_f/2, i.b/2, 100)
    Ixx_arr = np.zeros(len(y_arr))
    for j, y in enumerate(y_arr):
        i.compute_inertia(x=y)
        Ixx_arr[j] = i.Ixx

    plt.plot(y_arr, Ixx_arr)
    plt.show()

class Inertia(Constants):

    def __init__(self,n_str):
        super().__init__()
        self.t_str = 0.0025 # Stringer thickness [m]                                          #VARIABLE
        self.t_sp = 0.01  # Spar thickness [m]                                               #VARIABLE
        self.t_sk = 9.1/1000  # Skin thickness [m]                                           #FIXED
        self.h_str = 0.03   # Stringer height [m]                                             #VARIABLE
        self.w_str = 0.03   # Stringer width [m]
        self.h_sp_c = 0.091  # Height of the spar over local chord length                    #FIXED
        self.n_str = n_str # Number of stringers on top and bottom (n_str * 2 = total_n_str)
        self.w_sk_c = 0.43  # Width of the skin over the local chord length                  #FIXED

    def chord_inertia(self, x):
        return (self.c_tip - self.c_kink_out) / (0.5 * self.b_out) * (x - 0.5 * self.b_in) + self.c_kink_out

    def compute_inertia(self,x):
        c = self.chord_inertia(x=x)
        height = c*self.h_sp_c #sizing the spar height to the corresponding chord
        width = c*self.w_sk_c  #sizing the skin width to the corresponding chord

        ###Calculate the skin MOI

        skin_Ixx = 2 * (1/12) * pow(height,3) * self.t_sp + \
                        2 * (1/12) * pow(self.t_sk,3) * width + 2 * width * self.t_sk * (height/2-self.t_sk/2)**2

        skin_Iyy = 2 * (1/12) * pow(self.t_sp,3) * height + 2 * height * self.t_sp * (width/2 - self.t_sp/2)**2 +\
                        2 * (1/12) * pow(width,3) * self.t_sk

        ###MOI around for stringers I-beam

        str_Ixx = 1/12 * self.t_str * pow(self.h_str,3) + 2 * 1/12 * pow(self.t_str,3) * self.w_str # MOI for just oen stringer
        str_Ixx_steiner = self.w_str * self.t_str * pow(height/2-self.h_str + self.t_str/2,2) + \
                          self.w_str * self.t_str * pow(height/2-self.t_str/2,2) +\
                          self.h_str * self.t_str * pow(height/2 - self.h_str/2,2)

        str_Ixx = self.n_str * str_Ixx
        str_Ixx_steiner  = self.n_str * str_Ixx_steiner

        ###MOI around y for I-beams
        ###Here the distribution of the stringers is important as they change the distancce from x

        if self.n_str == 0:
            str_Iyy = 0
            str_Iyy_steiner = 0

        if self.n_str == 1:
            str_Iyy = 1 / 12 * pow(self.t_str, 3) * self.h_str + 2 * 1 / 12 * self.t_str * pow(self.w_str, 3)
            str_Iyy_steiner = 0

        elif self.n_str%2 ==0:
            str_Iyy = 1/12 * pow(self.t_str,3) * self.h_str + 2 * 1/12 * self.t_str * pow(self.w_str,3)
            spacing = width/(self.n_str-1) #distance between the stringers
            str_Iyy_steiner = 0
            for i in range(1,m.ceil(self.n_str/2)+1):
                str_Iyy_steiner = str_Iyy_steiner + 2 * self.w_str * self.t_str * (spacing * i)**2 + self.h_str * self.t_str * (spacing * i)**2
            str_Iyy = str_Iyy * self.n_str
            str_Iyy_steiner = str_Iyy_steiner * 2 # we already looked at one side, we need to multiply by 2

        elif self.n_str%2 == 1: #odd number of stringers
            str_Iyy = 1/12 * pow(self.t_str,3) * self.h_str + 2 * 1/12 * self.t_str * pow(self.w_str,3)
            spacing = width/(self.n_str-1) #distance between the stringers
            str_Iyy_steiner = 0
            for i in range(1, m.ceil((self.n_str-1) / 2) + 1):
                str_Iyy_steiner = str_Iyy_steiner + 2 * self.w_str * self.t_str * (
                            spacing * i) ** 2 + self.h_str * self.t_str * (spacing * i) ** 2
            str_Iyy = str_Iyy * self.n_str
            str_Iyy_steiner = str_Iyy_steiner * 2  # we already looked at one side, we need to multiply by 2

        self.Ixx_no_str = skin_Ixx
        self.Ixx = self.Ixx_no_str + str_Ixx + str_Ixx_steiner
        self.Iyy_no_str = skin_Iyy
        self.Iyy = self.Iyy_no_str + str_Iyy + str_Iyy_steiner

# if __name__ == '__main__':
#     MOI_initial = Inertia_initial(n_str = 10)
#     MOI_ibeam = Inertia(n_str = 10)
#     x = 6
#     MOI_initial.compute_inertia(x=x)
#     MOI_ibeam.compute_inertia(x=x)
#     print("Ixx no stringers: " ,"FIRST BEAM" , MOI_initial.Ixx_no_str , "IBEAM", MOI_ibeam.Ixx_no_str)
#     print("Iyy no stringers: ", "FIRST BEAM", MOI_initial.Iyy_no_str, "IBEAM", MOI_ibeam.Iyy_no_str)
#     print("Ixx with stringers: ", "FIRST BEAM", MOI_initial.Ixx, "IBEAM", MOI_ibeam.Ixx)
#     print("Iyy with stringers: ", "FIRST BEAM", MOI_initial.Iyy, "IBEAM", MOI_ibeam.Iyy)
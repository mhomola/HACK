import unittest
# from Subsystem_design.Wing import loads
import numpy as np
import matplotlib.pyplot as plt
from Subsystem_design.Wing import Stresses
# class TestLoads(unittest.TestCase):
#     def test_lift(self):
#         """
#         In this test we look graphically at the lift prime distribution. It is expected to be the
#         highest at the root and the smallest at the tip.
#         """
#         lw = loads.Loads_w()
#         lw.compute_loads()
#         x_arr = np.linspace(0.5*lw.width_f, lw.b/2, 1000)
#         L_arr = np.zeros(1000)
#
#         for i, x in enumerate(x_arr):
#             L_arr[i] = lw.Lift(x)
#         plt.plot(x_arr,L_arr/10**3)
#         plt.xlabel("Span distance [m]")
#         plt.ylabel("Lift[kN]")
#         plt.show()
#
#     # def test_drag(self):
#     #     """
#     #     In this test we look graphically at the drag distribution. It is expected to be the
#     #     highest at the root and the smallest at the tip.
#     #     """
#     #     lw = loads.Loads_w()
#     #     lw.compute_loads()
#     #     x_arr = np.linspace(0.5 * lw.width_f, lw.b / 2, 1000)
#     #     D_arr = np.zeros(1000)
#     #
#     #     for i, x in enumerate(x_arr):
#     #         D_arr[i] = lw.Draq_w(x)  ###Why does not this work :((
#     #
#     #     plt.plot(x_arr, D_arr / 10 ** 3)
#     #     plt.xlabel("Span distance [m]")
#     #     plt.ylabel("Drag[kN]")
#     #     plt.show()
#
#     def test_body_weight(self):
#         lw = loads.Loads_w()
#         lw.compute_loads()
#         x_arr = np.linspace(0.5 * lw.width_f, lw.b / 2, 1000)
#         Wb_arr = np.zeros(1000)
#
#         for i, x in enumerate(x_arr):
#             Wb_arr[i] = lw.body_weight(x)
#
#         plt.plot(x_arr, Wb_arr / 10 ** 3)
#         plt.xlabel("Span distance [m]")
#         plt.ylabel("Body Weight[kN]")
#         plt.show()

class TestStresses(unittest.TestCase):
    def test_shearX(self):
        """
        If we only have a Vx load applied, this means that at y=0 our stresses should be 0.
        So q2(s2/2) = 0 and q4(s4) = 0 and q5(0) = 0
        :return:
        """
        h = 2
        L = 5
        t = 0.002
        Ixx = L * h ** 3 / 12 - (L - t) * (h - t) ** 3 / 12
        Iyy = h * L ** 3 / 12 - (h - t) * (L - t) ** 3 / 12

        stress = Stresses.stresses(Ixx=Ixx,Iyy=Iyy,Ixx_str = Ixx, Iyy_str=Iyy,h=h,L=L,t_upper=t,t_spar1=t,t_spar2=t,t_lower=t)
        stress.shear_loads(Vx = 500,Vy = 0, T = 0)
        stress.compute_stresses()

        self.assertAlmostEqual(stress.q2_tot(stress.s2/2),0)
        self.assertAlmostEqual(stress.q4_tot(stress.s4),0)
        self.assertAlmostEqual(stress.q5_tot(0), 0)

    def test_shearY(self):
        """
        If we only have a Vy load applied, this means that at x=0 our stresses should be 0.
        So q1(0) = 0 and q3(s3/2) = 0 and q6(s6) = 0
        :return:
        """
        h = 2
        L = 5
        t = 0.002
        Ixx = L * h ** 3 / 12 - (L - t) * (h - t) ** 3 / 12
        Iyy = h * L ** 3 / 12 - (h - t) * (L - t) ** 3 / 12

        stress = Stresses.stresses(Ixx=Ixx,Iyy=Iyy,Ixx_str = Ixx, Iyy_str=Iyy,h=h,L=L,t_upper=t,t_spar1=t,t_spar2=t,t_lower=t)
        stress.shear_loads(Vx = 0,Vy = 500, T = 0)
        stress.compute_stresses()

        self.assertAlmostEqual(stress.q1_tot(0),0)
        self.assertAlmostEqual(stress.q3_tot(stress.s3/2),0)
        self.assertAlmostEqual(stress.q6_tot(stress.s6),0)

    def test_bending_Mx(self):
        """
        If we only have a Mx load applied, this means that the y + is in tension :
                    - sigma2(s2/2->s2) >0
                    - sigma3 >0
                    - sigma4>0
        and y- is in compression.
                    -sigma5<0
                    -sigma6<0
                    -sigma1<1
                    -sigma2(0->s2/2)
        :return:
        """

        h = 2
        L = 5
        t = 0.002
        Ixx = L * h ** 3 / 12 - (L - t) * (h - t) ** 3 / 12
        Iyy = h * L ** 3 / 12 - (h - t) * (L - t) ** 3 / 12

        stress = Stresses.stresses(Ixx=Ixx, Iyy=Iyy, Ixx_str=Ixx, Iyy_str=Iyy, h=h, L=L, t_upper=t, t_spar1=t,
                                   t_spar2=t, t_lower=t)
        stress.shear_loads(Vx=0, Vy=0, T=0)
        stress.bending_loads(Mx = 1000,My = 0)
        stress.compute_stresses()

        ###TENSION TESTS
        for i in np.arange(stress.s2/2*1.01,stress.s2,0.1):
            self.assertGreater(stress.sigma2(i),0)

        for i in np.arange(0,stress.s3,0.1):
            self.assertGreater(stress.sigma3(i),0)

        for i in np.arange(0,stress.s4,0.1):
            self.assertGreater(stress.sigma4(i),0)

        ###COMPRESSION TESTS
        for i in np.arange(0, stress.s5, 0.1):
            self.assertLess(stress.sigma5(i), 0)

        for i in np.arange(0, stress.s6, 0.1):
            self.assertLess(stress.sigma6(i), 0)

        for i in np.arange(0, stress.s1, 0.1):
            self.assertLess(stress.sigma1(i), 0)

        for i in np.arange(0, stress.s2/2*0.99, 0.1):
            self.assertLess(stress.sigma2(i), 0)

    def test_bending_My(self):
        """
        If we only have a Mx load applied, this means that the x + is in compression :
                    - sigma1<0
                    - sigma2<0
                    - sigma3(0->s3/2)<0
        and x- is in tension
                    -sigma3(s3/2->0)>0
                    -sigma4>0
                    -sigma5>0
                    -sigma6>0
        :return:
        """

        h = 2
        L = 5
        t = 0.002
        Ixx = L * h ** 3 / 12 - (L - t) * (h - t) ** 3 / 12
        Iyy = h * L ** 3 / 12 - (h - t) * (L - t) ** 3 / 12

        stress = Stresses.stresses(Ixx=Ixx, Iyy=Iyy, Ixx_str=Ixx, Iyy_str=Iyy, h=h, L=L, t_upper=t, t_spar1=t,
                                   t_spar2=t, t_lower=t)
        stress.shear_loads(Vx=0, Vy=0, T=0)
        stress.bending_loads(Mx = 0,My = 500)
        stress.compute_stresses()

        ###COMPRESSION TESTS

        for i in np.arange(0,stress.s1,0.1):
            self.assertLess(stress.sigma1(i),0)

        for i in np.arange(0,stress.s2,0.1):
            self.assertLess(stress.sigma2(i),0)

        for i in np.arange(0,stress.s3/2*0.99,0.1):
            self.assertLess(stress.sigma3(i),0)

        ###TENSION TESTS
        for i in np.arange(stress.s3/2*1.01, stress.s3, 0.1):
            self.assertGreater(stress.sigma3(i), 0)

        # for i in np.arange(0, stress.s4, 0.1):
        #     self.assertGreater(stress.sigma4(i), 0)
        #
        # for i in np.arange(0, stress.s5, 0.1):
        #     self.assertGreater(stress.sigma5(i), 0)
        #
        # for i in np.arange(0, stress.s6, 0.1):
        #     self.assertGreater(stress.sigma6(i), 0)








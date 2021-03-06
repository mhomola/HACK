import numpy as np
import scipy as sp
import scipy.integrate as spint
import math as m
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.collections as mcoll
import matplotlib.path as mpath

class stresses():
    def __init__(self,Ixx,Iyy,Ixx_str,Iyy_str,h,L,t_upper,t_spar1,t_spar2,t_lower):
        """
        :param Ixx: [m^4] MOI about x-axis(normal) no stringer
        :param Iyy: [m^4] MOI about y-axis(normal) no stringer
        :param Ixx_str: [m^4] MOI about x-axis(normal) with stringer
        :param Iyy_str: [m^4] MOI about y-axis(normal) with stringer
        :param h: [m] height of the wing box
        :param L: [m] width of the wingbox
        :param t_upper: [m] upper plate thickness
        :param t_spar1: [m] front spar thickness
        :param t_spar2: [m] rear spar thickness
        :param t_lower: [m] lower plate thickness
        """
        self.Ixx = Ixx
        self.Iyy = Iyy
        self.Ixx_str = Ixx_str
        self.Iyy_str = Iyy_str
        self.h = h
        self.L = L
        self.t_upper = t_upper
        self.t_spar1 = t_spar1
        self.t_spar2 = t_spar2
        self.t_lower = t_lower

        #maximum values of all the defined segments
        self.s1 = self.L/2
        self.s2 = self.h
        self.s3 = self.L
        self.s4 = self.h/2
        self.s5 = self.h/2
        self.s6 = self.L/2

    ### SHEAR FLOWS

    def shear_loads(self,Vx,Vy,T):
        """
        :param Vx: [N] horizontal shear force positive left
        :param Vy: [N] vertical shear force positive down
        :param T [N*m] torque shear flow positive is ccw
        :return:
        """
        self.Vx = Vx
        self.Vy = Vy
        self.T = T

    def shear_flows_y(self):
        """
        We make the cut along the symmetry line of y so we have qs0 = 0. The basic shear flows are the actual shear flows.
        :return:
        """

        def qb1(s): # y=ct
            return -self.Vy/self.Ixx * self.t_upper * (-self.h/2) * s

        def qb2(s): #y = -h/2+s
            return -self.Vy/self.Ixx * self.t_spar1 * (-self.h/2 * s + s**2/2) + qb1(self.s1)

        def qb3(s): #y = ct
            return -self.Vy/self.Ixx * self.t_lower *  (self.h/2) * s + qb2(self.s2)

        def qb4(s): # y from h/2 to 0
            return -self.Vy/self.Ixx * self.t_spar2 * (self.h/2 * s - s**2/2) + qb3(self.s3)

        def qb5(s): #y from o to -h/2
            return -self.Vy/self.Ixx * self.t_spar2 * (-s**2/2) + qb4(self.s4)

        def qb6(s): #y stays constant at -h/2
            return -self.Vy/self.Ixx * self.t_upper * (-self.h/2 * s) + qb5(self.s5)

        self.q1_Vy = qb1
        self.q2_Vy = qb2
        self.q3_Vy = qb3
        self.q4_Vy = qb4
        self.q5_Vy = qb5
        self.q6_Vy = qb6

    def shear_flows_x(self):
        """
        We make the cut along the symmetry line of x so we have qs0 = 0. The basic shear flows are the actual shear flows.
        Now we have a different order starting from the cur, that is why we start our analysis from s5.
        """

        def qb5(s):
            return -self.Vx/self.Iyy * self.t_spar2 * (-self.L/2*s)

        def qb6(s):
            return - self.Vx/self.Iyy * self.t_upper * (-self.L/2*s + s**2/2) + qb5(self.s5)

        def qb1(s):
            return - self.Vx/self.Iyy * self.t_upper * (s**2/2) + qb6(self.s6)

        def qb2(s):
            return -self.Vx/self.Iyy * self.t_spar1 * (self.L/2*s) + qb1(self.s1)

        def qb3(s):
            return -self.Vx/self.Iyy * self.t_lower * (self.L/2*s-s**2/2) + qb2(self.s2)

        def qb4(s):
            return -self.Vx/self.Iyy * self.t_spar2 * (-self.L/2*s) + qb3(self.s3)

        self.q1_Vx = qb1
        self.q2_Vx = qb2
        self.q3_Vx = qb3
        self.q4_Vx = qb4
        self.q5_Vx = qb5
        self.q6_Vx = qb6

    def shear_flows_torque(self):
        """
        q = T/(2*Am)
        :return:
        """
        def q_torque(s):
            return self.T/(2*self.L*self.h)

        self.q_torque = q_torque

    def shear_flows_total(self):

        def q1_tot(s):

            return self.q1_Vy(s) + self.q1_Vx(s) + self.q_torque(s)

        def q2_tot(s):
            return self.q2_Vy(s) + self.q2_Vx(s) + self.q_torque(s)

        def q3_tot(s):
            return self.q3_Vy(s) + self.q3_Vx(s) + self.q_torque(s)

        def q4_tot(s):
            return self.q4_Vy(s) + self.q4_Vx(s) + self.q_torque(s)

        def q5_tot(s):
            return self.q5_Vy(s) + self.q5_Vx(s) + self.q_torque(s)

        def q6_tot(s):
            return self.q6_Vy(s) + self.q6_Vx(s) + self.q_torque(s)

        self.q1_tot = q1_tot
        self.q2_tot = q2_tot
        self.q3_tot = q3_tot
        self.q4_tot = q4_tot
        self.q5_tot = q5_tot
        self.q6_tot = q6_tot

    def shear_flow_plotter(self,type,show):

        """
        Define for which type of shear you want the plotting
        - y for Vy
        - x for Vx
        -torque
        -total
        """

        if type == "Vy":
            q1_plot = self.q1_Vy
            q2_plot = self.q2_Vy
            q3_plot = self.q3_Vy
            q4_plot = self.q4_Vy
            q5_plot = self.q5_Vy
            q6_plot = self.q6_Vy

        elif type == "Vx":
            q1_plot = self.q1_Vx
            q2_plot = self.q2_Vx
            q3_plot = self.q3_Vx
            q4_plot = self.q4_Vx
            q5_plot = self.q5_Vx
            q6_plot = self.q6_Vx

        elif type == "torque":
            q1_plot = self.q_torque
            q2_plot = self.q_torque
            q3_plot = self.q_torque
            q4_plot = self.q_torque
            q5_plot = self.q_torque
            q6_plot = self.q_torque

        elif type == "total":
            q1_plot = self.q1_tot
            q2_plot = self.q2_tot
            q3_plot = self.q3_tot
            q4_plot = self.q4_tot
            q5_plot = self.q5_tot
            q6_plot = self.q6_tot

        n = 1000

        ###Region 1
        s1 = np.linspace(0,self.s1,num=n)
        x1 = np.linspace(0,self.L/2,num=n)
        y1 = -self.h/2 * np.ones(n)
        q1 = q1_plot(s1)
        shear_stress1 = q1/self.t_upper
        path = mpath.Path(np.column_stack([x1, y1]))
        verts = path.interpolated(steps=1).vertices
        x1, y1 = verts[:, 0, ], verts[:, 1]
        maxabs = np.max(np.abs(q1))

        ###Region 2
        s2 = np.linspace(0, self.s2, num=n)
        x2 = self.L / 2 * np.ones(n)
        y2 = np.linspace(-self.h/2, self.h/2, num=n)
        q2 = q2_plot(s2)
        shear_stress2 = q2 / self.t_spar1
        path = mpath.Path(np.column_stack([x2, y2]))
        verts = path.interpolated(steps=1).vertices
        x2, y2 = verts[:, 0, ], verts[:, 1]
        maxabs2 = np.max(np.abs(q2))
        maxabs = max(maxabs2, maxabs)

        ###Region 3
        s3 = np.linspace(0, self.s3, num=n)
        x3 = np.linspace(self.L/2, -self.L/2, num=n)
        y3 = self.h / 2 * np.ones(n)
        q3 = q3_plot(s3)
        shear_stress3 = q3/ self.t_lower
        path = mpath.Path(np.column_stack([x3, y3]))
        verts = path.interpolated(steps=1).vertices
        x3, y3 = verts[:, 0, ], verts[:, 1]
        maxabs3 = np.max(np.abs(q3))
        maxabs = max(maxabs3, maxabs)

        ###Region 4
        s4 = np.linspace(0, self.s4, num=n)
        x4 = -self.L/2 * np.ones(n)
        y4 = np.linspace(self.h/2, 0, num=n)
        q4 = q4_plot(s4)
        shear_stress4 = q4/ self.t_spar2
        path = mpath.Path(np.column_stack([x4, y4]))
        verts = path.interpolated(steps=1).vertices
        x4, y4 = verts[:, 0, ], verts[:, 1]
        maxabs4 = np.max(np.abs(q4))
        maxabs = max(maxabs4, maxabs)

        ###Region 5
        s5 = np.linspace(0, self.s5, num=n)
        x5 = -self.L / 2 * np.ones(n)
        y5 = np.linspace(0, -self.h/2, num=n)
        q5 = q5_plot(s5)
        shear_stress5 = q5 / self.t_spar2
        path = mpath.Path(np.column_stack([x5, y5]))
        verts = path.interpolated(steps=1).vertices
        x5, y5 = verts[:, 0, ], verts[:, 1]
        maxabs5 = np.max(np.abs(q5))
        maxabs = max(maxabs5, maxabs)

        ###Region 6
        s6 = np.linspace(0, self.s6, num=n)
        x6 = np.linspace(-self.L/2, 0, num=n)
        y6 = -self.h / 2 * np.ones(n)
        q6 = q6_plot(s6)
        shear_stress6 = q6 / self.t_upper
        path = mpath.Path(np.column_stack([x6, y6]))
        verts = path.interpolated(steps=1).vertices
        x6, y6 = verts[:, 0, ], verts[:, 1]
        maxabs6 = np.max(np.abs(q6))
        maxabs = max(maxabs6, maxabs)

        if show == True:
            fig = plt.figure(4)

            #the negative/positive comnvention only shows if the initial coord syst is respected
            q1 = np.absolute(q1)
            q2 = np.absolute(q2)
            q3 = np.absolute(q3)
            q4 = np.absolute(q4)
            q5 = np.absolute(q5)
            q6 = np.absolute(q6)

            colorline(x1, y1, q1, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=2)
            colorline(x2, y2, q2, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=2)
            colorline(x3, y3, q3, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=2)
            colorline(x4, y4, q4, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=2)
            colorline(x5, y5, q5, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=2)
            colorline(x6, y6, q6, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=2)

            sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'),
                                       norm=plt.Normalize(-maxabs, maxabs))
            sm.set_array([])

            plt.colorbar(sm, label=r'$q$ [N/m]', fraction=0.20, pad=0.04, orientation="horizontal")
            plt.xlim(-self.L - 0.1, 0.1)
            plt.ylim(self.h / 2 - 0.02, -self.h / 2 + 0.02)
            plt.axis('scaled')
            plt.gca().invert_xaxis()
            plt.xlabel(r'$z$ [m]')
            plt.ylabel(r'$y$ [m]')
            plt.title('Shear flow distribution')
            plt.show()
        self.q_max = maxabs
        self.shear_stress_max = max([max(shear_stress1), max(shear_stress2), max(shear_stress3), max(shear_stress4), max(shear_stress5), max(shear_stress6)])
    ### NORMAL STRESSES

    def bending_loads(self,Mx,My):
        """
        :param Mx: [N*m] positive defined with right hand rule
        :param My: [N*m] positive defined with right hand rule
        :return:
        """
        #Here we adjust the signs using the color of the first quadrant rule
        self.Mx = Mx
        self.My = -My

    def sigma_total(self):

        def sigma1(s):
            x = s
            y = -self.h/2
            return self.Mx * y/self.Ixx_str + self.My * x/ self.Iyy_str

        def sigma2(s):
            x = -self.L/2
            y = -self.h/2 + s
            return self.Mx * y/self.Ixx_str  + self.My * x/ self.Iyy_str

        def sigma3(s):
            x = self.L/2 - s
            y = self.h/2
            return self.Mx * y/self.Ixx_str  + self.My * x/ self.Iyy_str

        def sigma4(s):
            x = -self.L/2
            y = self.h/2-s
            return self.Mx * y / self.Ixx_str  + self.My * x / self.Iyy_str

        def sigma5(s):
            x = -self.L / 2
            y = - s
            return self.Mx * y / self.Ixx_str  + self.My * x / self.Iyy_str

        def sigma6(s):
            x = - self.L/2 + s
            y = - self.h/2
            return self.Mx * y / self.Ixx_str  + self.My * x / self.Iyy_str

        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.sigma3 = sigma3
        self.sigma4 = sigma4
        self.sigma5 = sigma5
        self.sigma6 = sigma6

    def sigma_plotter(self,show):

        n = 1000

        ###Region 1
        s1 = np.linspace(0, self.s1, num=n)
        x1 = np.linspace(0,self.L/2,num=n)
        y1 = -self.h/2 * np.ones(n)
        sigma1 = self.sigma1(s1)

        path = mpath.Path(np.column_stack([x1, y1]))
        verts = path.interpolated(steps=1).vertices
        x1, y1 = verts[:, 0, ], verts[:, 1]
        maxabs = np.max(np.abs(sigma1))

        ###Region 2
        s2 = np.linspace(0, self.s2, num=n)
        x2 = self.L / 2 * np.ones(n)
        y2 = np.linspace(-self.h/2, self.h/2, num=n)
        sigma2 = self.sigma2(s2)

        path = mpath.Path(np.column_stack([x2, y2]))
        verts = path.interpolated(steps=1).vertices
        x2, y2 = verts[:, 0, ], verts[:, 1]
        maxabs2 = np.max(np.abs(sigma2))
        maxabs = max(maxabs2, maxabs)

        ###Region 3
        s3 = np.linspace(0, self.s3, num=n)
        x3 = np.linspace(self.L/2, -self.L/2, num=n)
        y3 = self.h / 2 * np.ones(n)
        sigma3 = self.sigma3(s3)

        path = mpath.Path(np.column_stack([x3, y3]))
        verts = path.interpolated(steps=1).vertices
        x3, y3 = verts[:, 0, ], verts[:, 1]
        maxabs3 = np.max(np.abs(sigma3))
        maxabs = max(maxabs3, maxabs)

        ###Region 4
        s4 = np.linspace(0, self.s4, num=n)
        x4 = -self.L/2 * np.ones(n)
        y4 = np.linspace(self.h/2, 0, num=n)
        sigma4 = self.sigma4(s4)

        path = mpath.Path(np.column_stack([x4, y4]))
        verts = path.interpolated(steps=1).vertices
        x4, y4 = verts[:, 0, ], verts[:, 1]
        maxabs4 = np.max(np.abs(sigma4))
        maxabs = max(maxabs4, maxabs)

        ###Region 5
        s5 = np.linspace(0, self.s5, num=n)
        x5 = -self.L / 2 * np.ones(n)
        y5 = np.linspace(0, -self.h/2, num=n)
        sigma5 = self.sigma5(s5)
        path = mpath.Path(np.column_stack([x5, y5]))

        verts = path.interpolated(steps=1).vertices
        x5, y5 = verts[:, 0, ], verts[:, 1]
        maxabs5 = np.max(np.abs(sigma5))
        maxabs = max(maxabs5, maxabs)

        ###Region 6
        s6 = np.linspace(0, self.s6, num=n)
        x6 = np.linspace(-self.L/2, 0, num=n)
        y6 = -self.h / 2 * np.ones(n)
        sigma6 = self.sigma6(s6)
        path = mpath.Path(np.column_stack([x6, y6]))

        verts = path.interpolated(steps=1).vertices
        x6, y6 = verts[:, 0, ], verts[:, 1]
        maxabs6 = np.max(np.abs(sigma6))
        maxabs = max(maxabs6, maxabs)

        if show == True:
            fig = plt.figure(4)

            l_width = 8
            colorline(x1, y1, sigma1, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x2, y2, sigma2, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x3, y3, sigma3, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x4, y4, sigma4, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x5, y5, sigma5, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x6, y6, sigma6, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )

            sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'),
                                       norm=plt.Normalize(-maxabs, maxabs))
            sm.set_array([])

            plt.colorbar(sm, label=r'$Sigma$ [Pa]', fraction=0.20, pad=0.04, orientation="horizontal")
            plt.xlim(-self.L - 0.1, 0.1)
            plt.ylim(self.h / 2 - 0.02, -self.h / 2 + 0.02)
            plt.axis('scaled')
            plt.gca().invert_xaxis()
            plt.xlabel(r'$z$ [m]')
            plt.ylabel(r'$y$ [m]')
            plt.title('Normal stress distribution')
            plt.show()

        self.sigma_max = maxabs
        self.sigma_compressive = min([min(sigma1),min(sigma2),min(sigma3),min(sigma4),min(sigma5),min(sigma6)])
        self.sigma_tension = max([max(sigma1), max(sigma2), max(sigma3), max(sigma4), max(sigma5), max(sigma6)])


    def von_Misses(self):
        """
        Vm = square root ( nromal stress^2 + 3 * shear stress)
        shear stress = shear flow/thickness
        :return:
        """

        def vm1(s):
            return np.sqrt(self.sigma1(s)**2 + 3 * pow(self.q1_tot(s)/self.t_upper,2))

        def vm2(s):
            return np.sqrt(self.sigma2(s)**2 + 3 * pow(self.q2_tot(s)/self.t_spar1,2))

        def vm3(s):
            return np.sqrt(self.sigma3(s)**2 + 3 * pow(self.q3_tot(s)/self.t_lower,2))

        def vm4(s):
            return np.sqrt(self.sigma4(s) ** 2 + 3 * pow(self.q4_tot(s) / self.t_spar2,2))

        def vm5(s):
            return np.sqrt(self.sigma5(s) ** 2 + 3 * pow(self.q5_tot(s) / self.t_spar2,2))

        def vm6(s):
            return np.sqrt(self.sigma6(s) ** 2 + 3 * pow(self.q6_tot(s) / self.t_upper,2))

        self.vm1 = vm1
        self.vm2 = vm2
        self.vm3 = vm3
        self.vm4 = vm4
        self.vm5 = vm5
        self.vm6 = vm6

    def vm_plotter(self,show):

        n = 1000

        ###Region 1
        s1 = np.linspace(0, self.s1, num=n)
        x1 = np.linspace(0,self.L/2,num=n)
        y1 = -self.h/2 * np.ones(n)
        vm1 = self.vm1(s1)

        path = mpath.Path(np.column_stack([x1, y1]))
        verts = path.interpolated(steps=1).vertices
        x1, y1 = verts[:, 0, ], verts[:, 1]
        maxabs = np.max(np.abs(vm1))

        ###Region 2
        s2 = np.linspace(0, self.s2, num=n)
        x2 = self.L / 2 * np.ones(n)
        y2 = np.linspace(-self.h/2, self.h/2, num=n)
        vm2 = self.vm2(s2)

        path = mpath.Path(np.column_stack([x2, y2]))
        verts = path.interpolated(steps=1).vertices
        x2, y2 = verts[:, 0, ], verts[:, 1]
        maxabs2 = np.max(np.abs(vm2))
        maxabs = max(maxabs2, maxabs)

        ###Region 3
        s3 = np.linspace(0, self.s3, num=n)
        x3 = np.linspace(self.L/2, -self.L/2, num=n)
        y3 = self.h / 2 * np.ones(n)
        vm3 = self.vm3(s3)

        path = mpath.Path(np.column_stack([x3, y3]))
        verts = path.interpolated(steps=1).vertices
        x3, y3 = verts[:, 0, ], verts[:, 1]
        maxabs3 = np.max(np.abs(vm3))
        maxabs = max(maxabs3, maxabs)

        ###Region 4
        s4 = np.linspace(0, self.s4, num=n)
        x4 = -self.L/2 * np.ones(n)
        y4 = np.linspace(self.h/2, 0, num=n)
        vm4 = self.vm4(s4)

        path = mpath.Path(np.column_stack([x4, y4]))
        verts = path.interpolated(steps=1).vertices
        x4, y4 = verts[:, 0, ], verts[:, 1]
        maxabs4 = np.max(np.abs(vm4))
        maxabs = max(maxabs4, maxabs)

        ###Region 5
        s5 = np.linspace(0, self.s5, num=n)
        x5 = -self.L / 2 * np.ones(n)
        y5 = np.linspace(0, -self.h/2, num=n)
        vm5 = self.vm5(s5)
        path = mpath.Path(np.column_stack([x5, y5]))

        verts = path.interpolated(steps=1).vertices
        x5, y5 = verts[:, 0, ], verts[:, 1]
        maxabs5 = np.max(np.abs(vm5))
        maxabs = max(maxabs5, maxabs)

        ###Region 6
        s6 = np.linspace(0, self.s6, num=n)
        x6 = np.linspace(-self.L/2, 0, num=n)
        y6 = -self.h / 2 * np.ones(n)
        vm6 = self.vm6(s6)
        path = mpath.Path(np.column_stack([x6, y6]))

        verts = path.interpolated(steps=1).vertices
        x6, y6 = verts[:, 0, ], verts[:, 1]
        maxabs6 = np.max(np.abs(vm6))
        maxabs = max(maxabs6, maxabs)

        if show == True:
            fig = plt.figure(4)

            l_width = 8
            colorline(x1, y1, vm1, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x2, y2, vm2, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x3, y3, vm3, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x4, y4, vm4, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x5, y5, vm5, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )
            colorline(x6, y6, vm6, cmap=plt.get_cmap('jet'),
                      norm=plt.Normalize(-maxabs, maxabs), linewidth=l_width )

            sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'),
                                       norm=plt.Normalize(-maxabs, maxabs))
            sm.set_array([])

            plt.colorbar(sm, label=r'$Von Mises$ [Pa]', fraction=0.20, pad=0.04, orientation="horizontal")
            plt.xlim(-self.L - 0.1, 0.1)
            plt.ylim(self.h / 2 - 0.02, -self.h / 2 + 0.02)
            plt.axis('scaled')
            plt.gca().invert_xaxis()
            plt.xlabel(r'$z$ [m]')
            plt.ylabel(r'$y$ [m]')
            plt.title('Von Misses stress distribution')
            plt.show()

        self.vm_max = maxabs
        self.vm1 = vm1
        self.vm2 = vm2
        self.vm3 = vm3
        self.vm4 = vm4
        self.vm5 = vm5
        self.vm6 = vm6

    def compute_stresses(self):
        self.shear_flows_y()
        self.shear_flows_x()
        self.shear_flows_torque()
        self.shear_flows_total()
        self.sigma_total()
        self.von_Misses()

    def torsional_constant(self):
        line_int = self.s1 / self.t_upper + self.s2 / self.t_spar1 + self.s3 / self.t_lower + self.s4/self.t_spar2 + \
                    + self.s5/self.t_spar2 + self.s6/self.t_upper
        Ao = self.h * self.L
        self.J = 4*pow(Ao,2)/line_int



def colorline(x, y, z= None, cmap=plt.get_cmap('copper'),
                  norm=plt.Normalize(-1.0, 1.0), linewidth=5, alpha=1.0):
    """
    Source:
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    if not hasattr(z, "__iter__"):
        z = np.array([z])

    z = np.asarray(z)
    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)
    ax = plt.gca()
    ax.add_collection(lc)

    return lc

def make_segments(x, y):
    """
    Source:
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

if __name__ == '__main__':
    h = 2
    L = 5
    t = 0.002
    Ixx = L*h**3/12 - (L-t)*(h-t)**3/12
    Iyy = h * L ** 3 / 12 - (h - t) * (L - t) ** 3 / 12
    chord1 = stresses(Ixx=Ixx,Iyy=Iyy,Ixx_str = Ixx, Iyy_str=Iyy,h=h,L=L,t_upper=t,t_spar1=t,t_spar2=t,t_lower=t)
    #chord1.shear_loads(Vx=500,Vy=750,T=1000)
    chord1.shear_loads(Vx=508742.32303978031, Vy=-657923.2689075109, T=-587877.5473618811)
    chord1.bending_loads(Mx = 5550,My= 0)
    chord1.compute_stresses()
    chord1.shear_flow_plotter(type="total", show=True)
    #chord1.sigma_plotter()

    #chord1.vm_plotter(show=True)
    #print(chord1.vm_max)


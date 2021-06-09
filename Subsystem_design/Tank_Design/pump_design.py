import numpy as np
import matplotlib.pyplot as plt
import math as m

print('\n \n ---> PUMP DESIGN <--- ')

"""
So mass flow is constant and equal to m_ff = rho * A * V
I can therefore compute speed that the pump needs to provide. 
Also, the pump will have to provide a certain pressure for the hydrogen to get into the combustion chamber,
therefore I get the value of the pressure at the combustion chamber, I add the pressure loss due to the piping
and then I can know the pressure which the pump needs to deliver.
"""

class booster_pump():
      def __init__(self,L,D,mf,h1,h2):
            """

            :param L: [m] maximum pipe length that we must account for
            :param D: [m] diameter of the pipe
            :param mf:[kg/s] mass flow
            :param h1:[m] booster pump height
            :param h2:[m] final height
            """
            self.L = L
            self.D = D
            self.A = m.pi * (self.D/2)**2 #[m^2] cross-sectional area of the pipe
            self.mf = mf
            self.h1 = h1
            self.h2 = h2
            self.g = 9.81

      def flow_velocity(self):
            """
            Based on mf = rho * A * v equation
            """
            self.rho = 71  # [kg/m^3] density of liquid hydrogen
            self.v = self.mf / self.rho / self.A # [m/s] required velocity of the flow

      def pressure_loss(self):
            """
            Calculate pressure loss based on Darcy-Weisbach equation
            :return:
            """
            mu = 13.92 * 10**(-6) #dummy value
            self.Re = self.rho * self.v * self.D /mu
            f_d = 64/self.Re#Darcy Friction factor
            p_L = f_d *self.rho/2* self.v/self.D #pressure loss per meter
            self.p_loss = p_L * self.L
            #self.p_loss = 0

      def initial_pressure(self):
            """
            Based on the pernoulli equation ew can calculate the necessary initial pressure.
            p1 + 1/2 rho V1^2 + rho * g * h1 = p2 + 1/2 rho V2^2 + p_loss + rho * g * h2
            :return:
            """
            self.v1 = 0 #[m/s] we assume to start from stationary
            self.v2 = self.v
            self.p2 = 3.45 #[bar] = inlet pressure of high pressure pump brewer
            self.p1 = self.p2 * 10**5 + self.p_loss + 1/2 * self.rho *(self.v2**2-self.v1**2) + self.rho * self.g * (self.h2 - self.h1) #Pa

      def compute_booster(self):
            self.flow_velocity()
            self.pressure_loss()
            self.initial_pressure()

if __name__ == '__main__':
    boost_pump_1 =  booster_pump(L = 10 ,D = 0.02, mf= 0.121, h1 = 0.75, h2 = 0)
    boost_pump_1.compute_booster()
    print((boost_pump_1.p2*10**5-boost_pump_1.p1)) #pressure difference in Pascals
    print(boost_pump_1.Re)
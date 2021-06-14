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
            self.A = m.pi * (self.D/2)**2   # [m^2] cross-sectional area of the pipe
            self.mf = mf
            self.h1 = h1
            self.h2 = h2
            self.g = 9.81
            self.p_tank = 1.2 *10**5 #[Pa] pressure of the tank

      def flow_velocity(self):
            """
            Based on mf = rho * A * v equation
            """
            self.rho = 71.1                         # [kg/m^3] density of liquid hydrogen
            self.v = self.mf / self.rho / self.A    # [m/s] required velocity of the flow

      def f_factor(self):
            """
            Calculated with Swamee-Jain equation
            https://www.pumpfundamentals.com/images/tutorial/friction%20loss-pipe.pdf
            :return:
            """
            epsilon = 0.00015 * 0.3048 #from ft to m
            self.f = 0.25/pow((m.log10(epsilon/self.D/3.7+5.74/pow(self.Re,0.9))),2)
      def pressure_loss(self):
            """
            Calculate pressure loss based on Darcy-Weisbach equation
            :return:
            """
            mu = 13.92 * 10**(-6) #dummy value
            self.Re = self.rho * self.v * self.D /mu
            f_d = 64/self.Re#Darcy Friction factor
            self.f_factor()
            f_d = self.f
            p_L = f_d *self.rho/2* self.v/self.D #pressure loss per meter
            self.p_loss = p_L * self.L
            #self.p_loss = 0

      def initial_pressure(self):
            """
            Based on the Bernoulli equation ew can calculate the necessary initial pressure.
            p1 + 1/2 rho V1^2 + rho * g * h1 = p2 + 1/2 rho V2^2 + p_loss + rho * g * h2
            :return:
            """
            self.v1 = 0                   # [m/s] we assume to start from stationary
            self.v2 = self.v
            self.p2 = 3.45                # [bar] = inlet pressure of high pressure pump brewer
            self.p1 = self.p2 * 10**5 + self.p_loss + 1/2 * self.rho * (self.v2**2-self.v1**2) + \
                      self.rho * self.g * (self.h2 - self.h1) # Pa
            self.dp = -(self.p2*10**5-self.p1)

      def compute_booster(self):
            self.flow_velocity()
            self.pressure_loss()
            self.initial_pressure()

      def work(self):
            """
            Calculate the work done by the pump. Assume isothermal cycle
            https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node22.html
            """
            miu = self.mf/2 #no of moles
            R = 8.31
            T = 19.75 #[K]
            self.efficiency = 0.6 #efficiency of the pump
            self.Work = miu * R * T * m.log(self.p1/self.p_tank)
            #this is the wotk required for each second so we need a constant supply of
            self.Power = abs(self.Work/self.efficiency) #[W]

      def power_pump(self):
            """
            Hydraulic power calculation
            https://www.engineeringtoolbox.com/pumps-power-d_505.html
            :return:
            """
            vol = self.mf/71 #[m^3/s)
            self.flow = vol * 3600 #[m^3/L]
            self.power_pump = self.flow * (self.p1-self.p_tank)/3.6/10**6/self.efficiency





if __name__ == '__main__':
    boost_pump_1 = booster_pump(L=(0.55-0.35)*17*1.5, D=0.0127, mf=0.121, h1=0.38, h2=0)
    boost_pump_1.compute_booster()
    boost_pump_1.work()
    boost_pump_1.power_pump()
    print((boost_pump_1.p2*10**5-boost_pump_1.p1))  # pressure difference in Pascals
    print(boost_pump_1.p1/10**5)
    print(boost_pump_1.Power,"W")
    print("Boost pump power",boost_pump_1.power_pump, "kW")
    ###Calculate power required of the high pressure pump from brewer table 4-16
    p1 = 345 #KPa
    p2 = 5068 #KPa
    flow = boost_pump_1.flow
    eff = 0.602
    power_HP = flow * (p2-p1)*1000 /3.6/pow(10,6)/eff
    print("High Pressure Pump Power",power_HP,"[kW]")


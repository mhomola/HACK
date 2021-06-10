import scipy.integrate as integrate
from scipy.interpolate import interp1d
from math import exp,log,inf
import numpy as np
import matplotlib.pyplot as plt

#Model to determine the T change due to aviation

class Climate_assess():

    def __init__(self):
        self.H = 30                                                  # Span of years considered
        self.RF_2CO2 = 3.7                                             # Radiative forcing corresponding to   [W/m^2]
                                                                       # a doubling of the concentration
        self.t_0 = 2035
        #self.t = self.t_0 + self.H
        self.dt = 0.1
        self.sensitivity = 3.0                                         # Sensitivity                              [K]
        self.alpha_t = 0.595
        self.tau_t1 = 8.4                                                                                          #[years]
        self.tau_t2 = 409.5                                            #                                            [year]


        'Efficacies'
        self.Eff_CO2 = 1                                               # Efficacy of CO2                        [-]
        self.Eff_CH_4 = 1.18                                           # Efficacy of CH4                        [-]
        self.Eff_O3 = 1.37                                             # Efficacy of O3                         [-]
        self.Eff_H2O = 1.14                                            # Efficacy of H2O                        [-]
        self.Eff_SO4 = 0.9                                             # Efficacy of SO2                        [-]
        self.Eff_soot = 0.7                                            # Efficacy of soot                       [-]
        self.Eff_contrails = 0.59                                      # Efficacy of contrails                  [-]
        'Constants for RF of CO2'
        self.EI_CO2 = 3.16                                             # Emissions index of CO2                 [kg/kgkerosene]
        self.E_CO2 = 0                                                 # Absolute CO2 emissions                 [kg]
        self.X_CO2_0 = 380                                             # Background concentration of CO2        [ppmv]
        self.alpha_1 = 0.067
        self.alpha_2 = 0.1135
        self.alpha_3 = 0.152
        self.alpha_4 = 0.0970
        self.alpha_5 = 0.041
        self.tau_1 = inf
        self.tau_2 = 313.8
        self.tau_3 = 79.8
        self.tau_4 = 18.8
        self.tau_5 = 1.7

        self.A_CO2 = 1.8 * 10**(-15)                                   # W/m^2
        self.alpha_c1 = 0.259
        self.alpha_c2 = 0.338
        self.alpha_c3 = 0.186
        self.tau_c1 = 172.9                                            # year
        self.tau_c2 = 18.51                                            # year
        self.tau_c3 = 1.186                                            # year

        'Contants for RF of NOx'
        self.A_CH4 = -5.16 * 10**-13                                   # Coefficient for CH4 production        [W/m^2/kgNOx]
                                                                  # from NOx emissions
        self.A_O3_L = -1.21 * 10**-13                                  # Coefficient for long-lived O3         [W/m^2/kgNOx]
                                                                  # production from NOx emissions
        self.tau_n = 12                                                # Perturbation lifetime                 [years]
        self.si = 0
        self.RF_E_ref_NOx = 1.01 * 10**-11                             # Radiative forcing due to NOx-O3s      [(W/m^2)/kgNOx]
                                                                  # per unit of NOx
        'Constants for RF of Water,soot and sulfate'
        self.EI_H2O = 1.26                                             # Emission index of H2O                 [kgH2O/kgkerosene]
        self.EI_soot = 2.0 * 10**-4                                    # Emission index of soot                {kgsoot/kgkerosene]
        self.EI_SO4 = 4.0 * 10**-5                                     # Emission index of SO4                 [kgSO4/kgkerosene]
        self.RF_E_ref_H2O = 7.43 * 10**(-15)                             # Radiative forcing due to H2O          [(W/m^2)/kgNOx]
                                                                  # per unit of H2O
        self.RF_E_ref_SO4 = -1.0 * 10**(-10)                             # Radiative forcing due to SO4          [(W/m^2)/kgNOx]
                                                                  # per unit of SO4
        self.RF_E_ref_soot = 5. * 10**(-10)                              # Radiative forcing due to soot          [(W/m^2)/kgNOx]
                                                                  # per unit of soot
        'Constants for RF of contrails'
        self.RF_L_ref_AIC = 2.21 * 10**-12                             # Radiative forcing due to contrails     [(W/m^2)/n mile]

    def emissions(self,m_k):
        e_CO2 = m_k * self.EI_CO2
        e_SO4 = m_k * self.EI_SO4
        e_H2O = m_k * self.EI_H2O
        e_soot = m_k * self.EI_soot
        e_NOx = m_k * (0.0986 * (1420000/101325)**0.4 - exp(772/194.4))
        return e_CO2,e_SO4,e_H2O,e_soot,e_NOx



    def emissions_level(self,e_CO2,e_H2O,e_NOx,e_soot,e_sulfate):
        self.eCO2 = e_CO2
        self.eH2O = e_H2O
        self.eNOx = e_NOx
        self.esoot = e_soot
        self.esulfate = e_sulfate
    def E(self,e ,U,time):
        '''
        :param e: Emission of a certain species in [kg]
        :param U: Number of missions/number of flights per year
        :return: Time-variant annual emissions rate
        '''
        sigma = np.zeros(len(time))
        sigma[time<=(self.t_0 + self.H)] = 1.
        sigma[time>(self.t_0 + self.H)] = 0.
        # if time <= (self.t_0 + self.H):
        #     sigma = 1
        # elif time > (self.t_0 + self.H):
        #     sigma = 0

        E = e * U * sigma

        return E

    '''Start of functions to compute RF of CO2'''

    def G_X_CO2_second_method(self,tau):
        alphas = np.array([self.alpha_c1, self.alpha_c2, self.alpha_c3])
        taus = np.array([self.tau_c1, self.tau_c2, self.tau_c3])
        A_coeff = self.A_CO2 * np.ones(len(tau))

        G_CO2 = 0
        for i in range(len(taus)):
            G_CO2 += alphas[i] * (np.exp(-tau/taus[i])-1)

        G_CO2 = (G_CO2 + np.ones(len(G_CO2))) * A_coeff

        return G_CO2

    '''Start of functions to compute RF of NOx (CH4, O3L and O3S)'''

    def s(self, compound):
        '''
        :param compound:Name of the species we want the RF of.
                        Can be CO2, CH4,O3L,O3S,soot,sulfate,H2O and contrails
        :return: forcing factor function for CH4 or O3S or contrails
        to take into account the change in RF of NOx with altitude
        '''
        if compound == 'CH4':
            forcing_fact = np.array([0.864, 0.865, 0.920777, 0.9592, 0.9624, 0.945, 0.9448, 0.9346, 0.938, 0.9132,
                                     0.9205, 0.9345, 0.9728, 1.0433, 1.1278, 1.187, 1.2034, 1.1885, 1.179])

            altitude = 0.3048 * np.array([0.0, 17582.52, 19617.886, 21263.886, 23302.9859, 24712.7742, 25293.3739,
                                          26420.091, 27505.983, 28509.118, 29562.49, 31242.0659, 33433.873,
                                          34645.2834, 35564.9354, 36792.6272, 37579.1538, 39084.3820, 41220.932])
            s = interp1d(altitude, forcing_fact, kind='cubic')
        if compound == 'O3S':
            forcing_fact = np.array([0.475, 0.476, 0.545, 0.5764, 0.6149, 0.666, 0.7119, 0.7029,
                                     0.7104, 0.7681, 0.8217, 0.8794, 0.9579, 0.9979, 1.060,
                                     1.1311, 1.222, 1.314, 1.4057, 1.497, 1.59, 1.671, 1.7542, 1.8288, 1.899])

            altitude = 0.3048 * np.array([0.0, 17734.057, 19270.4087, 20162.2675, 21487.656, 22736.4962,
                                          23452.5, 24338.567, 25615.992, 26827.0408, 27833.806, 28856.331, 30201.358,
                                          31624.234,
                                          32571.11, 33550.85, 34199.759, 34816.349, 35466.1077, 36338.689, 37290.789,
                                          38163.418,
                                          39184.116, 40136.381, 41251.765])
            s = interp1d(altitude, forcing_fact, kind='cubic')

        if compound == 'contrails':
            forcing_fact = np.array(
                [0.0243, 0.0244, -0.0071, 0.0108, 0.0978, 0.184, 0.2763, 0.3677, 0.459, 0.551, 0.642,
                 0.734, 0.826, 0.9176, 1.007, 1.101, 1.192, 1.28, 1.37, 1.467, 1.559, 1.651, 1.742, 1.83,
                 1.926, 2.016, 2.088, 1.983, 1.89, 1.798, 1.706, 1.613, 1.521, 1.429,
                 1.336, 1.24, 1.1524, 1.060, 0.968, 0.8837, 0.799])

            altitude = 0.3048 * np.array([0.0, 17688.384, 19554.874, 21704.664, 22710.769, 23702.416, 24512.913,
                                          25320.0093, 25860.9077, 26328.665, 26787.919, 27251.425, 27655.397,
                                          28080.632, 28490.7, 28880.073, 29297.654, 29679.514, 30083.487,
                                          30491.712, 30891.433, 31311.565, 31712.136, 32192.651, 32651.904,
                                          33214.7455, 33467.646, 34199.192, 34901.396, 35557.107,
                                          36189.857, 36823.456, 37445.149, 37810.001, 38128.927, 38452.105,
                                          38779.536, 39106.966, 39410.290, 40242.343, 41221.233])
            s = interp1d(altitude, forcing_fact, kind='cubic')

        return s
    def G_NOx(self,tau,compound):
        '''
        :param t: Time until we want to check the climate impact [year]
        :compound: Name of the species we want the RF of.
                    Can be CO2, CH4,O3L,O3S,soot,sulfate,H2O and contrails
        :return: Response function of CH4
        '''
        if compound == 'CH4':
            GNOx = self.A_CH4 * np.exp(-(tau)/self.tau_n)
        if compound == 'O3L':
            GNOx = self.A_O3_L * np.exp(-(tau) / self.tau_n)
        return GNOx

    def RF_(self,t,h,e,U,compound):
        '''
        :param t: Time until we want to check the climate impact [year]
        :param h: Altitude or altitudes at which we are in [m]
        :param e: Emission of a certain species in [kg]
        :param U: Number of missions/number of flights per year
        :param compound: Name of the species we want the RF of.
                         Can be CO2, CH4,O3L,O3S,soot,sulfate,H2O and contrails
        :return: Radiative forcing  [W/m^2]
        '''
        t_prime = np.arange(self.t_0, t + 1,1)
        E = self.E(e,U,t_prime)

        if compound == 'CO2':
            #RF = ((1 / log(2)) * np.log((self.X_CO2_0 + self.delta_X_CO2(t,e,U)) / self.X_CO2_0))
            RF = []
            delta_T = []
            tprime_a = self.t_0
            tprime_b = self.t_0 + 1
            for anual_E in E:
                tprime = np.arange(tprime_a,tprime_b +self.dt, self.dt)
                RF_CO2 = self.G_X_CO2_second_method(t - tprime)
                RF_CO2_conv = anual_E * integrate.simps(RF_CO2, tprime)
                delta_T_CO2 = self.G_T(t-tprime) * RF_CO2_conv
                delta_T_CO2_conv = integrate.simps(delta_T_CO2)
                RF.append(RF_CO2_conv)
                delta_T.append(delta_T_CO2_conv)
                tprime_a +=1
                tprime_b +=1
            RF = np.array(RF)
            delta_T = np.array(delta_T)
            print('The RF of CO2 is:',RF,len(RF))
            plt.plot(t_prime,RF,label='Radiative forcing of'+str(compound))
            plt.legend()
            plt.show()
        if compound == 'CH4':
            s_CH4 = self.s(compound='CH4')
            #RF = s_CH4(h) * integrate.simps(self.G_NOx(t-t_prime,compound)*self.E(e,U,t_prime),t_prime)*np.ones(len(t_prime))
            #RF = s_CH4(h) * self.G_NOx(t-t_prime,compound) * self.E(e,U,t_prime)
            RF = []
            delta_T = []
            tprime_a = self.t_0
            tprime_b = self.t_0 + 1
            for anual_E in E:
                tprime = np.arange(tprime_a,tprime_b +self.dt, self.dt)
                RF_CH4 = self.G_NOx(t-tprime,compound)
                RF_CH4_conv = s_CH4(h) * anual_E * integrate.simps(RF_CH4, tprime)
                delta_T_CH4 = self.G_T(t - tprime) * RF_CH4_conv
                delta_T_CH4_conv = integrate.simps(delta_T_CH4)
                RF.append(RF_CH4_conv)
                delta_T.append(delta_T_CH4_conv)
                tprime_a +=1
                tprime_b +=1
            RF = np.array(RF)
            delta_T = np.array(delta_T)
            print('The RF of O3-CH4 is:',RF)
            plt.plot(t_prime, RF, label='Radiative forcing of' + str(compound))
            plt.legend()
            plt.show()
        if compound == 'O3L':
            s_O3L = self.s(compound='CH4')
            #RF = s_O3L(h) * integrate.simps(self.G_NOx(t - t_prime, compound) * self.E(e, U,t_prime), t_prime)*np.ones(len(t_prime))
            #RF = s_O3L(h) * self.G_NOx(t - t_prime, compound) * self.E(e, U, t_prime)
            RF = []
            delta_T = []
            tprime_a = self.t_0
            tprime_b = self.t_0 + 1
            for anual_E in E:
                tprime = np.arange(tprime_a,tprime_b +self.dt, self.dt)
                RF_O3L = self.G_NOx(t-tprime,compound)
                RF_O3L_conv = s_O3L(h) * anual_E * integrate.simps(RF_O3L, tprime)
                delta_T_O3L = self.G_T(t - tprime) * RF_O3L_conv
                delta_T_O3L_conv = integrate.simps(delta_T_O3L)
                RF.append(RF_O3L_conv)
                delta_T.append(delta_T_O3L_conv)
                tprime_a +=1
                tprime_b +=1
            RF = np.array(RF)
            delta_T = np.array(delta_T)
            print('The RF of O3L is:',RF)
            plt.plot(t_prime, RF, label='Radiative forcing of' + str(compound))
            plt.legend()
            plt.show()
        if compound == 'O3S':
            s_O3S = self.s(compound='O3S')
            RF = s_O3S(h) * self.RF_E_ref_NOx * self.E(e,U,t_prime)
            delta_T = []
            tprime_a = self.t_0
            tprime_b = self.t_0 + 1
            for i in RF:
                tprime = np.arange(tprime_a, tprime_b + self.dt, self.dt)
                delta_T_O3S = self.G_T(t - tprime) * i
                delta_T_O3S_conv = integrate.simps(delta_T_O3S)
                delta_T.append(delta_T_O3S_conv)
                tprime_a += 1
                tprime_b += 1
            delta_T = np.array(delta_T)
            print('The RF of O3s is:',RF)
            plt.plot(t_prime, RF, label='Radiative forcing of' + str(compound))
            plt.legend()
            plt.show()
        if compound == 'H2O':
            RF = self.RF_E_ref_H2O * self.E(e,U,t_prime)
            delta_T = []
            tprime_a = self.t_0
            tprime_b = self.t_0 + 1
            for i in RF:
                tprime = np.arange(tprime_a, tprime_b + self.dt, self.dt)
                delta_T_H2O = self.G_T(t - tprime) * i
                delta_T_H2O_conv = integrate.simps(delta_T_H2O)
                delta_T.append(delta_T_H2O_conv)
                tprime_a += 1
                tprime_b += 1
            delta_T = np.array(delta_T)
            print('The RF of H2O is:',RF)
            plt.plot(t_prime, RF, label='Radiative forcing of' + str(compound))
            plt.legend()
            plt.show()
        if compound == 'soot':
            RF = self.RF_E_ref_soot * self.E(e,U,t_prime)
            delta_T = []
            tprime_a = self.t_0
            tprime_b = self.t_0 + 1
            for i in RF:
                tprime = np.arange(tprime_a, tprime_b + self.dt, self.dt)
                delta_T_soot = self.G_T(t - tprime) * i
                delta_T_soot_conv = integrate.simps(delta_T_soot)
                delta_T.append(delta_T_soot_conv)
                tprime_a += 1
                tprime_b += 1
            delta_T = np.array(delta_T)
            print('The RF of soot is:', RF)
            plt.plot(t_prime, RF, label='Radiative forcing of' + str(compound))
            plt.legend()
            plt.show()
        if compound == 'sulfate':
            RF = self.RF_E_ref_SO4 * self.E(e,U,t_prime)
            delta_T = []
            tprime_a = self.t_0
            tprime_b = self.t_0 + 1
            for i in RF:
                tprime = np.arange(tprime_a, tprime_b + self.dt, self.dt)
                delta_T_SO4 = self.G_T(t - tprime) * i
                delta_T_SO4_conv = integrate.simps(delta_T_SO4)
                delta_T.append(delta_T_SO4_conv)
                tprime_a += 1
                tprime_b += 1
            delta_T = np.array(delta_T)
            print('The RF of SO4 is:', RF)
            plt.plot(t_prime, RF, label='Radiative forcing of' + str(compound))
            plt.legend()
            plt.show()
        #if compound == 'contrails':
            #s_contrails = self.s(compound='contrails')
            #RF = s_contrails(h) * self.RF_L_ref_AIC * L

        return RF,delta_T

    def G_T (self,tau):

        #GT = (2.246/36.8) * np.exp(-(tau)/36.8)
        G_T = self.sensitivity * ((self.alpha_t/self.tau_t1) * np.exp(-tau/self.tau_t1) + ((1-self.alpha_t)/self.tau_t2)*np.exp(-tau/self.tau_t2))
        # plt.plot(tau,GT,label='With first formula')
        # plt.plot(tau,G_T,label= 'With second formula')
        # plt.legend()
        # plt.show()
        return G_T

    def RF_norm(self,t,h,e_CO2, e_H2O, e_NOx, e_soot, e_sulfate,U):
        '''
        :param t: Time until we want to check the climate impact [year]
        :param h: Altitude or altitudes at which we are in [m]
        :param e_CO2: Emission of CO2 in [kg]
        :param e_H2O: Emission of H2O in [kg]
        :param e_NOx: Emission of NOx in [kg]
        :param e_soot: Emission of soot in [kg]
        :param e_sulfate: Emission of sulfate in [kg]
        :param U: Number of missions/number of flights per year
        :param compound: Name of the species we want the RF of.
                         Can be CO2, CH4,O3L,O3S,soot,sulfate,H2O and contrails
        :return: Normalized Radiative forcing  [-]
        '''
        self.emissions_level(e_CO2, e_H2O, e_NOx, e_soot, e_sulfate)

        RFnorm_O3S = self.RF_(t, h, self.eNOx, U, compound='O3S') * self.Eff_O3
        RFnorm_H2O = self.RF_(t, h, self.eH2O, U, compound='H2O') * self.Eff_H2O
        RFnorm_soot = self.RF_(t, h, self.esoot, U, compound='soot') * self.Eff_soot
        RFnorm_CO2 = self.RF_(t, h, self.eCO2, U, compound='CO2')
        RFnorm_CH4 = (self.RF_(t, h, self.eNOx, U, compound='CH4') * self.Eff_CH_4)
        RFnorm_O3L = (self.RF_(t, h, self.eNOx, U, compound='O3L') * self.Eff_O3)
        #RFnorm_sulfate = self.RF(t, h, e, U, compound='sulfate') * self.Eff_SO4
        #RFnorm_contrails = self.RF(t, h, e, U, compound='contrails') * self.Eff_contrails

        RF_norm = np.zeros(len(RFnorm_O3S))
        for i in range(len(RFnorm_O3S)):
            RF_norm[i] = RFnorm_O3S[i] + RFnorm_H2O[i] + RFnorm_soot[i] + RFnorm_CH4[i] + RFnorm_O3L[i] + RFnorm_CO2[i]


        '''To check compliance with requirements'''
        RFnorm = (1/self.RF_2CO2) * RF_norm
        t_prime = np.arange(self.t_0, t + 1, 1)

        return RFnorm

    def delta_T(self,t,h,e_CO2, e_H2O, e_NOx, e_soot, e_sulfate,U,plot=False):
        '''
        :param t: Time until we want to check the climate impact [year]
        :param h: Altitude or altitudes at which we are in [m]
        :param e_CO2: Emission of CO2 in [kg]
        :param e_H2O: Emission of H2O in [kg]
        :param e_NOx: Emission of NOx in [kg]
        :param e_soot: Emission of soot in [kg]
        :param e_sulfate: Emission of sulfate in [kg]
        :param U: Number of missions/number of flights per year
        :param plot: If we want to plot the change in Temperature per year
        :return: Change in temperature over the span of years selected [K]
        '''

        'To check with climate impact requirement'
        self.emissions_level(e_CO2, e_H2O, e_NOx, e_soot, e_sulfate)

        _,delta_T_O3S = self.RF_(t, h, self.eNOx, U, compound='O3S')
        _,delta_T_H2O = self.RF_(t, h, self.eH2O, U, compound='H2O')
        _,delta_T_soot = self.RF_(t, h, self.esoot, U, compound='soot')
        _,delta_T_CO2 = self.RF_(t, h, self.eCO2, U, compound='CO2')
        _,delta_T_CH4 = self.RF_(t, h, self.eNOx, U, compound='CH4')
        _,delta_T_O3L = self.RF_(t, h, self.eNOx, U, compound='O3L')
        delta_T_tot = (1/self.RF_2CO2)*(delta_T_CO2 + delta_T_CH4*self.Eff_CH_4 + delta_T_H2O*self.Eff_H2O + delta_T_O3L*self.Eff_O3 + delta_T_O3S*self.Eff_O3 + delta_T_soot*self.Eff_soot)
        print(len(delta_T_tot))
        t_prime = np.arange(self.t_0, t + 1,1)

        if plot == True:
            plt.plot(t_prime,delta_T_O3S,color = 'tab:red',label= 'O3S')
            plt.plot(t_prime, delta_T_H2O, color='tab:orange', label='H2O')
            plt.plot(t_prime, delta_T_soot, color='tab:blue', label='soot')
            plt.plot(t_prime, delta_T_CO2, color='tab:green', label='CO2')
            plt.plot(t_prime, delta_T_CH4, color='black', label='CH4')
            plt.plot(t_prime, delta_T_O3L, color='grey', label='O3L')
            plt.plot(t_prime, delta_T_tot, color='navy', label='Total')

            plt.xlabel('[years]')
            plt.ylabel('$\Delta_T$ [K]')
            plt.legend()
            plt.show()

        #del_T = integrate.simps(self.G_T(t-t_prime) * RF_norm, t_prime)

        return delta_T_tot

    def ATR(self,t,h,e_CO2, e_H2O, e_NOx, e_soot, e_sulfate,U,plot=False):
        '''
        :param t: Time until we want to check the climate impact [year]
        :param h: Altitude or altitudes at which we are in [m]
        :param e_CO2: Emission of CO2 in [kg]
        :param e_H2O: Emission of H2O in [kg]
        :param e_NOx: Emission of NOx in [kg]
        :param e_soot: Emission of soot in [kg]
        :param e_sulfate: Emission of sulfate in [kg]
        :param U: Number of missions/number of flights per year
        :return: Average temperature response for a certain flight phase [K]
        '''
        #print('The deltaT is:',self.delta_T(t, h, e_CO2, e_H2O, e_NOx, e_soot, e_sulfate, U))
        #H_array = np.arange(0.,self.H + self.dt,self.dt)
        #H_array = self.H
        #ATR = (1/self.H) * integrate.simps(self.delta_T(t,h,e_CO2, e_H2O, e_NOx, e_soot, e_sulfate,U),H_array)
        inte_interval = np.arange(0,self.H + 1, 1)
        ATR = (1/self.H)*integrate.simps(self.delta_T(t, h, e_CO2, e_H2O, e_NOx, e_soot, e_sulfate, U,plot),inte_interval)
        return ATR

if __name__ == '__main__':
    climate = Climate_assess()
    climate.s(compound='O3S')
    U = 50                 # Number of flights in a year

    print('Start analysis for LTO')
    'To plot the change in CO2 concentration in ppmv per year'

    #clim =climate.delta_X_CO2(t = 2135,e=50,U=U,plot= True)
    #print(clim)
    # tau = np.arange(climate.t_0, 2195 + climate.dt, climate.dt)
    # climate.G_T(2195- tau)

    e_CO2, e_SO4, e_H2O, e_soot, e_NOx = climate.emissions(13300)
    # print(e_CO2,e_SO4,e_H2O,e_soot,e_NOx)
    ATR = climate.ATR(t=2135,h=6000,e_CO2= e_CO2, e_H2O= e_H2O, e_NOx=1000,e_soot= e_soot,e_sulfate=e_SO4,U=U,plot = True)
    print('The average temperature response, A_100, for the LTO of the HACK is:',ATR,'[K]')
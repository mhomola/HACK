
from Subsystem_design.common_constants import Constants
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from math import exp,log,inf
import numpy as np
import matplotlib.pyplot as plt

#Model to determine the T change due to aviation

class Climate_assess():

    def __init__(self,t):
        self.H = 40                                                # Span of years considered(larger than 30 years)
        self.H_int = 100
        self.RF_2CO2 = 3.7                                             # Radiative forcing corresponding to   [W/m^2]
                                                                       # a doubling of the concentration
        self.t_0 = 2035#1940
        self.dt = 1
        self.t = t
        self.t_prime = np.arange(self.t_0, t + self.dt, self.dt)
        self.sensitivity = 3.0                                         # Sensitivity                              [K]
        self.alpha_t = 0.595
        self.tau_t1 = 8.4                                                                                          #[years]
        self.tau_t2 = 409.5                                            #                                            [year]
        self.growth_rate = 0.015                                      # Growth rate of number of flights between 2017-2040
        self.H2_flights_2035 = 0.08                                   # Percentage of flights on H2 in 2035
        self.H2_flights_2050 = 0.4                                     # Percentage of H2 flights in 2050
        self.flights_EU_2024 = 11411000                                # Number of flights in 2024 post corona, best scenario
        self.narrow_body_2025 = 0.64                                   # Percentage of flights narrow-body pass A/C take over in 2025
        self.narrow_body_2040 = 0.63                                   # Percentage of flights narrow-body pass A/C take over in 2040
        self.single_aisle_share = np.array([0.43,0.184,0.203,0.079,0.08,0.018]) #Assumed to stay constant over the years
        self.growth_rates = np.array([0.045,0.017,0.02,0.029,0.036,0.036,0.032])
        self.world_flights_2019 = 38.6 *10**6
        self.Narrow_body_flights_percent = 0.56

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
        self.X_CO2_0 = 311                                             # Background concentration of CO2        [ppmv]
        self.alpha_1 = 0.067/10**12#0.067                               #                                        [ppmv/kg]
        self.alpha_2 = 0.1135/10**12#0.1135                             #                                        [ppmv/kg]
        self.alpha_3 = 0.152/10**12#0.152                               #                                        [ppmv/kg]
        self.alpha_4 = 0.0970/10**12#0.0970                             #                                        [ppmv/kg]
        self.alpha_5 = 0.041/10**12#0.041
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
        self.EI_NOx = 0.014                                            # Emission index of NOx                 [kgH2O/kgkerosene]
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
        'In kg. 1kg = 1*10**-9 Tg'
        e_CO2 = m_k * self.EI_CO2
        e_SO4 = m_k * self.EI_SO4
        e_H2O = m_k * self.EI_H2O
        e_soot = m_k * self.EI_soot
        e_NOx = m_k * self.EI_NOx
        return e_CO2,e_SO4,e_H2O,e_soot,e_NOx



    def emissions_level(self,e_CO2,e_H2O,e_NOx,e_soot,e_sulfate):
        self.eCO2 = e_CO2
        self.eH2O = e_H2O
        self.eNOx = e_NOx
        self.esoot = e_soot
        self.esulfate = e_sulfate

    def narrow_body_percentage(self):
        percentages = np.array([self.narrow_body_2025,self.narrow_body_2040])
        time = np.array([2025.,2040.])

        narrow_body_percen = interp1d(time,percentages,fill_value='extrapolate')

        return narrow_body_percen

    def H2_percentage_trend(self):

        percentages= np.array([self.H2_flights_2035,0.1,0.18,self.H2_flights_2050,0.42,0.41,0.30,0.])
        hack_percentages = np.linspace(1.0, 0.0, len(self.t_prime[self.t_prime<=self.t_0 + self.H]))
        time = np.array([2035,2040,2045,2050,2055,2060,2065,self.t_0 + self.H])
        H2_perce = interp1d(time,percentages)
        H2_flights_perc_year = np.ones(len(self.t_prime))
        H2_flights_perc_year[self.t_prime<=self.t_0 + self.H] = H2_perce(self.t_prime[self.t_prime<=self.t_0 + self.H])* hack_percentages
        H2_flights_perc_year[self.t_prime > self.t_0 + self.H] = 0.
        #H2_flights_perc_year = H2_flights_perc_year * hack_percentages
        return H2_flights_perc_year

    def market_share_per_region(self):
        #Regions:
        # Asia/Pacific, North America, Europe, Latin America, Middle East, Africa
        share_2021 = np.array([0.3864,0.2272,0.2362,0.0571,0.0741,0.019])
        share_2016 = np.array([0.35,0.24,0.26,0.07,0.05,0.02])
        shares = np.vstack((share_2016,share_2021))
        time = np.array([2016,2021])
        share_func_Asia = interp1d(time,shares[:,0])
        share_func_NAmerica = interp1d(time, shares[:, 1])
        share_func_Europe = interp1d(time, shares[:, 2])
        share_func_Latin = interp1d(time, shares[:, 3])
        share_func_MiddleE = interp1d(time, shares[:, 4])
        share_func_Africa = interp1d(time, shares[:, 5])


        return share_func_Asia,share_func_NAmerica,share_func_Europe,share_func_Latin,share_func_MiddleE,share_func_Africa



    def total_flights_per_region_2035(self):
        # Get narrow_body flights in the world in 2019
        share_func_Asia, share_func_NAmerica, share_func_Europe, share_func_Latin, share_func_MiddleE, share_func_Africa = self.market_share_per_region()
        #Asia
        Asia_flights_2023 = share_func_Asia(2019) * self.world_flights_2019
        Asia_flights_2035 = Asia_flights_2023 * (1+self.growth_rates[0])**12

        #North America
        NA_flights_2023 = share_func_NAmerica(2019) * self.world_flights_2019
        NA_flights_2035 = NA_flights_2023 * (1 + self.growth_rates[1])**12

        #Europe
        Europe_flights_2024 = share_func_Europe(2019) * self.world_flights_2019
        Europe_flights_2035 = Europe_flights_2024 * (1 + self.growth_rates[2])**11

        #Latin America
        Latin_flights_2023 = share_func_Latin(2019) * self.world_flights_2019
        Latin_flights_2035 = Latin_flights_2023 * (1 + self.growth_rates[3])**12

        #Middle East
        MiddleE_flights_2024 = share_func_MiddleE(2019) * self.world_flights_2019
        MiddleE_flights_2035 = MiddleE_flights_2024 * (1 + self.growth_rates[4])**11

        #Africa
        Africa_flights_2024 = share_func_Africa(2019) * self.world_flights_2019
        Africa_flights_2035 = Africa_flights_2024 * (1 + self.growth_rates[5])**11

        #World
        World_flights_2035 = np.sum([Asia_flights_2035,NA_flights_2035,Europe_flights_2035,Latin_flights_2035,MiddleE_flights_2035,Africa_flights_2035])

        return Asia_flights_2035,NA_flights_2035,Europe_flights_2035,Latin_flights_2035,MiddleE_flights_2035,Africa_flights_2035, World_flights_2035

    def narrow_body_flights_per_region(self):
        _,_, share_func_Europe,_, _,_ = self.market_share_per_region()
        Asia_flights_2035, NA_flights_2035, Europe_flights_2035, Latin_flights_2035, MiddleE_flights_2035, Africa_flights_2035, World_flights_2035 = self.total_flights_per_region_2035()

        #Asia
        Asia_tot_flights = np.linspace(Asia_flights_2035,Asia_flights_2035 * (1 + self.growth_rates[0]) ** len(self.t_prime),len(self.t_prime))
        #North America
        NA_tot_flights = np.linspace(NA_flights_2035,NA_flights_2035 * (1 + self.growth_rates[1]) ** len(self.t_prime),len(self.t_prime))
        #Europe
        Europe_tot_flights = np.linspace(Europe_flights_2035, Europe_flights_2035 * (1 + self.growth_rates[2]) ** len(self.t_prime),
                                     len(self.t_prime))
        #Latin
        Latin_tot_flights = np.linspace(Latin_flights_2035, Latin_flights_2035 * (1 + self.growth_rates[3]) ** len(self.t_prime),
                                     len(self.t_prime))
        #Middle East
        MiddleE_tot_flights = np.linspace(MiddleE_flights_2035, MiddleE_flights_2035 * (1 + self.growth_rates[4]) ** len(self.t_prime),
                                     len(self.t_prime))
        #Africa
        Africa_tot_flights = np.linspace(Africa_flights_2035, Africa_flights_2035 * (1 + self.growth_rates[5]) ** len(self.t_prime),
                                     len(self.t_prime))
        #World
        World_tot_flights = np.sum(np.vstack((Asia_tot_flights,NA_tot_flights,Europe_tot_flights,Latin_tot_flights,MiddleE_tot_flights,Africa_tot_flights)),axis=0)
        print('World_tot_flights',World_tot_flights)
        #Narrow_body_flights
        #Narrow_body_flights_percent = (self.narrow_body_percentage()(2019) / self.single_aisle_share[2]) * share_func_Europe(2019)
        World_Narrow_body_flights = self.Narrow_body_flights_percent * World_tot_flights
        print('World_Narrow_body_flights',World_Narrow_body_flights)
        #Asia
        self.Asia_Narrow_body_flights = World_Narrow_body_flights * self.single_aisle_share[0]
        #North America
        self.NA_Narrow_body_flights = World_Narrow_body_flights * self.single_aisle_share[1]
        #Europe
        self.Europe_Narrow_body_flights = World_Narrow_body_flights * self.single_aisle_share[2]
        #Latin
        self.Latin_Narrow_body_flights = World_Narrow_body_flights * self.single_aisle_share[3]
        #Middle East
        self.MiddleE_Narrow_body_flights = World_Narrow_body_flights * self.single_aisle_share[4]
        #Africa
        self.Africa_Narrow_body_flights = World_Narrow_body_flights * self.single_aisle_share[5]

        return World_Narrow_body_flights


    def number_aircraft_H2(self):

        # tot_flights_2035 = self.flights_EU_2024 * (1+self.growth_rate)**11
        # tot_flights = np.linspace(tot_flights_2035,tot_flights_2035*(1+self.growth_rate)**len(self.t_prime),len(self.t_prime))
        #
        # narrow_body_percent = self.narrow_body_percentage()
        # narrow_body_flights = tot_flights * narrow_body_percent(self.t_prime)

        H2_percentage = self.H2_percentage_trend()
        N_aircraft_H2 = self.narrow_body_flights_per_region() * H2_percentage

        return N_aircraft_H2


    def E(self,e ,U,time):
        '''
        :param e: Emission of a certain species in [kg]
        :param U: Number of missions/number of flights per year
        :return: Time-variant annual emissions rate
        '''
        sigma = np.zeros(len(time))
        sigma[time<=self.t_0 + self.H] = 1.
        sigma[time>self.t_0 + self.H] = 0.


        E = e * U * sigma

        return E

    def verification(self):
        E_CO2 = 10**9 * np.array([7.7,8.3,9.0,9.7,10.5,11.3,12.3,13.2,14.3,15.4,16.7,18.0,19.4,21.0,22.7,24.5,26.5,28.6,
                          30.9,33.3,36.,39.5,43.5,45.9,48.,51.3,55.6,65.6,74.3,77.8,78.0,90.,96.0,96.1,96.4,102.1,
                          105.7,110.1,110.9,109.3,110.5,112.0,119.5,123.4,129.9,135.6,141.4,146.5,146.9,143.4,142.0,
                          144.1,150.0,154.3])           #kg
        CO2_concentration = np.array([0.05,0.35,1.13,1.22,1.37])                                #[ppmv]
        years = np.array([1950,1970,1990,1992,1995])
        RF_CO2 = np.array([0.001,0.007,0.021,0.022,0.024])
        delta_T = np.array([0.0,0.001,0.003,0.004,0.004])

        return E_CO2,CO2_concentration,years,RF_CO2,delta_T

    def errors(self,e,U):
        '''Errors in CO2 concentration'''
        Conc_study = self.verification()[1]
        Conc_model = self.delta_X_CO2(e,U)[self.t_prime==self.verification()[2]]

        Errors_conc = np.absolute(Conc_study - Conc_model)/Conc_study
        Error_conc = np.max(Errors_conc)




        return Error_conc

    '''Start of functions to compute RF of CO2'''
    def G_X_CO2(self,tau):
        '''
        :param t: Time until we want to check the climate impact [year]
        :return: Impulse responses function G_XCO2
        '''
        alphas = np.array([self.alpha_1, self.alpha_2, self.alpha_3, self.alpha_4, self.alpha_5])
        taus = np.array([self.tau_1, self.tau_2, self.tau_3, self.tau_4, self.tau_5])
        G_XCO2 = 0

        for i in range(len(taus)):
            G_XCO2 += alphas[i] * np.exp(-(tau) / taus[i])


        return G_XCO2


    def delta_X_CO2(self,e,U,plot=False):
        '''
        :param t: Time until we want to check the climate impact [year]
        :return: Change in atmospheric concentration of CO2
        '''


        delt_XCO2_per_year = np.convolve(self.E(e,U,self.t_prime),self.G_X_CO2(self.t_prime-self.t_0),'full')
        delt_XCO2_per_year = delt_XCO2_per_year[:self.t_prime.size]
        #delt_XCO2_per_year = np.convolve(self.verification()[0], self.G_X_CO2(self.t_prime - self.t_0), 'full')
        #delt_XCO2_per_year = delt_XCO2_per_year[:self.t_prime.size]
        # if plot == True:
        #     #plt.plot(self.t_prime,self.E(e,U,self.t_prime),label = 'Emissions per year')
        #     #plt.plot(self.t_prime,self.G_X_CO2(self.t_prime-self.t_0),label='Impulse function of CO2')
        #     plt.plot(self.t_prime, delt_XCO2_per_year,color = 'tab:red',label='Change in atmospheric CO2 concentration from model')
        #     #plt.plot(self.verification()[2],self.verification()[1],marker= '+',label='Change in atmospheric CO2 concentration from study')
        #     plt.xlabel('years')
        #     plt.ylabel('$\Delta_{CO2 concentration}$ [ppmv]')
        #     plt.legend()
        #     plt.show()

        #print('$delta_{XCO2} with quad rule:',delta_XCO2_1)
        #delta_XCO2 = integrate.simps(self.G_X_CO2(t - t_prime) * self.E(e, U, t_prime), t_prime)

        return delt_XCO2_per_year

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
    def s(self,compound):
        '''
        :param compound:Name of the species we want the RF of.
                        Can be CO2, CH4,O3L,O3S,soot,sulfate,H2O and contrails
        :return: forcing factor function for CH4 or O3S or contrails
        to take into account the change in RF of NOx with altitude
        '''
        if compound == 'CH4':
            forcing_fact = np.array([0.864,0.865,0.920777,0.9592,0.9624,0.945,0.9448,0.9346,0.938,0.9132,
                                 0.9205,0.9345,0.9728,1.0433,1.1278,1.187,1.2034,1.1885,1.179])

            altitude = 0.3048 * np.array([0.0,17582.52,19617.886,21263.886,23302.9859,24712.7742,25293.3739,
                             26420.091,27505.983,28509.118,29562.49,31242.0659,33433.873,
                             34645.2834,35564.9354,36792.6272,37579.1538,39084.3820,41220.932])
            s = interp1d(altitude,forcing_fact,kind='cubic')
        if compound == 'O3S':
            forcing_fact = np.array([0.475,0.476,0.545,0.5764,0.6149,0.666,0.7119,0.7029,
                                     0.7104,0.7681,0.8217,0.8794,0.9579,0.9979,1.060,
                                     1.1311,1.222,1.314,1.4057,1.497,1.59,1.671,1.7542,1.8288,1.899])

            altitude = 0.3048 * np.array([0.0,17734.057,19270.4087,20162.2675,21487.656,22736.4962,
                                 23452.5,24338.567,25615.992,26827.0408,27833.806,28856.331,30201.358,31624.234,
                                 32571.11,33550.85,34199.759,34816.349,35466.1077,36338.689,37290.789,38163.418,
                                 39184.116,40136.381,41251.765])
            s = interp1d(altitude, forcing_fact, kind='cubic')

        if compound == 'contrails':
            forcing_fact = np.array([0.0243,0.0244,-0.0071,0.0108,0.0978,0.184,0.2763,0.3677,0.459,0.551,0.642,
                                     0.734,0.826,0.9176,1.007,1.101,1.192,1.28,1.37,1.467,1.559,1.651,1.742,1.83,
                                     1.926,2.016,2.088,1.983,1.89,1.798,1.706,1.613,1.521,1.429,
                                     1.336,1.24,1.1524,1.060,0.968,0.8837,0.799])

            altitude = 0.3048 * np.array([0.0,17688.384,19554.874,21704.664,22710.769,23702.416,24512.913,
                                 25320.0093,25860.9077,26328.665,26787.919,27251.425,27655.397,
                                 28080.632,28490.7,28880.073,29297.654,29679.514,30083.487,
                                 30491.712,30891.433,31311.565,31712.136,32192.651,32651.904,
                                 33214.7455,33467.646,34199.192,34901.396,35557.107,
                                 36189.857,36823.456,37445.149,37810.001,38128.927,38452.105,
                                 38779.536,39106.966,39410.290,40242.343,41221.233])
            s = interp1d(altitude,forcing_fact,kind='cubic')

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


    def RF_(self,h,e,U,compound):
        '''
        :param t: Time until we want to check the climate impact [year]
        :param h: Altitude or altitudes at which we are in [m]
        :param e: Emission of a certain species in [kg]
        :param U: Number of missions/number of flights per year
        :param compound: Name of the species we want the RF of.
                         Can be CO2, CH4,O3L,O3S,soot,sulfate,H2O and contrails
        :return: Radiative forcing  [W/m^2]
        '''

        if compound == 'CO2':
            RF = ((1 / log(2)) * np.log((self.X_CO2_0 + self.delta_X_CO2(e,U)) / self.X_CO2_0))

            #print('The RF of CO2 is:',RF,len(RF))

            #plt.plot(self.t_prime,RF,label='Radiative forcing of '+ str(compound) +' from model')
            #plt.plot(self.verification()[2],self.verification()[3],marker='+',label='Radiative forcing of ' + str(compound) + ' from study')
            #plt.xlabel('years')
            # plt.ylabel('RF [W/m^2]')
            # plt.legend()
            # plt.show()
        if compound == 'CH4':
            s_CH4 = self.s(compound='CH4')

            RF = s_CH4(h) * np.convolve(self.E(e,U,self.t_prime),self.G_NOx(self.t_prime - self.t_0, compound),'full')
            RF = RF[:self.t_prime.size]

            #print('The RF of O3-CH4 is:',RF,len(RF))
            # plt.plot(self.t_prime, RF, label='Radiative forcing of' + str(compound))
            # plt.legend()
            # plt.show()
        if compound == 'O3L':
            s_O3L = self.s(compound='CH4')
            RF = s_O3L(h) * np.convolve(self.E(e,U,self.t_prime),self.G_NOx(self.t_prime - self.t_0, compound),'full')
            RF = RF[:self.t_prime.size]

            #print('The RF of O3L is:',RF,len(RF))
            # plt.plot(self.t_prime, RF, label='Radiative forcing of' + str(compound))
            # plt.legend()
            # plt.show()
        if compound == 'O3S':
            s_O3S = self.s(compound='O3S')
            RF = s_O3S(h) * self.RF_E_ref_NOx * self.E(e,U,self.t_prime)

            #print('The RF of O3s is:',RF)
            # plt.plot(self.t_prime, RF, label='Radiative forcing of' + str(compound))
            # plt.legend()
            # plt.show()
        if compound == 'H2O':
            RF = self.RF_E_ref_H2O * self.E(e,U,self.t_prime)
            #print('The RF of H2O is:',RF)
            # plt.plot(self.t_prime, RF, label='Radiative forcing of' + str(compound))
            # plt.legend()
            # plt.show()
        if compound == 'soot':
            RF = self.RF_E_ref_soot * self.E(e,U,self.t_prime)
            #print('The RF of soot is:', RF)
            # plt.plot(self.t_prime, RF, label='Radiative forcing of' + str(compound))
            # plt.legend()
            # plt.show()
        if compound == 'sulfate':
            RF = self.RF_E_ref_SO4 * self.E(e,U,self.t_prime)
            #print('The RF of SO4 is:', RF)
            # plt.plot(self.t_prime, RF, label='Radiative forcing of' + str(compound))
            # plt.legend()
            # plt.show()
        #if compound == 'contrails':
            #s_contrails = self.s(compound='contrails')
            #RF = s_contrails(h) * self.RF_L_ref_AIC * L

        return RF

    def G_T(self,tau):

        GT = (2.246/36.8) * np.exp(-(tau)/36.8)
        #G_T = self.sensitivity * ((self.alpha_t/self.tau_t1) * np.exp(-tau/self.tau_t1) + ((1-self.alpha_t)/self.tau_t2)*np.exp(-tau/self.tau_t2))
        # plt.plot(tau,GT,label='With first formula')
        # plt.plot(tau,G_T,label= 'With second formula')
        # plt.legend()
        # plt.show()
        return GT

    def RF_norm(self,h,e_CO2, e_H2O, e_NOx, e_soot, e_sulfate,U):
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

        RFnorm_O3S = self.RF_(h, self.eNOx, U, compound='O3S') * self.Eff_O3
        RFnorm_H2O = self.RF_(h, self.eH2O, U, compound='H2O') * self.Eff_H2O
        RFnorm_soot = self.RF_(h, self.esoot, U, compound='soot') * self.Eff_soot
        RFnorm_CO2 = self.RF_(h, self.eCO2, U, compound='CO2')
        RFnorm_CH4 = (self.RF_(h, self.eNOx, U, compound='CH4') * self.Eff_CH_4)
        RFnorm_O3L = (self.RF_(h, self.eNOx, U, compound='O3L') * self.Eff_O3)
        #RFnorm_sulfate = self.RF(t, h, e, U, compound='sulfate') * self.Eff_SO4
        #RFnorm_contrails = self.RF(t, h, e, U, compound='contrails') * self.Eff_contrails
        RF_norm = np.zeros(len(RFnorm_O3S))
        for i in range(len(RFnorm_O3S)):
            RF_norm[i] = RFnorm_O3S[i] + RFnorm_H2O[i] + RFnorm_soot[i] + RFnorm_CH4[i] + RFnorm_O3L[i]

        '''To check compliance with requirements'''
        RFnorm = RFnorm_CO2 + ((1/self.RF_2CO2) * RF_norm)

        return RFnorm

    def delta_T(self,h,e_CO2, e_H2O, e_NOx, e_soot, e_sulfate,U,plot=False):
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

        self.emissions_level(e_CO2, e_H2O, e_NOx, e_soot, e_sulfate)

        delta_T_O3S = np.convolve(self.RF_(h, self.eNOx, U, compound='O3S'),self.G_T(self.t_prime - self.t_0),'full')
        delta_T_O3S = (self.Eff_O3/self.RF_2CO2)*delta_T_O3S[:self.t_prime.size]

        delta_T_H2O = np.convolve(self.RF_(h, self.eH2O, U, compound='H2O'),self.G_T(self.t_prime - self.t_0),'full')
        delta_T_H2O = (self.Eff_H2O/self.RF_2CO2)*delta_T_H2O[:self.t_prime.size]

        delta_T_soot = np.convolve(self.RF_(h, self.esoot, U, compound='soot'),self.G_T(self.t_prime - self.t_0),'full')
        delta_T_soot =(self.Eff_soot/self.RF_2CO2)*delta_T_soot[:self.t_prime.size]

        delta_T_CO2 = np.convolve(self.RF_(h, self.eCO2, U, compound='CO2'),self.G_T(self.t_prime - self.t_0),'full')
        delta_T_CO2 = delta_T_CO2[:self.t_prime.size]

        delta_T_CH4 = np.convolve(self.RF_(h, self.eNOx, U, compound='CH4'),self.G_T(self.t_prime - self.t_0),'full')
        delta_T_CH4 = (self.Eff_CH_4/self.RF_2CO2)*delta_T_CH4[:self.t_prime.size]

        delta_T_O3L = np.convolve(self.RF_(h, self.eNOx, U, compound='O3L'),self.G_T(self.t_prime - self.t_0),'full')
        delta_T_O3L = (self.Eff_O3/self.RF_2CO2)* delta_T_O3L[:self.t_prime.size]

        delta_T_tot = delta_T_CO2 + delta_T_CH4+ delta_T_H2O  + delta_T_O3L  + delta_T_O3S  + delta_T_soot


        'To check with climate impact requirement'
        #t_prime_2 = np.linspace(self.t_0, self.t, len(delta_T_O3S))
        #t_prime_3 = np.linspace(self.t_0,self.t,len(delta_T_CO2))
        #RF_norm = self.RF_norm(t,h,e_CO2, e_H2O, e_NOx, e_soot, e_sulfate,U)

        #delta_T = np.convolve(RF_norm,self.G_T(t-t_prime),'full')

        if plot == True:
            plt.plot(self.t_prime, delta_T_O3S, color='tab:red', label='O3S')
            plt.plot(self.t_prime, delta_T_H2O, color='tab:orange', label='H2O')
            plt.plot(self.t_prime, delta_T_soot, color='tab:blue', label='soot')
            plt.plot(self.t_prime, delta_T_CO2, color='tab:green', label='CO2_$\Delta_T$ from model')
           # plt.plot(self.verification()[2],self.verification()[4],marker='+',label='CO2_$\Delta_T$ from study')
            plt.plot(self.t_prime, delta_T_CH4, color='black', label='CH4')
            plt.plot(self.t_prime, delta_T_O3L, color='grey', label='O3L')
            plt.plot(self.t_prime,delta_T_tot,color = 'yellow',label= 'Change in temperature per year')
            plt.xlabel('years')
            plt.ylabel('$\Delta_T$ [K]')
            plt.legend()
            plt.show()

        #del_T = integrate.simps(self.G_T(t-t_prime) * RF_norm, t_prime)

        return delta_T_tot

    def ATR(self,h,e_CO2, e_H2O, e_NOx, e_soot, e_sulfate,U,plot=False):
        '''
        :param t: Time until we want to check the climate impact [year]
        :param h: Altitude or altitudes at which we are in [m]
        :param e_CO2: Emission of CO2 in [Tg]
        :param e_H2O: Emission of H2O in [Tg]
        :param e_NOx: Emission of NOx in [Tg]
        :param e_soot: Emission of soot in [Tg]
        :param e_sulfate: Emission of sulfate in [Tg]
        :param U: Number of missions/number of flights per year
        :return: Average temperature response for a certain flight phase [K]
        '''

        ATR = (1/self.H_int) * integrate.simps(self.delta_T(h,e_CO2, e_H2O, e_NOx, e_soot, e_sulfate,U,plot),self.t_prime-self.t_0)

        return ATR




if __name__ == '__main__':
    climate = Climate_assess(t = 2135) #1995
    climate.s(compound='O3S')
    e_CO2, e_SO4, e_H2O, e_soot, e_NOx = climate.emissions(13300)
    climate.narrow_body_flights_per_region()
    # print('Asia_Narrow_body_flights',climate.Asia_Narrow_body_flights)
    # print('NA_Narrow_body_flights',climate.NA_Narrow_body_flights)
    # print('Europe_Narrow_body_flights',climate.Europe_Narrow_body_flights)
    # print('Latin_Narrow_body_flights',climate.Latin_Narrow_body_flights)
    # print('MiddleE_Narrow_body_flights',climate.MiddleE_Narrow_body_flights)
    # print('Africa_Narrow_body_flights',climate.Africa_Narrow_body_flights)
    H2_flights = climate.H2_percentage_trend()
    U_H2 = climate.number_aircraft_H2()                 # Number of flights in a year
    U_ker = U_H2
    #E_H2 = climate.E(e=e_CO2 ,U=U_H2,time=climate.t_prime)
    #E_ker = climate.E(e=e_CO2, U=U_ker, time=climate.t_prime)
    #percentages= climate.H2_percentage_trend()
    # plt.subplot(131)
    # plt.plot(climate.t_prime,percentages,label='Percentage of H2 flights')
    # plt.legend()
    # plt.subplot(121)
    # plt.plot(climate.t_prime,H2_flights, label='Percentage of H2 flights')
    # plt.xlabel('years')
    # plt.ylabel('Percentage [%]')
    # plt.legend()
    # plt.subplot(122)
    plt.plot(climate.t_prime,U_H2,label='Utilization rate',color='tab:red')
    plt.xlabel('years')
    plt.ylabel('Number of missions per year')
    plt.xlim(2035,2085)
    plt.legend()

    # plt.subplot(133)
    # plt.plot(climate.t_prime, E_H2,label='Emissions per year')
    # plt.legend()
    plt.show()
    # print('Start analysis for LTO')
    # 'To plot the change in CO2 concentration in ppmv per year'

    # clim =climate.delta_X_CO2(e=e_CO2,U=10,plot= True)
    # print(clim)
    #
    #
    # ATR = climate.ATR(h=11000,e_CO2= e_CO2, e_H2O= e_H2O, e_NOx=e_NOx,e_soot= e_soot,e_sulfate=e_SO4,U=U_ker,plot = True)
    # print('The average temperature response, A_100, for the LTO of the HACK is:',ATR,'[K]')
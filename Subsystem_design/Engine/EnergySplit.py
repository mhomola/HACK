""" Here we find the USABLE volume of each fuel type based on the energy they need to provide """

import numpy as np


class Energy_Split():
    def __init__(self):
        super().__init__()

        self.V_fuel_neo = 23859 # [L]

        self.rho_ker = 0.81 # [kg/L]
        self.rho_h2 = 0.07111 # [kg/L]
        self.LHV_ker = 43 # [MJ]
        self.LHV_h2 = 120 # [MJ]


        self.m_fuel_neo = self.rho_ker * self.V_fuel_neo # [kg]
        self.energy_tot_neo = self.m_fuel_neo * self.LHV_ker # [MJ]

        self.energy_h2 = 1/3 * self.energy_tot_neo # [MJ]
        self.energy_ker = self.energy_tot_neo - self.energy_h2 # [MJ]

        self.V_h2 = self.energy_h2 / self.LHV_h2 / self.rho_h2 # [L]
        self.V_ker = self.energy_ker / self.LHV_ker / self.rho_ker # [L]

        self.m_h2 = self.energy_h2 / self.LHV_h2 # [kg]
        self.m_ker = self.energy_ker / self.LHV_ker # [kg]


        self.m_h2_DPU = 250 # [kg]
        self.energy_h2_DPU = self.LHV_h2 * self.m_h2_DPU # [MJ]
        self.V_h2_DPU = self.m_h2_DPU / self.rho_h2

        # self.idle = 1229 # [s]
        self.t_taxi_out = 7.5 * 60 # [s]
        self.t_taxi_in = 8 * 60 # [s]
        self.fuelflow_neo_to = 0.373 # 0.250196 # [kg/s]
        self.fuelflow_neo_ti = 0.373 # 0.137205 # [kg/s]
        self.energy_h2_to = self.LHV_h2 * self.fuelflow_neo_to * self.t_taxi_out # [MJ]
        self.energy_h2_ti = self.LHV_h2 * self.fuelflow_neo_ti * self.t_taxi_in  # [MJ]
        self.m_h2_to = self.energy_h2_to / self.LHV_h2 # [kg]
        self.m_h2_ti = self.energy_h2_ti / self.LHV_h2  # [kg]


        # ======= ENERGY RATIO ======= #
        self.energy_h2_flight = self.energy_h2 - self.energy_h2_DPU - self.energy_h2_to - self.energy_h2_ti
        self.energy_ker_flight = self.energy_ker

        self.ER_h2_flight = self.energy_h2_flight / ( self.energy_h2_flight + self.energy_ker_flight )
        self.ER_ker_flight = self.energy_ker_flight / ( self.energy_h2_flight + self.energy_ker_flight )

        self.ER_h2 = np.array([ 1, self.ER_h2_flight, self.ER_h2_flight, self.ER_h2_flight, self.ER_h2_flight, 1, 1 ]) # idle, taxi out, take-off, climb, cruise, approach, taxi in
        self.ER_ker = np.array([ 0, self.ER_ker_flight, self.ER_ker_flight, self.ER_ker_flight, self.ER_ker_flight, 0, 0])  # idle, taxi out, take-off, climb, cruise, approach, taxi in

        # find LHV_f for each phase, according to mass fractions
        self.LHV_hack = self.ER_h2*self.LHV_h2 + self.ER_ker*self.LHV_ker  # [MJ/kg], array of LHV for each phase


        # ======= MASS RATIO ======= #
        self.m_h2_flight = self.m_h2 - self.m_h2_DPU - self.m_h2_to - self.m_h2_ti
        self.m_ker_flight = self.m_ker
        self.MR_h2_flight = self.m_h2_flight / ( self.m_h2_flight + self.m_ker_flight )
        self.MR_ker_flight = self.m_ker_flight / ( self.m_h2_flight + self.m_ker_flight )
        self.MR_h2 = np.array([ 1, self.MR_h2_flight, self.MR_h2_flight, self.MR_h2_flight, self.MR_h2_flight, 1, 1])  # idle, taxi out, take-off, climb, cruise, approach, taxi in
        self.MR_ker = np.array([0, self.MR_ker_flight, self.MR_ker_flight, self.MR_ker_flight, self.MR_ker_flight, 0, 0])  # idle, taxi out, take-off, climb, cruise, approach, ta


es = Energy_Split()
LHV_hack = es.LHV_hack
ER_h2, MR_h2 = es.ER_h2, es.MR_h2
ER_ker, MR_ker = es.ER_ker, es.MR_ker

if __name__ == '__main__':
    print('Phases: ', np.array(['taxi out', 'takeoff', 'climb', 'cruise', 'approach', 'taxi in', 'idle']))
    print('\nEnergy ratio of H2: ', es.ER_h2)
    print('Energy ratio of kerosene: ', es.ER_ker)
    print('\nMass ratio of H2: ', es.MR_h2)
    print('Mass ratio of kerosene: ', es.MR_ker)
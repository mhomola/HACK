from Subsystem_design.common_constants import Constants
import matplotlib.pyplot as plt
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
import numpy as np
from math import log
from Subsystem_design.Engine.Thrust_Required import thrust_req
from Subsystem_design.aerodynamic_subsys import cd0clean, wingar
from Subsystem_design.Engine.EnergySplit import Energy_Split

class Compute_weight(Constants):
    def __init__(self):
        super().__init__()


    def Tank_mass(self):
        self.tank_mass = 2*self.pod_tank_mass
    def Feeding_sys_m(self):
        self.Feeding_mass = 22                          # Feeding system from tanks to engine + tank to tank
    def Struc_m(self):
        self.struc_mass = self.Wing_Weight_320HACK-self.Wing_Weight_320neo + self.pylon_weight
    def Fuel_cell_m(self):
        self.FC_m = 781                                 # Entire system [kg]
    def weight_break_down_HACK(self):
        self.Tank_mass()
        self.Feeding_sys_m()
        self.Struc_m()
        self.Fuel_cell_m()

        self.OEW_HACK = self.OEW_320neo +self.struc_mass +self.Feeding_mass + self.FC_m + self.tank_mass

        self.Max_fuel_mass_capacity_HACK = self.W_kerosene + 2* self.pod_H2_mass

        self.MPLW_HACK = self.MPLW_320neo                           # Maximum Payload weight of HACK                  [kg]
        self.MZFW_HACK = self.MPLW_HACK + self.OEW_HACK             # Maximum zero fuel weight of HACK                [kg]
        self.MTOW_HACK = self.MTOW_320neo                           # Maximum take-off weight of Hack                 [kg]

        self.Max_fuel_at_max_PL_HACK = self.MTOW_HACK - self.MPLW_HACK - self.OEW_HACK  # Maximum fuel at maximum payload[kg]

        self.MRW_HACK = self.MRW_320HACK

class performance(Compute_weight):
    def flight_profile_weights(self,Mf,Mto):
        #Starts with MTOW, obtain Wf/Wto:
        Wf_Wto = Mf/Mto                                         # Fuel fraction between maximum fuel mass          [-]
                                                                # @ max. payload, and MTOW
        M_ff = 1 - Wf_Wto
        W5_W4 = M_ff  * (1/self.W4_W3) * (1/self.W6_W5) * (1/self.W7_W6)   # Fuel ratio between the end            [-]
                                                                           # and beginning of cruise
        return W5_W4
    def Range(self,L_D_ratio,cruise_f_ratio,SFC):

        R = (self.V_cruise/(SFC*self.g_0)) * (L_D_ratio) *log(1/cruise_f_ratio)
        return R
    def payload_range_diagram(self,L_over_D,SFC,mission,phase_durations):

        # Function to plot the payload-range diagram
        if mission == 'hack':
            self.weight_break_down_HACK()
            MTOW = self.MTOW_HACK
            MPLW = self.MPLW_HACK
            Max_fuel_MPLW = self.Max_fuel_at_max_PL_HACK
            Max_fuel_capacity = self.Max_fuel_mass_capacity_HACK
            OEW = self.OEW_HACK
        if mission == 'neo':
            MTOW = self.MTOW_320neo
            MPLW = self.MPLW_320neo
            Max_fuel_MPLW = self.MTOW_320neo - self.MPLW_320neo - self.OEW_320neo
            Max_fuel_capacity = self.Max_fuel_mass_capacity_320neo
            OEW = self.OEW_320neo

        #Point A: R = 0 and W_Pl = MPLW
        R_A = 0
        W_Pl_A = MPLW                                               # Payload mass at point A             [kg]

        #Point B: R = ?? and W_Pl = MPLW
        Mto = MTOW
        W_Pl_B = MPLW                                               # Payload mass at point B             [kg]
        Mf = Max_fuel_MPLW                                          # Maximum fuel at maximum payload    [kg]
        W5W4_B = self.flight_profile_weights(Mf,Mto)                          # Fuel ratio during cruise for point B [-]
        R_B = self.Range(L_over_D,W5W4_B,SFC)                                     # Range at point B                     [m]


        # Point C: R = ?? and Wf = MFW, Wpl = MTOW - OEM - Maxfuel
        if mission == 'hack':
            self.Energy_split(Max_fuel_capacity)                                          # Maximum fuel capacity of HACK      [kg]
            massH2, massker = self.get_fuel_before_take_off(self.mass_h2,self.mass_k,phase_durations,mission,True)
            Mf = massH2 + massker
        else:
            massH2, massker = self.get_fuel_before_take_off(0.,Max_fuel_capacity, phase_durations, mission, True)
            Mf = massH2 + massker

        W_Pl_C = Mto - Mf - OEW                                           #Payload mass at point C             [kg]
        W5W4_C = self.flight_profile_weights(Mf, Mto)                               # Fuel ratio during cruise for point C [-]
        R_C = self.Range(L_over_D, W5W4_C,SFC)                                      # Range at point C                     [m]


        # Point D: R = ?? and W_Pl = 0--> W_fuel = Maxfuel
        W_PL_D = 0                                                                  # Payload mass at point D            [kg]
        Mto = OEW + Mf + W_PL_D                                           # Take off weight at point D         [kg]
        W5W4_D = self.flight_profile_weights(Mf, Mto)                               # Fuel ratio during cruise for point D [-]
        R_D = self.Range(L_over_D, W5W4_D,SFC)                                      # Range at point D                     [m]

        self.Range_array = np.array([R_A*0.001,R_B*0.001,R_C*0.001,R_D*0.001])
        self.Payload_array = np.array([W_Pl_A,W_Pl_B,W_Pl_C,W_PL_D])

    def read_files(self,name_of_file):
        main_file = 'C:\\Users\\daf6111\\Documents\\universidade\\Third year\\DSE\\HACK\\Subsystem_design\\Engine'
        file = open(main_file + '\\' + name_of_file,'r')
        file_data = file.readlines()
        mf_h2 = float(file_data[20].split('\t')[1])
        mf_ker = float(file_data[21].split('\t')[1])
        return mf_h2,mf_ker

    def Energy_split(self,m_tot):
        E_k = (m_tot * self.LHV_ker * self.LHV_h2)/(self.LHV_h2+self.LHV_ker)
        E_h2 = 0.5 * E_k
        self.mass_k = E_k / self.LHV_ker
        self.mass_h2 = E_h2/ self.LHV_h2

    def get_fuel_before_take_off(self,m_h2,m_k,time_phases,mission,take_off):
        # Note: must read taxi-out and take-off files to get mass of fuel at beggining of mission.
        mf_h2_taxiout, mf_k_taxiout = self.read_files(mission + '_taxi_out.txt')
        mf_h2_takeoff1, mf_k_takeoff1 = self.read_files(mission + '_take_off1.txt')
        mf_h2_takeoff2, mf_k_takeoff2 = self.read_files(mission + '_take_off2.txt')
        if take_off == True:
            m_h2 = m_h2 - mf_h2_taxiout * time_phases[0] - mf_h2_takeoff1 * time_phases[1]- mf_h2_takeoff2 * time_phases[2]
            m_k = m_k - mf_k_taxiout * time_phases[0] - mf_k_takeoff1 * time_phases[1] - mf_k_takeoff2 * time_phases[2]
        else:
            m_h2 = m_h2 + mf_h2_taxiout * time_phases[0] + mf_h2_takeoff1 * time_phases[1] + mf_h2_takeoff2 * time_phases[2]
            m_k = m_k + mf_k_taxiout * time_phases[0] + mf_k_takeoff1 * time_phases[1] + mf_k_takeoff2 * time_phases[2]

        return m_h2,m_k

    def mission_profile(self,phase_durations,mission):

        if mission ==  'hack':
            # m_f is the fuel weigth available at Take-off
            self.weight_break_down_HACK()
            m_f = self.Max_fuel_at_max_PL_HACK
            self.Energy_split(m_tot=m_f)
            m_h2 = self.mass_h2
            m_k = self.mass_k

        if mission == 'neo':
            # m_f is the fuel weigth available at Take-off
            m_f = self.MTOW_320neo - self.MPLW_320neo - self.OEW_320neo
            m_h2 = 0
            m_k = m_f

        time_phases = phase_durations
        m_h2, m_k = self.get_fuel_before_take_off(m_h2,m_k,time_phases,mission,take_off = False)

        if mission == 'hack':
            m_f = m_h2 + m_k
            self.Energy_split(m_tot=m_f)
            m_h2 = self.mass_h2
            m_k = self.mass_k

        'Taxi-out'
        mf_h2_taxiout, mf_k_taxiout = self.read_files(mission + '_taxi_out.txt')
        fs_h2_taxi_out = np.linspace(m_h2, m_h2 - mf_h2_taxiout * time_phases[0], 50)
        fs_k_taxi_out = np.linspace(m_k, m_k - mf_k_taxiout * time_phases[0], 50)
        'Take-off1'
        mf_h2_takeoff1, mf_k_takeoff1 = self.read_files(mission + '_take_off1.txt')
        fs_h2_take_off1 = np.linspace(fs_h2_taxi_out[-1], fs_h2_taxi_out[-1] - mf_h2_takeoff1 * time_phases[1], 50)
        fs_k_take_off1 = np.linspace(fs_k_taxi_out[-1], fs_k_taxi_out[-1] - mf_k_takeoff1 * time_phases[1], 50)
        'Take-off2'
        mf_h2_takeoff2, mf_k_takeoff2 = self.read_files(mission + '_take_off2.txt')
        fs_h2_take_off2 = np.linspace(fs_h2_take_off1[-1], fs_h2_take_off1[-1] - mf_h2_takeoff2 * time_phases[2], 50)
        fs_k_take_off2 = np.linspace(fs_k_take_off1[-1], fs_k_take_off1[-1] - mf_k_takeoff2 * time_phases[2], 50)
        'Take-off3'
        mf_h2_takeoff3, mf_k_takeoff3 = self.read_files(mission + '_take_off3.txt')
        fs_h2_take_off3 = np.linspace(fs_h2_take_off2[-1], fs_h2_take_off2[-1] - mf_h2_takeoff3 * time_phases[3], 50)
        fs_k_take_off3 = np.linspace(fs_k_take_off2[-1], fs_k_take_off2[-1] - mf_k_takeoff3 * time_phases[3], 50)
        'Climb1'
        mf_h2_climb1, mf_k_climb1 = self.read_files(mission +'_climb1.txt')
        fs_h2_climb1 = np.linspace(fs_h2_take_off3[-1], fs_h2_take_off3[-1] - mf_h2_climb1 * time_phases[4], 100)
        fs_k_climb1 = np.linspace(fs_k_take_off3[-1], fs_k_take_off3[-1] - mf_k_climb1 * time_phases[4], 100)
        'Climb2'
        mf_h2_climb2, mf_k_climb2 = self.read_files(mission + '_climb2.txt')
        fs_h2_climb2 = np.linspace(fs_h2_climb1[-1], fs_h2_climb1[-1] - mf_h2_climb2 * time_phases[5], 100)
        fs_k_climb2 = np.linspace(fs_k_climb1[-1], fs_k_climb1[-1] - mf_k_climb2 * time_phases[5], 100)
        'Cruise1'
        mf_h2_cruise1, mf_k_cruise1 = self.read_files(mission +'_cruise1.txt')
        fs_h2_cruise1 = np.linspace(fs_h2_climb2[-1], fs_h2_climb2[-1] - mf_h2_cruise1 * time_phases[6],200)
        fs_k_cruise1 = np.linspace(fs_k_climb2[-1], fs_k_climb2[-1] - mf_k_cruise1 * time_phases[6], 200)
        'Cruise2'
        mf_h2_cruise2, mf_k_cruise2 = self.read_files(mission + '_cruise1.txt')
        fs_h2_cruise2 = np.linspace(fs_h2_cruise1[-1], fs_h2_cruise1[-1] - mf_h2_cruise2 * time_phases[7], 50)
        fs_k_cruise2 = np.linspace(fs_k_cruise1[-1], fs_k_cruise1[-1] - mf_k_cruise2 * time_phases[7], 50)
        'Approach1'
        mf_h2_approach1, mf_k_approach1 = self.read_files(mission +'_approach1.txt')
        fs_h2_approach1 = np.linspace(fs_h2_cruise2[-1], fs_h2_cruise2[-1] - mf_h2_approach1 * time_phases[8], 50)
        fs_k_approach1 = np.linspace(fs_k_cruise2[-1], fs_k_cruise2[-1] - mf_k_approach1 * time_phases[8], 50)
        'Approach2'
        mf_h2_approach2, mf_k_approach2 = self.read_files(mission + '_approach2.txt')
        fs_h2_approach2 = np.linspace(fs_h2_approach1[-1], fs_h2_approach1[-1] - mf_h2_approach2 * time_phases[9], 50)
        fs_k_approach2 = np.linspace(fs_k_approach1[-1], fs_k_approach1[-1] - mf_k_approach2 * time_phases[9], 50)
        'Taxi-in'
        mf_h2_taxiin, mf_k_taxiin = self.read_files(mission +'_taxi_in.txt')
        fs_h2_taxi_in = np.linspace(fs_h2_approach2[-1], fs_h2_approach2[-1] - mf_h2_taxiin * time_phases[10], 50)
        fs_k_taxi_in = np.linspace(fs_k_approach2[-1], fs_k_approach2[-1] - mf_k_taxiin * time_phases[10], 50)

        fs_h2_total = np.hstack((fs_h2_taxi_out,fs_h2_take_off1,fs_h2_take_off2,fs_h2_take_off3,fs_h2_climb1,fs_h2_climb2,fs_h2_cruise1,fs_h2_cruise2,fs_h2_approach1,fs_h2_approach2,fs_h2_taxi_in))
        fs_k_total = np.hstack((fs_k_taxi_out,fs_k_take_off1,fs_k_take_off2,fs_k_take_off3,fs_k_climb1,fs_k_climb2,fs_k_cruise1,fs_k_cruise2,fs_k_approach1,fs_k_approach2,fs_k_taxi_in))

        return fs_h2_total,fs_k_total


if __name__ == '__main__':
    const = Constants()
    """Weights of A320-HACK"""

    AC_weights = Compute_weight()                                               # Initiallize class of weight estimation
    AC_weights.weight_break_down_HACK()

    Aerodynamic_charac = AerodynamicCharacteristics()
    Aerodynamic_charac.L_over_D_cruise()

    Performance = performance()

    T = thrust_req(cd0clean, wingar)

    '''Plotting fuel consumption during flight'''
    fs_h2_total_hack,fs_k_total_hack = Performance.mission_profile(phase_durations= T.durations, mission='hack')
    fs_h2_total_neo, fs_k_total_neo = Performance.mission_profile(phase_durations= T.durations, mission='neo')
    time = np.linspace(0,np.sum(T.durations),800)
    plt.plot(time,fs_h2_total_hack,label= 'Hydrogen')
    plt.plot(time, fs_k_total_hack, label='Kerosene')
    plt.plot(time, fs_h2_total_hack+fs_k_total_hack, label='Total')
    plt.ylabel('Mass of fuel on board [kg]',fontsize = 15)
    plt.xlabel('Time [s]',fontsize = 15)
    plt.legend(fontsize = 15)
    plt.show()

    '''Plotting payload range diagrams'''
    #HACK
    main_file = 'C:\\Users\\daf6111\\Documents\\universidade\\Third year\\DSE\\HACK\\Subsystem_design\\Engine'
    file1 = open(main_file + '\\' + 'hack_cruise.txt', 'r')
    file_data1 = file1.readlines()
    TSFC_cruise1 = float(file_data1[39].split('\t')[1])

    #NEO
    file2 = open(main_file + '\\' + 'neo_cruise.txt', 'r')
    file_data2 = file2.readlines()
    TSFC_cruise2 = float(file_data2[39].split('\t')[1])

    Performance.payload_range_diagram(L_over_D=Aerodynamic_charac.L_D_ratio_HACK,SFC= TSFC_cruise1*10**-6,mission='hack',phase_durations=T.durations)
    Range_HACK = Performance.Range_array
    Payload_HACK = Performance.Payload_array

    Performance.payload_range_diagram(L_over_D=Aerodynamic_charac.L_D_ratio_neo,SFC = TSFC_cruise2*10**-6,mission='neo',phase_durations=T.durations)
    Range_neo = Performance.Range_array
    Payload_neo = Performance.Payload_array

    plt.plot(Range_HACK, Payload_HACK, marker='*', color='tab:red',label='A320-HACK')
    plt.plot(Range_neo,Payload_neo,marker='o',color='navy',label='A320neo')
    plt.xlabel('Range [km]')
    plt.ylabel('Payload Mass [kg]')
    plt.legend()
    plt.show()



    print('The OEW of the A320HACK is:',AC_weights.OEW_HACK)
    print('The MPLW of the A320HACK is:', AC_weights.MPLW_HACK)
    print('The max fuel mass of the A320HACK is:', AC_weights.Max_fuel_mass_capacity_HACK)
    print('The max fuel @ MPLW of the A320HACK is:', AC_weights.Max_fuel_at_max_PL_HACK)
    print('The MZFW of the A320HACK is:', AC_weights.MZFW_HACK)
    print('The MTOW of the A320HACK is:', AC_weights.MTOW_HACK)
    print('The extra structural mass of the A320HACK is:', AC_weights.struc_mass)
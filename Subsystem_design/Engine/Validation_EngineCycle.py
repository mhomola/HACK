from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.EngineCycle import Engine_Cycle
# from Subsystem_design.Engine.CoolEngine import *
import numpy as np

class Valid_Engine(Constants):
    def __init__(self):
        super().__init__()

    def true_vals_cruise(self, p0, W_HPC, mf_airfuel):
        self.v0_real = 230.17168   # [m/s]
        self.T00_real = 243.01197 # [K]
        self.p00_real = 0.32366 * 10**5 # [Pa]
        self.T02_real = 243.01197 # [Pa]
        self.p02_real = 0.31395 * 10**5 # [Pa]
        self.T021_real = 271.40302 # [K]
        self.p021_real = 0.446 * 10**5 # [Pa]
        self.T013_real = 271.40125 # [K]
        self.p013_real = 0.44599 * 10**5 # [Pa]
        self.T025_real = 369.93091 # [K]
        self.p025_real = 1.2016 * 10**5 # [Pa]
        self.T03_real = 724.90893 # [K]
        self.p03_real = 11.70102 * 10**5 # [Pa]
        self.T04_real = 1459.30433 # [K]
        self.p04_real = 10.99347 * 10**5 # [Pa]
        self.T045_real = 1016.76997 # [K] named TT49
        self.p045_real = 2.87838 * 10**5 # [Pa] named PT49
        self.T05_real = 625.36361 # [K] named TT6 (there's a T54)
        self.p05_real = 0.35812 * 10**5 # [Pa] named PT6 (there's a PT54)
        self.T8_real = 546.56551 # [K]
        self.p8_real = 0.21654 * 10**5 # [Pa]
        self.v8_real = 411.09753 # [m/s]
        self.T18_real = 226.07871 # [K]
        self.p18_real = 0.23257 * 10**5 # [Pa]
        self.v18_real = 301.54515 # [m/s]

        self.OPR_real = 37.27031

        self.eta_inlet_real = ( ( (self.p02_real/p0)**( (self.k_air-1)/self.k_air ) ) - 1 ) / ( (self.k_air-1)/2*(self.M0[i]**2) )
        self.W_HPT_real = ( self.T04_real - self.T045_real ) * mf_airfuel * self.cp_gas
        self.eta_mech_H_real = W_HPC / self.W_HPT_real
        self.PR_cc_real = self.p04_real/self.p03_real

        self.eta_HPT_real = ( 1 - self.T045_real/self.T04_real ) / ( 1 - (self.p045_real/self.p04_real) ** ((self.k_gas-1)/self.k_gas) )
        self.PR_nozzle_real = self.p05_real/self.p8_real
        self.eta_nozzle_real = (self.k_gas-1)/(self.k_gas+1)  *  1/( 1 - (1/self.PR_nozzle_real)**((self.k_gas-1)/self.k_gas) )

        # self.eta_nozzle_real = (1 - self.T8_real/self.T05_real) / ( 1 - (self.p8_real/self.p05_real)**((self.k_gas-1)/self.k_gas) )


if __name__ == "__main__":
    cycle, const, val = Engine_Cycle(), Constants(), Valid_Engine()
    index = -3
    print('Phase:',const.phases[index])
    cycle.cycle_analysis(aircraft='neo', i=index)
    e = 0.1
    i = index
    val.true_vals_cruise(const.p0[i], cycle.W_HPC, cycle.mf_airfuel)

    print( 'Error in v0: ', (const.v0[i] - val.v0_real)/val.v0_real )
    print('Error in T00: ', (cycle.T00 - val.T00_real) / val.T00_real)
    print('Error in p00: ', (cycle.p00 - val.p00_real) / val.p00_real)
    print('Error in T02: ', (cycle.T02 - val.T02_real) / val.T02_real)
    print('Error in p02: ', (cycle.p02 - val.p02_real) / val.p02_real)
    print('Error in T021: ', (cycle.T021 - val.T021_real) / val.T021_real)
    print('Error in p021: ', (cycle.p021 - val.p021_real) / val.p021_real)
    # print('Error in T013: ', (cycle.T013 - val.T013_real) / val.T013_real)
    # print('Error in p013: ', (cycle.p013 - val.p013_real) / val.p013_real)
    print('Error in T025: ', (cycle.T025 - val.T025_real) / val.T025_real)
    print('Error in p025: ', (cycle.p025 - val.p025_real) / val.p025_real)
    print('Error in T03: ', (cycle.T03 - val.T03_real) / val.T03_real)
    print('Error in p03: ', (cycle.p03 - val.p03_real) / val.p03_real)
    print('Error in T04: ', (cycle.T04 - val.T04_real) / val.T04_real)
    print('Error in p04: ', (cycle.p04 - val.p04_real) / val.p04_real)
    print('Error in T045: ', (cycle.T045 - val.T045_real) / val.T045_real)
    print('Error in p045: ', (cycle.p045 - val.p045_real) / val.p045_real)
    print('Error in T05: ', (cycle.T05 - val.T05_real) / val.T05_real)
    print('Error in p05: ', (cycle.p05 - val.p05_real) / val.p05_real)
    print('Error in T8: ', (cycle.T8 - val.T8_real) / val.T8_real)
    print('Error in p8: ', (cycle.p8 - val.p8_real) / val.p8_real)
    print('Error in v8: ', (cycle.v8 - val.v8_real) / val.v8_real)
    print('Error in T18: ', (cycle.T18 - val.T18_real) / val.T18_real)
    print('Error in p18: ', (cycle.p18 - val.p18_real) / val.p18_real)
    print('Error in v18: ', (cycle.v18 - val.v18_real) / val.v18_real)
    print('Error in OPR: ', (cycle.OPR - val.OPR_real) / val.OPR_real)

    print('\nREAL VALUES\nReal eta_inlet: ', val.eta_inlet_real)
    print('Real W_HPT:', val.W_HPT_real, 'Computed W_HPT:', cycle.W_HPT)
    print('Real eta_mech HPT:', val.eta_mech_H_real)
    print('Real PR_cc:', val.PR_cc_real, 'Real eta_HPT:', val.eta_HPT_real)
    print('Real eta_nozzle:', val.eta_nozzle_real, 'PR nozzle real:', val.PR_nozzle_real)



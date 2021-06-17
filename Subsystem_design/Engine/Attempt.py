from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.EngineCycle import Engine_Cycle
from Subsystem_design.Engine.CoolEngine import Engine_Cool
import numpy as np

cycle = Engine_Cycle()
const = Constants()
cool = Engine_Cool()

cycle.cycle_analysis('neo', 'take_off')

# FIND MF_FUEL
cool.integral(cool.N2_cp_data, cycle.T03, cycle.T04)
# A = cool.cp_integral * cycle.mf_hot
A = 1150 * (cycle.T04 - cycle.T03) * cycle.mf_hot
cool.cp_first_last(cool.N2_cp_data, cycle.T04)
cpf = cool.cp
print(cpf)
mf_fuel = A / (cycle.eta_cc * cycle.LHV_f * 10**6 - cycle.T04 * 1150 ) #cpf)
print(cycle.mf_fuel, 'vs.', mf_fuel)
# mf_fuel = cycle.mf_fuel

TPZ = cycle.T04 + 1
err = 20
save_err, save_TPZ, save_avg = np.array([]), np.array([]), np.array([])

while abs(err) > 0.05 and TPZ < 3000:
    print('TPZ = ', TPZ)

    cool.integral(cool.N2_cp_data, cycle.T03, cycle.T04)
    # L1 = cool.cp_integral * cycle.mf_hot
    L1 = 1150 * (cycle.T04 - cycle.T03) * cycle.mf_hot

    cool.integral(cool.C12H26_cp_data, cycle.T04, TPZ)
    # L2 = cool.cp_integral * mf_fuel
    L2 = 1150 * (TPZ - cycle.T04) * mf_fuel

    cool.integral(cool.N2_cp_data, cycle.T04, TPZ)
    # L3 = cool.cp_integral
    L3 = 1150 * (TPZ - cycle.T04)

    cool.integral(cool.N2_cp_data, cycle.T03, cycle.T04)
    # L4 = cool.cp_integral
    L4 = 1000 * (cycle.T04 - cycle.T03)

    R5 = cycle.eta_cc * mf_fuel * cycle.LHV_f * 10**6

    cool.cp_first_last(cool.C12H26_cp_data, TPZ)
    print(cool.cp)
    # R6 = cool.cp * mf_fuel * TPZ
    # R6 = 1150 * mf_fuel * TPZ
    # cool.integral(cool.C12H26_cp_data, cycle.T03, TPZ)
    # R6 = cool.cp_integral * mf_fuel
    R6 = 1150 * (TPZ - cycle.T03) * mf_fuel

    cool.integral(cool.N2_cp_data, cycle.T03, TPZ)
    # R7 = cool.cp_integral * mf_fuel
    R7 = 1150 * (TPZ - cycle.T03) * mf_fuel

    lhs = (L1 - L2) / (L3 + L4)
    rhs = (R5 - R6) / R7

    err = (lhs - rhs) * 2 / (lhs + rhs)
    print('LHS:', lhs, '; RHS:', rhs, '; err:', err)
    print('Mass flow of air for combustion:', lhs / cycle.mf_hot, 'vs.', rhs / cycle.mf_hot)
    save_err = np.append(save_err, abs(err))
    save_TPZ = np.append(save_TPZ, TPZ)
    save_avg = np.append(save_avg, (lhs + lhs) / 2)

    TPZ += 10



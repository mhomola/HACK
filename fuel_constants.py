def gas_law(p, T):
    return p*M_H2/(R*T)

H2_ed = 33.5                                      # Energy density of hydrogen                         [kWh/kg]
k_ed = 12.0                                       # Energy density of kerosene                         [kWh/kg]

# Assumption: storing H2 at 3 bar
p_tank = 300000                                   # Pressure in tank for liquid H2                     [Pa]
p_tank_g = 35000000                               # Pressure in tank for gas H2                        [Pa]
T_tank = 20                                       # Temperature in tank for liquid H2                  [K]
T_tank_g = 293.15                                   # Assuming H2 is stored at 20 C for gas              [K]
LH2_d = 65.0                                      # Mass density of liquid hydrogen                    [kg/m^3]
# H2_d = 38.0                                     # Mass density of gaseous hydrogen stored at 70 MPa  [kg/m^3]
M_H2 = 2.016*10**(-3)                             # Molar mass of H2                                   [kg/mole]
R = 8.314                                         # Gas constant                                       [J/(mol*K)]
GH2_d = gas_law(p_tank, T_tank)                   # Density of gas H2 (boil off caused)
GH2_d_g = gas_law(p_tank, T_tank)                 # Density of gas H2 (fully gaseous tanks)
k_d = 810.0                                       # Mass density of kerosene                           [kg/m^3]
y = 0.95                                          # LH2 volume ratio (from Winnefeld et al.)           [-]
fuel_capacity_a320neo = 26730                     # Maximum Fuel capacity of A320neo                   [l]


#Engine performance Constants (estimate)

eta_th = 0.55                                     # Engine thermal efficiency                          [-]
eta_prop = 0.42                                   # Engine propulsive efficiency                       [-]
eta_comb = 0.99                                   # Engine combustion efficiency                       [-]
eta_tot = eta_th * eta_prop * eta_comb            # Engine total efficiency                            [-]

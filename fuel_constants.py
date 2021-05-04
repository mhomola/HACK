def gas_law(p,T):
    return p*M_H2/(R*T)

H2_ed = 33.5        # Energy density of hydrogen                           [kWh/kg]
k_ed = 12.0         # Energy density of kerosene                           [kWh/kg]

# Assumption: storing H2 at 3 bar
p_tank = 300000 # [Pa]
T_tank = 20 # [K]
LH2_d = 65.0        # Mass density of liquid hydrogen                      [kg/m^3]
# H2_d = 38.0         # Mass density of gaseous hydrogen stored at 70 MPa    [kg/m^3]
M_H2 = 2.016*10**(-3) # Molar mass of H2                                    [kg/mole]
R = 8.314           # Gas constant                                          [J/(mol*K)]
GH2_d = gas_law(p_tank,T_tank)
k_d = 810.0         # Mass density of kerosene                             [kg/m^3]
y = 0.95            # LH2 volume ratio (from Winnefeld et al.)             [-]
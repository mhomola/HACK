import matlab.engine
eng = matlab.engine.start_matlab()
from EngineCycle import Engine_Cycle

ec = Engine_Cycle()
ec.cycle_analysis(aircraft='neo', i=-3)
print(ec.T03, ec.p03)

P_input = round(ec.p03,3)
T_input = round(ec.T03,3)

TPZ = eng.reactor1('kerosene', float(ec.p03), float(ec.T03), 1)
# TPZ = eng.reactor1('kerosene', 3283120.11, 811.11, 1)
print(TPZ - 10000)


# tri = eng.triarea(1.5,3)
# print(tri)
# print(tri*4)





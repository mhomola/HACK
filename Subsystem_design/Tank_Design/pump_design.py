import numpy as np
import matplotlib.pyplot as plt
from Fuel_Masses_Estimation import kff, h2ff

print('\n \n ---> PUMP DESIGN <--- ')

# Fuel flow of hydrogen into the engine:

F_flow_H2 = [h2ff[0], h2ff[2], h2ff[4], h2ff[6], h2ff[8], h2ff[10]]
flight_phases = ['Idle', 'Taxi out', 'Climb', 'Cruise', 'Descent', 'Taxi in']

print('\n The fuel flows during: ', flight_phases,
      '\n are                  : ', F_flow_H2)
plt.bar(flight_phases, F_flow_H2)
plt.show()

max_ff = max(F_flow_H2)

"""
So mass flow is constant and equal to m_ff = rho * A * V
I can therefore compute speed that the pump needs to provide. 
Also, the pump will have to provide a certain pressure for the hydrogen to get into the combustion chamber,
therefore I get the value of the pressure at the combustion chamber, I add the pressure loss due to the piping
and then I can know the pressure which the pump needs to deliver.
"""

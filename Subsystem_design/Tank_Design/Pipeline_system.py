import numpy as np
import math as m
path1 = (0.55-0.34)*17 #path from wing pod to engine
path2 = 0.34*17 +  6   #path from central tank to the engine
L = (path1+path2)*1.2  #20% to account for curves and all that

t_steel = 0.041/100  #[m]
t_al = 0.041/100     #[m]
t_foam = 3.81/100    #[m]

rho_steel = 8090 #kg/m^3 https://www.aksteel.com/sites/default/files/2018-01/321201512.pdf
rho_foam = 28.03 #kg/m^3 https://tigerfoam.com/sprayfoaminsulation/open-cell-vs-closed-cell-foam-which-should-i-choose/#:~:text=Closed%20cell%20foam%20is%20much%20denser%20than%20open%20cell%20foam,per%20cubic%20foot%20or%20more.
rho_al = 2720 #kg/m^3 https://www.engineersedge.com/materials/densities_of_metals_and_elements_table_13976.htm
d = 1.27/100 # [m]

V_steel = (m.pi * (d/2 + t_steel)**2 - m.pi*(d/2)**2)*L #[m^3]
V_foam = (m.pi * (d/2 + t_steel + t_foam )**2 - m.pi*(d/2 + t_steel)**2)*L #[m^3]
V_al = (m.pi * (d/2 + t_steel + t_foam +t_al)**2 - m.pi*(d/2 + t_steel + t_foam)**2)*L #[m^3]

m_steel = rho_steel* V_steel
m_foam = rho_foam * V_foam
m_al = rho_al * V_al

m_total = m_steel + m_foam + m_al

print(m_steel)
print(m_foam)
print(m_al)
print(m_total)
print(L)
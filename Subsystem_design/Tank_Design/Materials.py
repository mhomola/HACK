
class Material():
    def __init__(self, E, yield_strength, ult_strength, density, conductivity):
        self.E = E                            # GPa
        self.yield_strength = yield_strength  # MPa
        self.ult_strength = ult_strength      # MPa
        self.density = density                # kg / m^3
        self.conductivity = conductivity      # W / (m·K)


# Define the materials you will be using
# Al_8090_T8151 = Material(E=88, yield_strength=453, ult_strength=756, density=2540, conductivity=95.3)
Al_7075_T6 = Material(E=68.9, yield_strength=324, ult_strength=414, density=2700, conductivity=167)
Al_2219_T37 = Material(E=79, yield_strength=463, ult_strength=668, density=2840, conductivity=121)
Al_2090_T81 = Material(E=89, yield_strength=663, ult_strength=764, density=2590, conductivity=88)
WL049_T851 = Material(E=87, yield_strength=715, ult_strength=853, density=2540, conductivity=95.3) # will not be chosen
TI = Material(E=120, yield_strength=711, ult_strength=784, density=4490, conductivity=7.8)

# due to lack of information


class Insulation_Material():
    def __init__(self, conductivity, diffusivity, density):
        self.conductivity = conductivity  # W/(m·K)
        self.diffusivity = diffusivity    # m^2/s
        self.density = density            # kg/m^3


# Preliminary values, change those values for more specific

MLI = Insulation_Material(conductivity=3*10**-5, diffusivity=10**-4, density=35.24)

Rohacell_Foam_20K = Insulation_Material(conductivity=5*10**-3, diffusivity=10**-6, density=51.1)    # ciet
Rohacell_Foam_315K = Insulation_Material(conductivity=35*10**-3, diffusivity=10**-6, density=51.1)  # at 315K



class Material():
    def __init__(self, E, yield_strength, ult_strength, density, conductivity):
        self.E = E                            # GPa
        self.yield_strength = yield_strength  # MPa
        self.ult_strength = ult_strength      # MPa
        self.density = density                # kg / m^3
        self.conductivity = conductivity      # W / (m·K)


# Define the materials you will be using
# Al_7068 = Material(E=73, yield_strength=683, ult_strength=710, density=2850, conductivity=190)
Al_6061_T6 = Material(E=68.9, yield_strength=324, ult_strength=414, density=2700, conductivity=167)
Al_2219_T87 = Material(E=73.1, yield_strength=393, ult_strength=476, density=2840, conductivity=121)
Al_2090_T86 = Material(E=76, yield_strength=520, ult_strength=550, density=2590, conductivity=88)
Al_8090_T8151 = Material(E=77, yield_strength=450, ult_strength=370, density=2540, conductivity=95.3)

class Insulation_Material():
    def __init__(self, conductivity, diffusivity, density):
        self.conductivity = conductivity  # W/(m·K)
        self.diffusivity = diffusivity    # m^2/s
        self.density = density            # kg/m^3


# Preliminary values, change those values for more specific
MLI = Insulation_Material(conductivity=10**-3, diffusivity=10**-4, density=35.24)
Rohacell_Foam_best = Insulation_Material(conductivity=5*10**-3, diffusivity=10**-6, density=35.24)
Rohacell_Foam_worst = Insulation_Material(conductivity=35*10**-3, diffusivity=10**-6, density=35.24)
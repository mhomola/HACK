class Material():
    def __init__(self, E, yield_strength, ult_strength, density, conductivity):
        self.E = E                            # GPa
        self.yield_strength = yield_strength  # MPa
        self.ult_strength = ult_strength      # MPa
        self.density = density                # kg / m^3
        self.conductivity = conductivity      # W / (m·K)


# Define the materials you will be using
Al_2090_T81 = Material(E=89, yield_strength=663*0.7929837, ult_strength=764, density=2590, conductivity=88)


class Insulation_Material():
    def __init__(self, conductivity, diffusivity, density):
        self.conductivity = conductivity  # W/(m·K)
        self.diffusivity = diffusivity    # m^2/s
        self.density = density            # kg/m^3


# Preliminary values, change those values for more specific
MLI = Insulation_Material(conductivity=3*10**-5, diffusivity=10**-4, density=50)




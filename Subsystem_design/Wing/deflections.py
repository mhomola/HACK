from Subsystem_design.Wing.loads import Loads_w
from Subsystem_design.common_constants import Constants
import scipy.integrate as spint
from Subsystem_design.Wing.inertia import Inertia_normal, Inertia_shear

class Deflections(Constants):

    def __init__(self, E, Ixx, Iyy):
        super(Deflections, self).__init__()
        self.E = E
        self.Ixx = Ixx
        self.Iyy = Iyy

    def vertical_deflection(self, z):
        I = Inertia_normal(n_str=11)

        Ixx =
        v_slope = -1/self.E * spint.simps(y=Mx/Ixx, x=z)

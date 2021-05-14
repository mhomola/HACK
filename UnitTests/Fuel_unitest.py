import unittest
import fuel_constants
from Subsystem_design import fuel_required


class TestConstant(unittest.TestCase):

    def test_H2_ed_constant(self):
        """
        This function investigates the change in the calculated volume if the energy density of the H2 is changed to the liquid value.
        """
        V_k, V_H2 = fuel_required.fuel_volume_calc(d_LH2=fuel_constants.LH2_d, d_GH2=fuel_constants.GH2_d, d_GH2_g=fuel_constants.GH2_d_g, d_k=fuel_constants.k_d,
                                                   Ed_H2=fuel_constants.H2_ed, Ed_k=fuel_constants.k_ed,
                                                   tot_vol_k=fuel_constants.fuel_capacity_a320neo,
                                                   e_ratio=fuel_required.E_ratio, state='liquid')
        H2_ed = 37.5
        V_k, V_H2_l = fuel_required.fuel_volume_calc(d_LH2=fuel_constants.LH2_d, d_GH2=fuel_constants.GH2_d,
                                                     d_GH2_g=fuel_constants.GH2_d_g, d_k=fuel_constants.k_d,
                                                     Ed_H2=H2_ed, Ed_k=fuel_constants.k_ed,
                                                     tot_vol_k=fuel_constants.fuel_capacity_a320neo,
                                                     e_ratio=fuel_required.E_ratio, state='liquid')
        print(V_H2) #volume with energy desnity constant of gas H2
        print(V_H2_l) #volume with energy desnity constant of LH2
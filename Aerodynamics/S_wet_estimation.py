"""
This file contains some helper functions for a quick estimation of the surface wetted area of the design options.
"""

import math as m
import numpy as np

class S_wet_estimation_standard():
    def __init__(self, l_cockpit, l_cabin, l_tail, df1, df2):
        """

        :param l_cockpit: cockpit length in [m]
        :param l_cabin: cabin length in [m]
        :param l_tail: tail length in [m]
        :param df1: diameter 1 fuselage-vertical [m]
        :param df2: diameter 2 fuselage-horizontal [m]
        """
        self.l_cockpit = l_cockpit
        self.l_cabin = l_cabin
        self.l_tail = l_tail
        self.df1 = df1
        self.df2 = df2

    def cockpit_volume(self):
        """
        The cockpit is modelled as half an ellipsoid cut along the vertical plane.
        :return: Cockpit volume in [m^3]
        """
        self.V_cockpit = (4/3*m.pi*(self.l_cockpit/2)*(self.df1/2)*(self.df2/2))/2

    def cabin_volume(self):
        """
        The cabin is modelled as a cylinder with ellipse base
        :return: Cabin volume in [m^3]
        """
        self.V_cabin = m.pi * (self.df1/2)*(self.df2/2) * self.l_cabin

    def tail_volume(self):
        """
        The tail is modelled as a cone
        :return: Tail volume in [m^3]
        """
        self.V_tail = 1.3*m.pi*(self.df1/2)*(self.df2/2)*self.l_tail

    def calculate_volume(self):
        self.cockpit_volume()
        self.cabin_volume()
        self.tail_volume()
        self.volume = self.V_cockpit + self.V_cabin + self.V_tail

    def S_wet(self):
        """
        We use linear regression of the graph which shows the relation between the volume of the aircraft and the wet surface area.
        The relation uses feet dimensions so we need to convert.
        """
        ft = 0.3048 #1 ft in [m]
        self.S_wet_fus = 13.6 * (self.volume/ft**3)**0.668 * ft**2 # in [m2]


class S_wet_estimation_belly():
    def __init__(self,l_cockpit,l_cabin,l_tail,df1,dh,rb):
        """

        :param l_cockpit: cockpit length in [m]
        :param l_cabin: cabin length in [m]
        :param l_tail: tail length in [m]
        :param df1: diameter 1 fuselage-vertical [m]
        :param dh: distance between belly and lower part of circular fuselage [m]
        :param rb: the radius of the curved part at the belly [m]
        """
        self.l_cockpit = l_cockpit
        self.l_cabin = l_cabin
        self.l_tail = l_tail
        self.df1 = df1
        self.dh = dh
        self.rb = rb

    def cockpit_volume(self):
        """
        The cockpit is modelled as half an ellipsoid cut along the vertical plane.
        :return: Cockpit volume in [m^3]
        """
        self.V_cockpit = (4/3*m.pi*(self.l_cockpit/2)*(self.df1/2)**2)/2

    def cabin_volume(self):
        """
        The cabin is modelled as a cylinder with ellipse base
        :return: Cabin volume in [m^3]
        """
        self.cross_section = 0.5 * np.pi * (self.df1/2)**2 + 2 * (self.df1/2 * (self.df1/2 + self.dh) -
                                                                  self.rb**2 + 0.25 * np.pi * self.rb**2)
        self.V_cabin = self.cross_section * self.l_cabin

    def tail_volume(self):
        """
        The tail is modelled as a cone
        :return: Tail volume in [m^3]
        """
        self.V_tail = 1.3*m.pi*(self.df1/2)**2*self.l_tail

    def calculate_volume(self):
        self.cockpit_volume()
        self.cabin_volume()
        self.tail_volume()
        self.volume = self.V_cockpit + self.V_cabin + self.V_tail

    def S_wet(self):
        """
        We use linear regression of the graph which shows the relation between the volume of the aircraft and the wet surface area.
        The relation uses feet dimensions so we need to convert.
        """
        ft = 0.3048 #1 ft in [m]
        self.S_wet_fus = 13.6 * (self.volume/ft**3)**0.668 * ft**2 # in [m2]

class S_wet_estimation_beluga():
    def __init__(self, l_cockpit, l_cabin, l_tail, beluga, df, dfb):
        """

        :param l_cockpit: cockpit length in [m]
        :param l_cabin: cabin length in [m]
        :param l_tail: tail length in [m]
        :param beluga: beluga percentage of cabin length[%]
        :param df: diameter fuselage [m]
        :param dfb: diameter beluga[m]
        """
        self.l_cockpit = l_cockpit
        self.l_cabin = l_cabin
        self.l_tail = l_tail
        self.beluga = beluga
        self.l_beluga = self.l_cabin*beluga/100
        self.df = df
        self.dfb = dfb

    def beluga_area(self):
        """
        :return: Cross-sectional area of Beluga section.
        """
        factor = 0.6 #how much of the upper beluga circle is not overlapped with the main fuselage
        self.S_beluga = m.pi * (self.df/2)**2 + m.pi * (self.dfb/2)**2 * factor

    def cockpit_volume(self):
        """
        The cockpit is modelled as half an ellipsoid cut along the vertical plane. The base of the ellispoid is assumed
        to be a circle whose readius results in the surface of the beluga cross-section.
        :return: Cockpit volume in [m^3]
        """
        self.V_cockpit = (4/3*m.pi*(self.l_cockpit/2)*self.S_beluga/m.pi)/2

    def beluga_cabin_volume(self):
        """
        The beluga cabin is modelled as a cylinder whose base is equal to the surface of the beluga cross-section.
        :return: beluga cabin volume in [m^3]
        """
        self.V_beluga_cabin = self.S_beluga * self.l_beluga

    def cabin_volume(self):
        """
        The cabin is modelled as a cylinder with circular base df.
        :return: Cabin volume in [m^3]
        """
        self.V_cabin = m.pi * (self.df/2)**2 * (self.l_cabin-self.l_beluga)

    def tail_volume(self):
        """
        The tail is modelled as a cone with circular base.
        :return: Tail volume in [m^3]
        """
        self.V_tail = 1.3*m.pi*(self.df/2)**2*self.l_tail

    def calculate_volume(self):
        self.beluga_area()
        self.cockpit_volume()
        self.beluga_cabin_volume()
        self.cabin_volume()
        self.tail_volume()
        self.volume = self.V_cockpit + self.V_beluga_cabin + self.V_cabin + self.V_tail

    def S_wet(self):
        """
        We use linear regression of the graph which shows the relation between the volume of the aircraft and the wet surface area.
        The relation uses feet dimensions so we need to convert.
        """
        ft = 0.3048 #1 ft in [m]
        self.S_wet_fus = 13.6 * (self.volume/ft**3)**0.668 * ft**2 # in [m2]

if __name__ == '__main__':
    #trial for a normal configuration
    design1 = S_wet_estimation_standard(l_cockpit=5.04,l_cabin=24.49,l_tail=8.04,df1=4.14,df2=4.14)
    design1.calculate_volume()
    print(design1.volume)
    design1.S_wet()
    design2 = S_wet_estimation_beluga(l_cockpit=5.04,l_cabin=24.49,l_tail=8.04,beluga=100,df=4.14,dfb=1.968)
    design2.calculate_volume()
    print(design2.volume)
    design2.S_wet()
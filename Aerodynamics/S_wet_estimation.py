"""
This file contains some helper functions for a quick estimation of the surface wetted area of the design options.
"""

import math as m
import numpy as np

class S_wet_estimation():
    def __init__(self,l_cockpit,l_cabin,l_tail,df1,df2):
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

if __name__ == '__main__':
    #trial for a normal configuration
    #if you have additional volumes such as a beluga put them as an add on to the volume attribute after running th calculate_volume method
    design1 = S_wet_estimation(l_cockpit=5.04,l_cabin=24.49,l_tail=8.04,df1=4.14,df2=4.14)
    design1.calculate_volume()
    print(design1.volume)
    design1.S_wet()


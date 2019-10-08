import numpy as np 

class Receiver:
    def __init__(self, x,y,z):
        self.x = x
        self.y = y
        self.z = z
        self.location= np.array([self.x, self.y, self.z])
        self.type = "FOV"
        self.TYPES = ("FOV", "ADR", "IM")
        self.area = 1 #cm2
        self.Concentrator_gain = 1
        self.Filter_gain = 1

    def calc_power_at_receiver(self):
        pass

    def set_FOV(self, FOV):
        self.FOV = FOV

    def set_ADR_parameters(self, Azimuth_angles, Elevation_angles):
        self.Azimuth_angles = Azimuth_angles
        self.Elevation_angles = Elevation_angles

    def set_reception_area(self, area):
        self.area = area

    def set_type(self, type):
        if(type in self.TYPES):
            print("exists in types ")
            self.type = type

    def change_height(self, z):
        self.z = z
        self.location= np.array([self.x, self.y, self.z])
    def change_x(self, x):
        self.x = x
        self.location= np.array([self.x, self.y, self.z])
    def change_y(self, y):
        self.y = y
        self.location= np.array([self.x, self.y, self.z])


    def __repr__(self):
        return "Receiver at location {} with type {}".format(self.location, self.type)

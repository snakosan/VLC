"""TODO:
    This class defines the room dimensions and calculates distnace and angles with respect to two points

    Test this unit and document it

"""

class Room_Geometry:
    def __init__(self, height, width, depth ):
        self.source_point = 0
        self.destination_point = 0

        self.x = depth
        self.y = width
        self.z = height
        self.dimensions = [self.x, self.y, self.z]
        self.Lambertian_order = 1

    def set_points(self, source_point, destination_point ):
        self.source_point = source_point
        self.destination_point = destination_point

    def calculate_distance(self):
        distance = np.linalg.norm(self.source_point - self.destination_point)
        return distance

    def calculate_height_difference(self):
        return abs(self.source_point[2] - self.destination_point[2])

    def calculate_width_difference(self):
        return abs(self.source_point[1] - self.destination_point[1])

    def calculate_depth_difference(self):
        return abs(self.source_point[0] - self.destination_point[0])

    def calc_Tx_angle(self):
        distance = self.calculate_distance(self.source_point, self.destination_point)
        height = self.calculate_height_difference(self.source_point, self.destination_point)
        Tx_angle = np.arccos(height/distance)
        return Tx_angle

    def calc_power_at_dest(self, Tx_power, Element_area):
        Tx_angle = self.calc_Tx_angle()
        Rx_angle = self.Rx_angle()
        distance = calculate_distance()
        Power_at_surface = Tx_power * (self.Lambertian_order+1)* Element_area * \
        np.power(np.cos(Tx_angle),Lambertian_order) * np.cos(Rx_angle) / (4*np.pi * np.power(distance,2))

    def __str__(self):
        return "Room dimensions are {} and the two points are src: {} \
         and dest:{}" .format(self.dimensions, self.source_point, self.destination_point)

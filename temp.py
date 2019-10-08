
import numpy as np

Lambertian_order = 1
Element_area = 5*5  # 5cm*5cm
reflection_coeff = 0.9
Receiver_area = 1  #1 cm^2
Transmit_power = 1e-3  # 1mW
Concentrator_gain = 1
Filter_gain = 1
Speed_of_light = 3e8
Field_Of_View = 90  #degrees

Rx = np.array([2,3,1])
Tx = np.array([3,4,4])
system_settings = {
    "room_dimenstions" : [8,4,3],
    "reflectivity" : 0.7
}

transmitter = {
"location" : np.array([2,3,1]),
"tx_power" :1
}
receiver = {
    "location":np.array([3,4,4]),
    "collection_area":4e-4,
    "Concentrator_gain":1,
    "filter_gain":1,
    "type":0, #FOV
    "FOV": 90
}

Wall = {
    "width":800,
    "height":300,
    "reflectivity":1,
    "area_element":5*5e-4
}

def main():
    #Start with the case of the Rx and Tx on the same floor
    Access_point = np.array([200,400,100])
    Receiver = np.array([200,400,100])

    z=300
    count=0



    Elevation = 0
    #Ceiling
    for x in range(0,400,40): #Step is in centimeters
        for y in range(0,800,40):
            count=count+1
            #Get distance form source
            current_element = np.array([x,y,z])
            distance = np.linalg.norm(Access_point - current_element)
            Total_distance = distance #to keep track of the total distance
            height = Access_point[2] - current_element[2]
            Tx_angle = np.arccos(height/distance)
            Rx_angle = Tx_angle

            delta_x = current_element[0] - Receiver[0]
            delta_y = current_element[1] - Receiver[1]
            Angle = np.arctan(delta_x/delta_y) * 180/np.pi
            #With the dimensions we have , Dx is in the y axis, Dy is in the x_axis
            if(delta_x> = 0):
                if(delta_y> = 0):
                    pass  #First quadrant
                else: #Dy <0 second quadrant
                    Angle = Angle + 180
            else:
                if(delta_y>=0): # dx <0, dy >0 Fourth Quadrant
                    Angle = Angle + 360
                else: #Dy <0 Third quadrant quadrant
                    Angle = Angle + 180

            if(Access_point[1] >= current_element[1]): #The y coordiantes
                    #half of elevation angles (right to the receiver) are
                Elevation = 90 - Field_Of_View * np.pi / 180
            else:
                #and the other half (left from the receiver) is:
                Elevation =  90 + Field_Of_View * np.pi / 180
            print("Element: ", current_element, " Azimuth, ", Angle, " Elevation: ", Elevation)


if __name__=="__main__":
     main()




""" Testing the ADR angles (Az, El) work well is as follows:
# Testing the Elevation angle
1- check that at a given x position, varying the y, The elevation is the same
2- stop the y, changing the x to another one, the elevation increases the further you go away from the prev x
    and gets less the closer you are to the receiver.

#Testing the Az angle:
1- take one side, say right side wall:
    - since Rx is in center, Az will be from -theta to + theta (330 to 30 degres for example )
    - change x. of course y is fixed as the element is in the wall is fixed check the Az changed well

2- Repeat to other walls , moving to back wall, left wall, and then front wall.
They should be increaing to fill 360
"""

"""
TODO:

1- the receiver receives all the rays at a given iteration and then produces some calculations.
so feed the powers and their locations and angles and the receiver calculates the power based on their type

2- Wide FOV evaluation at differnent angles. the narrows the less power, the widest is the upper limit.

3- ADR implementation, which should be between the FOV lower and upper limit. plot that

4- The imaging reciver

Separate the global main caller and the receiver functionality

 {Main}-----> {Reciver}

 The main caller gives the receiver the different rays, their received power value and their angles and times. The receiver calcs based on its type
 - needs a lot of memory, but well organized.
"""



"""This version contains the single matrix of the program,
3D matrix where the first dimension is the wall number
the matrix will have the size of the biggest 2 dimensions of the system. the smaller dimension have zero values
and also use the indexes to index the smaller walls so not all elements are visited.
easier in the code, at the expense of more memory. same computation as only the relevant portions are considered
"""

import numpy as np
import matplotlib.pyplot as plt

Lambertian_order = 1
Element_area = 5*5  # 5cm*5cm
reflection_coeff = 0.9
Receiver_area = 1  #1 cm^2
Transmit_power = 1e-3  # 1mW
Concentrator_gain = 1
Filter_gain = 1
Speed_of_light = 3e8
Field_Of_View = 90  #degrees
Total_received_power = 0

def get_wall_prameters(wall_number):
    dimensions = [400,800] #default for the ceiling
    if(wall_number==1 or wall_number==2): #left and right
        dimensions = [400,300]
    if(wall_number==3 or wall_number==4): #front and back
        dimensions = [300,800]
    return dimensions

def calc_impulse_response(Timing_contribution, Power_contribution, time_bins):
    impulse_response = np.zeros(time_bins)
    for ind  in range(time_bins):
        print(ind)
        grouped_times = np.argwhere((Timing_contribution*1e9<ind+1) & (Timing_contribution*1e9>ind))
        for i in grouped_times:
            impulse_response[ ind]=  impulse_response[ind] + Power_contribution[tuple(i)]
    return impulse_response


def get_Azimuth(current_element, Receiver):
    delta_x = current_element[0] - Receiver[0]
    delta_y = current_element[1] - Receiver[1]

    if(delta_y!=0):
        Angle = np.arctan(delta_x/delta_y) * 180/np.pi
    else:
        Angle=0
    print(current_element, Receiver, delta_x, delta_y, Angle)
    #With the dimensions we have , Dx is in the y axis, Dy is in the x_axis
    if(delta_x>=0):
        if(delta_y>=0):
            pass  #First quadrant
        else: #Dy <0 second quadrant
            Angle = Angle + 180
    else:
        if(delta_y>=0): # dx <0, dy >0 Fourth Quadrant
            Angle = Angle + 360
        else: #Dy <0 Third quadrant quadrant
            Angle = Angle + 180
        pass
    return Angle

def get_elevation(current_element, Receiver):

    #Elevation is defined here as the arc sin of the Z over the R

    R = np.linalg.norm(np.array(current_element) - np.array(Receiver))
    delta_y = current_element[1] - Receiver[1]
    delta_z = current_element[2] - Receiver[2]

    if(R != 0):
        elevation = np.arcsin(delta_z/R) * 180/np.pi
    else:
        elevation=0

    if(delta_y>=0):#element y > rx y (element to the right of receiver )
        pass
    else:
        elevation = elevation+ 90
    print(current_element, Receiver, R, delta_z, elevation)
    return elevation

def calc_power_at_position(receiver_x,receiver_y):
    print("function called at X:", receiver_x, "and Y: ",receiver_y)
    Access_point = np.array([200,400,100])
    Receiver = np.array([receiver_x,receiver_y,100])
    #5 surfaces. biggest is 400*800
    Power_contribution = np.zeros((5,400,800))
    Timing_contribution = np.zeros((5,400,800))

    z=300

    for wall in range(5):
        params = get_wall_prameters(wall)
        [first_dimension, second_dimension] = get_wall_prameters(wall)
        print(params)
        element_width = 5
        for x in range(0,first_dimension,element_width): #Step is in centimeters
            for y in range(0,second_dimension,element_width):
                #Get distance form source
                current_element = np.array([x,y,z])
                distance = np.linalg.norm(Access_point - current_element)
                Total_distance = distance #to keep track of the total distance
                height = abs(Access_point[2] - current_element[2])

                Tx_angle = np.arccos(height/distance)
                #print("height: ", height, "distance 1:", distance, "Tx_angle1: ", Tx_angle*180/np.pi)

                Rx_angle = Tx_angle
                Power_at_surface = Transmit_power * (Lambertian_order+1)* Element_area * \
                np.power(np.cos(Tx_angle),Lambertian_order) * np.cos(Rx_angle) / (4*np.pi * np.power(distance,2))

                #we need to calculate distance from wall element to the receiver
                distance = np.linalg.norm(current_element - Receiver)
                Total_distance = Total_distance+ distance
                Tx_angle = np.arccos(height/distance)
                Rx_angle = Tx_angle
                #print("distance 2:", distance, "Tx_angle2: ", Tx_angle*180/np.pi)

                Power_at_receiver = Power_at_surface * reflection_coeff * Receiver_area * (Lambertian_order+1)* Element_area * \
                    np.power(np.cos(Tx_angle),Lambertian_order) * np.cos(Rx_angle) / (4 * np.pi * np.power(distance,2))
                Power_at_receiver = Power_at_receiver * Concentrator_gain * Filter_gain
                Power_at_receiver = Power_at_receiver*1e10
                #print("power at the receiver: ", Power_at_receiver*1e10)
                #print("reception angle ", Rx_angle)
                #Need to check for the angle is within the field of view.

                if (Rx_angle < Field_Of_View * np.pi/180):
                    #print("angle within FOV", Rx_angle*180/np.pi)
                    Power_contribution[wall,x,y] = Power_at_receiver
                else:
                    #print("angle OUT FOV",Rx_angle*180/np.pi)
                    Power_contribution[wall,x,y] = 0


                    # TODO:
                    #add the azimuth and Elevation classification

                Timing_contribution[wall,x,y] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first

    Total_received_power = np.sum(Power_contribution)
    return [Total_received_power, Timing_contribution, Power_contribution]


def main():
    #Start with the case of the Rx and Tx on the same floor
    """
    [power,time_mat,power_mat]=calc_power_at_position(100,100)
    print("Total_receied power is: ", power)
    imp_resp = calc_impulse_response(time_mat, power_mat, 80)
    plt.plot(imp_resp)
    plt.show()"""
    count_x=0
    count_y=0
    power_xy=np.zeros((10,20))
    x=400
    for y in range(0,900,100):
        print(get_elevation([x,y,300],[200,40,100]))
        #[Angle, Elevation]=get_elevation_and_azimuth_angles([x,y,300],[200,0,100])
    """
    for i in range(10):
        for j in range(20):
            print(i,j)
            [power, Time, powers]=calc_power_at_position(i*40,j*40)
            power_xy[i,j]=power


    plt.imshow(power_xy, cmap='hot', interpolation='nearest')
    plt.show()
    """
    print("main function called ")




if __name__=="__main__":
     main()


"""Testing the program works:
1- check first that get wall parameters returns the correct values
2- check that inside the loop of the get values, each wall has the same range of dimensions

3- test for small values of x,y,z maybe only x,y,z, 400,800,300 with step of 100
4- check the Tx_angle changes according to the calculated values
5- check the azimuth and elevation are correct
"""

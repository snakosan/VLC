
"""
TODO:

1- add FOV check at the receiver,  power received and its timing  ------DONE
2- plot the received power meshgrid
3- calculate the delay spread
4- clean the code, use the dictionary, put all walls into functions
    - possibly make a single function that receives the wall type and substitute the dimensions. since code is the same

5- analyze the power at different elements. is it alot from walls, or from ceiling. if we change the reflection coefficient how does it change it

6- make it more accurate by not considering dimensions below the communication floor
"""



"""
channel calculation:
    1- distnace between the transmitter and receiver
    2- angle of trasmitter and receiver
    3- other equation parameters
    4- reflection types considered
"""
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
    Access_point = np.array([2,4,1])
    Receiver = np.array([1,1,1])
    Power_contribution_ceiling = np.zeros((400,800))
    Power_contribution_left = np.zeros((400,300))
    Power_contribution_right = np.zeros((400,300))
    Power_contribution_front = np.zeros((800,300))
    Power_contribution_back = np.zeros((800,300))


    Timing_contribution_ceiling = np.zeros((400,800))
    Timing_contribution_left = np.zeros((400,300))
    Timing_contribution_right = np.zeros((400,300))
    Timing_contribution_front = np.zeros((800,300))
    Timing_contribution_back = np.zeros((800,300))
    z=3
    count=0


    #Ceiling: x,y
    #left right : x,z
    #front-back : y,z

    first_dimension = 400
    second_dimension = 800
    element_width = 5
    wall="Cieling"
    for x in range(0,first_dimension,element_width): #Step is in centimeters
        for y in range(0,second_dimension,element_width):
            count=count+1
            #Get distance form source
            current_element =np.array([x,y,z])
            distance = np.linalg.norm(Access_point - current_element)
            Total_distance = distance #to keep track of the total distance
            height = Access_point[2] - current_element[2]
            Tx_angle = np.arccos(height/distance)
            Rx_angle=Tx_angle
            Power_at_surface = Transmit_power * (Lambertian_order+1)* Element_area * \
            np.power(np.cos(Tx_angle),Lambertian_order) * np.cos(Rx_angle) / (4*np.pi * np.power(distance,2))

            #we need to calculate distance from wall element to the receiver
            distance = np.linalg.norm(current_element - Receiver)
            Total_distance = Total_distance+ distance
            Tx_angle = np.arccos(height/distance)
            Rx_angle=Tx_angle
            Power_at_receiver = Power_at_surface * reflection_coeff * Receiver_area * (Lambertian_order+1)* Element_area * \
                np.power(np.cos(Tx_angle),Lambertian_order) * np.cos(Rx_angle) / (4 * np.pi * np.power(distance,2))
            Power_at_receiver = Power_at_receiver * Concentrator_gain * Filter_gain
            #contribution form the ceiling  at x and y position

            #Need to check for the angle is within the field of view.
            #print(Rx_angle*180/np.pi)
            if (Rx_angle < Field_Of_View * np.pi/180):
                if(wall=="Cieling"):
                    Power_contribution_ceiling[x,y] = Power_at_receiver
                    Timing_contribution_ceiling[x,y] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first
                if(wall=="left"):
                    Power_contribution_right[x,z] = Power_at_receiver
                    Timing_contribution_ceiling[x,z] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first
                if(wall=="left"):
                    Power_contribution_left[x,z] = Power_at_receiver
                    Timing_contribution_ceiling[x,z] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first
                if(wall=="front"):
                    Power_contribution_front[y,z] = Power_at_receiver
                    Timing_contribution_ceiling[y,z] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first
                if(wall=="back"):
                    Power_contribution_back[y,z] = Power_at_receiver
                    Timing_contribution_ceiling[y,z] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first

            else:
                if(wall=="Cieling"):
                    Power_contribution_ceiling[x,y] = 0
                    Timing_contribution_ceiling[x,y] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first
                if(wall=="left"):
                    Power_contribution_right[x,z] = 0
                    Timing_contribution_ceiling[x,z] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first
                if(wall=="left"):
                    Power_contribution_left[x,z] = 0
                    Timing_contribution_ceiling[x,z] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first
                if(wall=="front"):
                    Power_contribution_front[y,z] = 0
                    Timing_contribution_ceiling[y,z] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first
                if(wall=="back"):
                    Power_contribution_back[y,z] = 0
                    Timing_contribution_ceiling[y,z] = (Total_distance* 10**-2)  / Speed_of_light  #converting from cm to m first

        
    print("main function called ")


if __name__=="__main__":
     main()

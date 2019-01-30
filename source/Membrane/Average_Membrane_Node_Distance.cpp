#include "Membrane.h"
#include "General_functions.hpp"

double Membrane::Average_Node_Distance()
{
    double average_membrane_Node_distance=0.0;
    double temp[3];
    double length;
    int Node_A, Node_B;
    for (int i=0; i<Num_of_Triangle_Pairs; i++)
    {
        Node_A=Node_Bond_list[i][0];
        Node_B=Node_Bond_list[i][1];
        temp[0]=Node_Position[Node_A][0]-Node_Position[Node_B][0];
        temp[1]=Node_Position[Node_A][1]-Node_Position[Node_B][1];
        temp[2]=Node_Position[Node_A][2]-Node_Position[Node_B][2];
        length= vector_length(temp);
        average_membrane_Node_distance+=length;
    }
    average_membrane_Node_distance=average_membrane_Node_distance/Num_of_Triangle_Pairs;
    return(average_membrane_Node_distance);
    
}
 

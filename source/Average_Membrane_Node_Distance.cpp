#include "Membrane.h"
#include "General_functions.hpp"

double Membrane::Average_Membrane_Node_Distance()
{
	double average_membrane_Node_distance=0.0;
	double temp[3];
	double temp_length;
	int temp_Node_A,temp_Node_B;
	for (int i=0; i<Num_of_Triangle_Pairs; i++)
	{
		temp_Node_A=Membrane_Edges[i][0];
		temp_Node_B=Membrane_Edges[i][1];
		temp[0]=Membrane_Node_Position[temp_Node_A][0]-Membrane_Node_Position[temp_Node_B][0];
        temp[1]=Membrane_Node_Position[temp_Node_A][1]-Membrane_Node_Position[temp_Node_B][1];
        temp[2]=Membrane_Node_Position[temp_Node_A][2]-Membrane_Node_Position[temp_Node_B][2];
		temp_length= vectorlength(temp);
		average_membrane_Node_distance+=temp_length;
	}
	average_membrane_Node_distance=average_membrane_Node_distance/Num_of_Triangle_Pairs;
	return(average_membrane_Node_distance);
	
}
 

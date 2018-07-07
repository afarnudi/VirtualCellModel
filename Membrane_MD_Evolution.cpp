#include "Membrane.h"


 void Membrane::Membrane_MD_Evolution()
 {
	         for(int j=0 ; j<Membrane_num_of_Nodes ; j++)
        {
            Membrane_Node_Position[j][0] += Membrane_Node_Velocity[j][0]*MD_Time_Step - Membrane_Node_Force[j][0]*MD_Time_Step*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Position[j][1] += Membrane_Node_Velocity[j][1]*MD_Time_Step - Membrane_Node_Force[j][1]*MD_Time_Step*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Position[j][2] += Membrane_Node_Velocity[j][2]*MD_Time_Step - Membrane_Node_Force[j][2]*MD_Time_Step*MD_Time_Step/(Membrane_Node_Mass*2.0);
        }
//
        
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)  // loop to count every particle  and update its velocity
        {
            Membrane_Node_Velocity[j][0] += - Membrane_Node_Force[j][0]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][1] += - Membrane_Node_Force[j][1]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][2] += - Membrane_Node_Force[j][2]*MD_Time_Step/(Membrane_Node_Mass*2.0);
        }
//
        
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)
        {
            Membrane_Node_Force[j][0]=0.0;
            Membrane_Node_Force[j][1]=0.0;
            Membrane_Node_Force[j][2]=0.0;
        }
 }
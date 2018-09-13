#include "Membrane.h"


 void Membrane::MD_Evolution(){

        for(int j=0 ; j<Num_of_Nodes ; j++){
            Node_Position[j][0] += Node_Velocity[j][0]*MD_Time_Step - Node_Force[j][0]*MD_Time_Step*MD_Time_Step/(Node_Mass*2.0);
            Node_Position[j][1] += Node_Velocity[j][1]*MD_Time_Step - Node_Force[j][1]*MD_Time_Step*MD_Time_Step/(Node_Mass*2.0);
            Node_Position[j][2] += Node_Velocity[j][2]*MD_Time_Step - Node_Force[j][2]*MD_Time_Step*MD_Time_Step/(Node_Mass*2.0);


            Node_Velocity[j][0] += - Node_Force[j][0]*MD_Time_Step/(Node_Mass*2.0);
            Node_Velocity[j][1] += - Node_Force[j][1]*MD_Time_Step/(Node_Mass*2.0);
            Node_Velocity[j][2] += - Node_Force[j][2]*MD_Time_Step/(Node_Mass*2.0);


            Node_Force[j][0]=0.0;
            Node_Force[j][1]=0.0;
            Node_Force[j][2]=0.0;
        }

 }

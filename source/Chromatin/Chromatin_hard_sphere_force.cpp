//
//  Membrane_spring_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Chromatin.h"


void Chromatin::hard_sphere (void){
    double le1,lmin;
    double deltax,deltay,deltaz,temp_Node_distance,temp_force;
    double temp_potential_energy = 0.0;
    int temp_Node_A, temp_Node_B;
   
    le1=2.3*Node_radius;
    lmin=2.*Node_radius;
    
    for (int k=0 ; k< Num_of_Nodes-2 ; k++)
    {
        temp_Node_B=k;
        for (int j=k+2; j<Num_of_Nodes; j++) {
            
            temp_Node_A=j;
            
            deltax=Node_Position[temp_Node_A][0]-Node_Position[temp_Node_B][0];
            deltay=Node_Position[temp_Node_A][1]-Node_Position[temp_Node_B][1];
            deltaz=Node_Position[temp_Node_A][2]-Node_Position[temp_Node_B][2];
            
            
            temp_Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
            temp_force=0.0;
            double temp_exp_le1=exp(1.0/(temp_Node_distance-le1));
            
            if(temp_Node_distance < le1   &&  temp_Node_distance > lmin  )  // repulsive force
            {
                temp_force= -(Spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin))*( 1.0/(temp_Node_distance-lmin) + 1.0/((temp_Node_distance-le1)*(temp_Node_distance-le1)));                 // force on i th from j
                temp_potential_energy= Spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin);
            }
            if(temp_force<-1000.05   ||  temp_Node_distance<lmin )
            {
                
                temp_force =-1000.05-Spring_force_cutt_off* ( 0.671965*Node_radius - temp_Node_distance );
                temp_potential_energy = 1.85038 + 1005.05 * ( 0.671965*Node_radius - temp_Node_distance )+0.5*Spring_force_cutt_off*( 0.671965*Node_radius - temp_Node_distance )*( 0.671965*Node_radius - temp_Node_distance );
                //            cout<<"temp_force_after_cut_off"<<temp_force<<endl;
            }
            
            Total_Potential_Energy += temp_potential_energy;
            
            // implimentation of forces:
            Node_Force[temp_Node_A][0] += -temp_force*deltax/temp_Node_distance-Damping_coefficient*(Node_Velocity[temp_Node_A][0]-Node_Velocity[temp_Node_B][0]);
            Node_Force[temp_Node_A][1] += -temp_force*deltay/temp_Node_distance-Damping_coefficient*(Node_Velocity[temp_Node_A][1]-Node_Velocity[temp_Node_B][1]);
            Node_Force[temp_Node_A][2] += -temp_force*deltaz/temp_Node_distance-Damping_coefficient*(Node_Velocity[temp_Node_A][2]-Node_Velocity[temp_Node_B][2]);
            
            Node_Force[temp_Node_B][0] += temp_force*deltax/temp_Node_distance+Damping_coefficient*(Node_Velocity[temp_Node_A][0]-Node_Velocity[temp_Node_B][0]); //from j  to i
            Node_Force[temp_Node_B][1] += temp_force*deltay/temp_Node_distance+Damping_coefficient*(Node_Velocity[temp_Node_A][1]-Node_Velocity[temp_Node_B][1]);
            Node_Force[temp_Node_B][2] += temp_force*deltaz/temp_Node_distance+Damping_coefficient*(Node_Velocity[temp_Node_A][2]-Node_Velocity[temp_Node_B][2]);
        }
    }
    // End of Membrane Node Pair forces
    
}

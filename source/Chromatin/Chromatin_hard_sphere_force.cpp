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
    double deltax,deltay,deltaz,Node_distance,temp_force;
    double temp_potential_energy = 0.0;
    int Node_A, Node_B;
   
    le1=2.3*Node_radius;
    lmin=2.*Node_radius;
    
    for (int k=0 ; k< Num_of_Nodes-2 ; k++)
    {
        Node_B=k;
        for (int j=k+2; j<Num_of_Nodes; j++) {
            
            Node_A=j;
            
            deltax=Node_Position[Node_A][0]-Node_Position[Node_B][0];
            deltay=Node_Position[Node_A][1]-Node_Position[Node_B][1];
            deltaz=Node_Position[Node_A][2]-Node_Position[Node_B][2];
            
            Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
            temp_force=0.0;
            
            if(Node_distance < le1   &&  Node_distance > lmin  )  // repulsive force
            {
                double exp_le1=exp(1.0/(Node_distance-le1));
                temp_force = ( (Spring_coefficient*exp_le1)/(Node_distance-lmin) )*( 1/(Node_distance-lmin)+1/( (Node_distance-le1)*(Node_distance-le1) ) );
                temp_potential_energy = Spring_coefficient*exp_le1/(Node_distance-lmin);
            }
            if(temp_force<-1000   ||  Node_distance<lmin )
            {
                double c=1000;
                temp_force = -2.0*c*Node_distance/(Node_radius) + 3*c;
                temp_potential_energy=   c*Node_distance*Node_distance/(Node_radius) - 3*c*Node_distance;
            }
            
            Total_Potential_Energy += temp_potential_energy;
            
            //damper
            if (Damping_coefficient>0.00001) {
                double delta_v[3] = {Node_Velocity[Node_A][0]-Node_Velocity[Node_B][0],
                    Node_Velocity[Node_A][1]-Node_Velocity[Node_B][1],
                    Node_Velocity[Node_A][2]-Node_Velocity[Node_B][2]};
                double delta_r[3]={deltax, deltay, deltaz};
                double temp_damp=innerproduct(delta_r, delta_v);
                double delta_r_2=vector_length_squared(delta_r);
                temp_damp=temp_damp/delta_r_2;
                
                Node_Force[Node_A][0] +=  - Damping_coefficient*temp_damp*deltax;
                Node_Force[Node_A][1] +=  - Damping_coefficient*temp_damp*deltay;
                Node_Force[Node_A][2] +=  - Damping_coefficient*temp_damp*deltaz;
                
                Node_Force[Node_B][0] += Damping_coefficient*temp_damp*deltax;
                Node_Force[Node_B][1] += Damping_coefficient*temp_damp*deltay;
                Node_Force[Node_B][2] += Damping_coefficient*temp_damp*deltaz;
            }
            
            temp_force=temp_force/Node_distance;
            // implimentation of forces:
            Node_Force[Node_A][0] +=  temp_force*deltax;
            Node_Force[Node_A][1] +=  temp_force*deltay;
            Node_Force[Node_A][2] +=  temp_force*deltaz;
            
            Node_Force[Node_B][0] += -temp_force*deltax;
            Node_Force[Node_B][1] += -temp_force*deltay;
            Node_Force[Node_B][2] += -temp_force*deltaz;
        }
    }
    
}

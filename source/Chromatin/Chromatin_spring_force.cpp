//
//  Membrane_spring_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Chromatin.h"


void Chromatin::potential_1 (void){
    double le0,le1,lmax,lmin;
    double deltax,deltay,deltaz,Node_distance,temp_force;
    double temp_potential_energy = 0.0;
    int Node_A, Node_B;
    
    lmin=2.*Node_radius;
    lmax=2.4*Node_radius;
    
    le1=2.1*Node_radius;
    le0=2.3*Node_radius;
    
    
    Total_Potential_Energy=0.0;
    
    for (int k=0 ; k< Num_of_Nodes-1 ; k++)
    {
        Node_B=k;
        Node_A=k+1;
        
        deltax=Node_Position[Node_A][0]-Node_Position[Node_B][0];
        deltay=Node_Position[Node_A][1]-Node_Position[Node_B][1];
        deltaz=Node_Position[Node_A][2]-Node_Position[Node_B][2];
        
        Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;
        
        if(Node_distance > le1  && Node_distance <lmax )  //bondforce
        {
            double exp_le0=exp(1.0/(le0-Node_distance));
            temp_force = ( (Spring_coefficient*exp_le0)/(Node_distance-lmax) )*( 1/(lmax-Node_distance)+1/( (le0-Node_distance)*(le0-Node_distance) ) );
            temp_potential_energy = Spring_coefficient*exp_le0/(lmax-Node_distance);
            
        } else if(Node_distance < le0   &&  Node_distance > lmin  )  // repulsive force
        {
            double exp_le1=exp(1.0/(Node_distance-le1));
            temp_force = ( (Spring_coefficient*exp_le1)/(Node_distance-lmin) )*( 1/(Node_distance-lmin)+1/( (Node_distance-le1)*(Node_distance-le1) ) );
            temp_potential_energy = Spring_coefficient*exp_le1/(Node_distance-lmin);
        }
        if(Node_distance>lmax || Node_distance<lmin)
        {
            cout<<"k == "<<k<<"\tNode distance = "<<Node_distance<<"\tF = "<<temp_force<<endl;
            cout<<"The potential between the Chromatin nodes is too weak for the current temperture of the system. Or the node potential cannot sustain the applied stress. As a result, a node pair distance has exceed the allowed regien defined by the node pairwise potential. Please adjust the configuration of the springs and restart the run.\n";
            exit(EXIT_FAILURE);
        }
        
        Total_Potential_Energy += temp_potential_energy;
        
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
        
        // implimentation of forces:
        Node_Force[Node_A][0] +=  temp_force*deltax/Node_distance;
        Node_Force[Node_A][1] +=  temp_force*deltay/Node_distance;
        Node_Force[Node_A][2] +=  temp_force*deltaz/Node_distance;
        
        Node_Force[Node_B][0] += -temp_force*deltax/Node_distance;
        Node_Force[Node_B][1] += -temp_force*deltay/Node_distance;
        Node_Force[Node_B][2] += -temp_force*deltaz/Node_distance;
    }
    // End of Membrane Node Pair forces
    
}


void Chromatin::FENE(void){
    
    double equi_point, delta_r_max, epsilon;
    double deltax,deltay,deltaz,temp_Node_distance,temp_force;
    //    int pos1,pos2,pos3,pos4;  // to making calculation of surface force easier
    double temp_potential_energy = 0.0;
    //    double temp_p1[3], temp_p2[3], temp_p3[3], temp_p4[3];
    //    double  N1[3], N2[3], N3[3], p3p1[3], p3p2[3], p4p2[3], p4p1[3], Ep2p1[3], sinus, F0, F1[3], F2[3], F3[3], F4[3];// for exmple p3p1 is p3-p1 and ....
    
    /// calculate network force:
    int temp_Node_A, temp_Node_B;
    // le0= 0.994;
    //lmax=1.170;
    //le1=0.654;
    //lmin=0.478;
    
    equi_point=1.05*Node_radius;
    delta_r_max=Node_radius;
    epsilon=Spring_coefficient;
    //Total_Potential_Energy=0.0;
    
    for (int k=1 ; k< Num_of_Nodes-1 ; k++)
    {
        temp_Node_B=k-1;
        temp_Node_A=k;
        
        deltax=Node_Position[temp_Node_A][0]-Node_Position[temp_Node_B][0];
        deltay=Node_Position[temp_Node_A][1]-Node_Position[temp_Node_B][1];
        deltaz=Node_Position[temp_Node_A][2]-Node_Position[temp_Node_B][2];
        
        
        temp_Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;
        
        if (temp_Node_distance<(equi_point-delta_r_max) || temp_Node_distance>(equi_point+delta_r_max)){
            cout<<"Node distance out of bounds of the FENE cut off.\nNode numbers "<<temp_Node_A<<" and "<<temp_Node_B<<endl;
            cout<<"Node distance "<<temp_Node_distance<<endl;
            exit(EXIT_FAILURE);
            
        } else {
            double temp_R=(temp_Node_distance-equi_point)/delta_r_max;
            temp_potential_energy=-(0.5)*epsilon*delta_r_max*delta_r_max*log(1-temp_R*temp_R);
            temp_force=epsilon*(temp_Node_distance-equi_point)/(1-temp_R*temp_R);
            
            
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

void Chromatin::Strong_spring (void){
    
    double deltax,deltay,deltaz,temp_Node_distance,temp_force;
    
    int temp_Node_A, temp_Node_B;
    
    for (int k=0 ; k< Num_of_Nodes-1 ; k++)
    {
        temp_Node_B=k;
        temp_Node_A=k+1;
        
        deltax=Node_Position[temp_Node_A][0]-Node_Position[temp_Node_B][0];
        deltay=Node_Position[temp_Node_A][1]-Node_Position[temp_Node_B][1];
        deltaz=Node_Position[temp_Node_A][2]-Node_Position[temp_Node_B][2];
        
        
        temp_Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
       
        temp_force=Spring_coefficient*(temp_Node_distance-2.1*Node_radius);
        
        
        // implimentation of forces:
        Node_Force[temp_Node_A][0] += -temp_force*deltax/temp_Node_distance+Damping_coefficient*(Node_Velocity[temp_Node_A][0]-Node_Velocity[temp_Node_B][0]);
        Node_Force[temp_Node_A][1] += -temp_force*deltay/temp_Node_distance+Damping_coefficient*(Node_Velocity[temp_Node_A][1]-Node_Velocity[temp_Node_B][1]);
        Node_Force[temp_Node_A][2] += -temp_force*deltaz/temp_Node_distance+Damping_coefficient*(Node_Velocity[temp_Node_A][2]-Node_Velocity[temp_Node_B][2]);
        
        Node_Force[temp_Node_B][0] += temp_force*deltax/temp_Node_distance-Damping_coefficient*(Node_Velocity[temp_Node_A][0]-Node_Velocity[temp_Node_B][0]); //from j  to i
        Node_Force[temp_Node_B][1] += temp_force*deltay/temp_Node_distance-Damping_coefficient*(Node_Velocity[temp_Node_A][1]-Node_Velocity[temp_Node_B][1]);
        Node_Force[temp_Node_B][2] += temp_force*deltaz/temp_Node_distance-Damping_coefficient*(Node_Velocity[temp_Node_A][2]-Node_Velocity[temp_Node_B][2]);
    }
    // End of Membrane Node Pair forces
}

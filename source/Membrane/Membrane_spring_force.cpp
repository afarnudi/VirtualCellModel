//
//  Membrane_spring_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"


void Membrane::potential_1 (void){
    double le0,le1,lmax,lmin;
    double deltax, deltay, deltaz, Node_distance, temp_force;
    double temp_potential_energy = 0.0;
    
    int Node_A, Node_B;
    
    lmax=Max_node_pair_length;
    lmin=Min_node_pair_length;
    
    le0=lmin+3*(lmax-lmin)/4.0;
    le1=lmin+(lmax-lmin)/4.0;
    
    Total_Potential_Energy=0.0;
    
    for (int k=0 ; k< Num_of_Node_Pairs ; k++)
    {
        Node_B=Node_Bond_list[k][0];
        Node_A=Node_Bond_list[k][1];
        
        deltax=Node_Position[Node_A][0]-Node_Position[Node_B][0];
        deltay=Node_Position[Node_A][1]-Node_Position[Node_B][1];
        deltaz=Node_Position[Node_A][2]-Node_Position[Node_B][2];
        
        Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;
        
        if(Node_distance >le1  & Node_distance < le0 )  //free zone
        {
            temp_potential_energy=0 ; // free zone
        } else if(Node_distance > le0  & Node_distance <lmax )  //bondforce
        {
            double exp_le0=exp(1.0/(le0-Node_distance));
            temp_force = ( (Spring_coefficient*exp_le0)/(Node_distance-lmax) )*( 1/(lmax-Node_distance)+1/( (le0-Node_distance)*(le0-Node_distance) ) );
            temp_potential_energy = Spring_coefficient*exp_le0/(lmax-Node_distance);
            
        } else if(Node_distance < le1   &  Node_distance > lmin  )  // repulsive force
        {
            double exp_le1=exp(1.0/(Node_distance-le1));
            temp_force = ( (Spring_coefficient*exp_le1)/(Node_distance-lmin) )*( 1/(Node_distance-lmin)+1/( (Node_distance-le1)*(Node_distance-le1) ) );
            temp_potential_energy = Spring_coefficient*exp_le1/(Node_distance-lmin);
        }
        
        if(Node_distance>lmax || Node_distance<lmin)
        {
            cout<<"k == "<<k<<"\tNode distance = "<<Node_distance<<"\tF = "<<temp_force<<endl;
            cout<<"The potential between the Membrane nodes is too weak for the current temperture of the system. Or the node potential cannot sustain the applied stress. As a result, a node pair distance has exceed the allowed regien defined by the node pairwise potential. Please adjust the configuration of the springs and restart the run.\n";
            exit(EXIT_FAILURE);
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
        
        // implimentation of forces:
        Node_Force[Node_A][0] +=  temp_force*deltax/Node_distance;
        Node_Force[Node_A][1] +=  temp_force*deltay/Node_distance;
        Node_Force[Node_A][2] +=  temp_force*deltaz/Node_distance;
        
        Node_Force[Node_B][0] += -temp_force*deltax/Node_distance;
        Node_Force[Node_B][1] += -temp_force*deltay/Node_distance;
        Node_Force[Node_B][2] += -temp_force*deltaz/Node_distance;
    }
    
}
void Membrane::potential_2 (void){
    
    //cout<<"Houkian Spring"<<endl;
    
    double deltax,deltay,deltaz,Node_distance,temp_force;
    
    int Node_A, Node_B;
    
    for (int k=0 ; k< Num_of_Node_Pairs ; k++)
    {
        Node_B=Node_Bond_list[k][0];
        Node_A=Node_Bond_list[k][1];
        
        deltax=Node_Position[Node_A][0]-Node_Position[Node_B][0];
        deltay=Node_Position[Node_A][1]-Node_Position[Node_B][1];
        deltaz=Node_Position[Node_A][2]-Node_Position[Node_B][2];
        
        
        Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        
        temp_force=-Spring_coefficient*(Node_distance-Average_node_pair_length);
        
        
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
        
        
        // implimentation of forces:
        Node_Force[Node_A][0] +=  temp_force*deltax/Node_distance;
        Node_Force[Node_A][1] +=  temp_force*deltay/Node_distance;
        Node_Force[Node_A][2] +=  temp_force*deltaz/Node_distance;
        
        Node_Force[Node_B][0] += -temp_force*deltax/Node_distance;
        Node_Force[Node_B][1] += -temp_force*deltay/Node_distance;
        Node_Force[Node_B][2] += -temp_force*deltaz/Node_distance;
    }
    
}

void Membrane::FENE(void){
    
    double equi_point, delta_r_max, epsilon;
    double deltax,deltay,deltaz,Node_distance,temp_force;
    double temp_potential_energy = 0.0;
   
    
    /// calculate network force:
    int Node_A, Node_B;
    
    equi_point=(Max_node_pair_length + Min_node_pair_length)/2;
    delta_r_max=Max_node_pair_length - equi_point;
    delta_r_max+=delta_r_max*0.2;
    epsilon=Spring_coefficient;
    //Total_Potential_Energy=0.0;
    
    for (int k=0 ; k< Num_of_Node_Pairs ; k++)
    {
        Node_B=Node_Bond_list[k][0];
        Node_A=Node_Bond_list[k][1];
        
        deltax=Node_Position[Node_A][0]-Node_Position[Node_B][0];
        deltay=Node_Position[Node_A][1]-Node_Position[Node_B][1];
        deltaz=Node_Position[Node_A][2]-Node_Position[Node_B][2];
        
        
        Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;
        
        if (Node_distance<(equi_point-delta_r_max) || Node_distance>(equi_point+delta_r_max)){
            
            cout<<"Node distance out of bounds of the FENE cut off.\nNode numbers "<<Node_A<<" and "<<Node_B<<endl;
            cout<<"Node distance "<<Node_distance<<endl;
            exit(EXIT_FAILURE);
            
        } else {
            double temp_R=(Node_distance-equi_point)/delta_r_max;
            temp_potential_energy=-(0.5)*epsilon*delta_r_max*delta_r_max*log(1-temp_R*temp_R);
            temp_force=-epsilon*(Node_distance-equi_point)/(1-temp_R*temp_R);
            
            
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
            
            
            // implimentation of forces:
            Node_Force[Node_A][0] +=  temp_force*deltax/Node_distance;
            Node_Force[Node_A][1] +=  temp_force*deltay/Node_distance;
            Node_Force[Node_A][2] +=  temp_force*deltaz/Node_distance;
            
            Node_Force[Node_B][0] += -temp_force*deltax/Node_distance;
            Node_Force[Node_B][1] += -temp_force*deltay/Node_distance;
            Node_Force[Node_B][2] += -temp_force*deltaz/Node_distance;
        }
    }
    
}

void Membrane::Relaxation_potential (void){
    
    double le0,le1,lmax,lmin;
    double deltax, deltay, deltaz, Node_distance, temp_force;
    double temp_potential_energy = 0.0;
    int Node_A, Node_B;
    
    lmax=Max_node_pair_length - 0.08*Max_node_pair_length;
    lmin=Min_node_pair_length + 0.16*Min_node_pair_length;
    
    le0=lmin+3*(lmax-lmin)/4.0;
    le1=lmin+(lmax-lmin)/4.0;
    
    Total_Potential_Energy=0.0;
    
    for (int k=0 ; k< Num_of_Node_Pairs ; k++)
    {
        Node_B=Node_Bond_list[k][0];
        Node_A=Node_Bond_list[k][1];
        
        deltax=Node_Position[Node_A][0]-Node_Position[Node_B][0];
        deltay=Node_Position[Node_A][1]-Node_Position[Node_B][1];
        deltaz=Node_Position[Node_A][2]-Node_Position[Node_B][2];
        
        
        Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;
        
        if(Node_distance >le1  & Node_distance < le0 )  //free zone
        {
            temp_potential_energy=0 ; // free zone
        } else if(Node_distance > le0  & Node_distance <lmax )  //bondforce
        {
            double exp_le0=exp(1.0/(le0-Node_distance));
            temp_force = ( (Spring_coefficient*exp_le0)/(Node_distance-lmax) )*( 1/(lmax-Node_distance)+1/( (le0-Node_distance)*(le0-Node_distance) ) );
            temp_potential_energy = Spring_coefficient*exp_le0/(lmax-Node_distance);
            
        } else if(Node_distance < le1   &  Node_distance > lmin  )  // repulsive force
        {
            double exp_le1=exp(1.0/(Node_distance-le1));
            temp_force = ( (Spring_coefficient*exp_le1)/(Node_distance-lmin) )*( 1/(Node_distance-lmin)+1/( (Node_distance-le1)*(Node_distance-le1) ) );
            temp_potential_energy = Spring_coefficient*exp_le1/(Node_distance-lmin);
        }
        /// my cutoff for force amplitute and for avoiding leting particle scape from force trap
        if(temp_force>1000  || Node_distance>lmax )
        {
            double c=1500;
            temp_force = -2.0*c*Node_distance/(lmax-lmin)+( (lmin+lmax)/(lmax-lmin) )*c;
            temp_potential_energy=   c*Node_distance*Node_distance/(lmax-lmin)+( (lmin+lmax)/(lmin-lmax) )*c*Node_distance;
        }
        
        
        if(temp_force<-1000   ||  Node_distance<lmin )
        {
            double c=1500;
            temp_force = -2.0*c*Node_distance/(lmax-lmin)+( (lmin+lmax)/(lmax-lmin) )*c;
            temp_potential_energy=   c*Node_distance*Node_distance/(lmax-lmin)+( (lmin+lmax)/(lmin-lmax) )*c*Node_distance;
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
        
        
        // implimentation of forces:
        Node_Force[Node_A][0] +=  temp_force*deltax/Node_distance;
        Node_Force[Node_A][1] +=  temp_force*deltay/Node_distance;
        Node_Force[Node_A][2] +=  temp_force*deltaz/Node_distance;
        
        Node_Force[Node_B][0] += -temp_force*deltax/Node_distance;
        Node_Force[Node_B][1] += -temp_force*deltay/Node_distance;
        Node_Force[Node_B][2] += -temp_force*deltaz/Node_distance;
    }
//    cout<<"here"<<endl;
}

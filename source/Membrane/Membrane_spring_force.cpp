//
//  Membrane_spring_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"
#include "General_functions.hpp"

void Membrane::FENE_log (void){
    double le0, le1, lmax, lmin;
    double deltax, deltay, deltaz, Node_distance, temp_force;
    double temp_potential_energy = 0.0;

    int Node_A, Node_B;

    double width= 2*Average_node_pair_length/3; // it defines a minimum limit for weil width. 2/3 is an important factor: it dose not let the bondes to go beyond the allowed zone for bending force.
    double delta_length=Max_node_pair_length-Min_node_pair_length;
    double width_scaling =0.2; // this variable adjusts the log_barrier potential width. if the edges have almost the same length, then it shlould be  redefined to have an appropriate weil width
    if ((1+ 2*width_scaling)*delta_length < width ){ //this condition shows the case when the mesh is almost ordered.
        width_scaling= (width-delta_length)/(2*delta_length); //in this case the width tuned in a way that the weil width becomes exactly 0.66*Average_node_pair_lenght
//        cout<<"\n\nif == true\n\n";
    }
    else{
        width=  (1+ 2*width_scaling)*delta_length; //for disordered meshes it sets the weil witdh by scaling factor 0.2. in this case, the width may become larger than 0.66Average which was the minimum limit.
//        cout<<"\n\nif == false\n\n";
    }
 
    lmin= Min_node_pair_length - width_scaling*delta_length;
    lmax= Max_node_pair_length +  width_scaling*delta_length;

    le0 = lmin + 3*(lmax-lmin)/4;
    le1 = lmin +   (lmax-lmin)/4;
    
//    cout<<"\nlmin = "<<lmin<<"lmax = "<<lmax<<"le0 = "<<le0<<"le1 = "<<le1<<endl<<endl;
    
    Total_Potential_Energy=0.0;

    for (int k=0 ; k< Num_of_Node_Pairs ; k++)
    {
        Node_B=Node_Bond_list[k][0];
        Node_A=Node_Bond_list[k][1];
        //        DamperCheck[k]= 0.0;
        deltax=Node_Position[Node_A][0]-Node_Position[Node_B][0];
        deltay=Node_Position[Node_A][1]-Node_Position[Node_B][1];
        deltaz=Node_Position[Node_A][2]-Node_Position[Node_B][2];

        Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;

        if(Node_distance > le0  & Node_distance <lmax )  //bondforce
        {
            double exp_le0=exp(1.0/(le0-Node_distance));
            temp_force = ( (Spring_coefficient*exp_le0)/(Node_distance-lmax) )*( 1.0/(lmax-Node_distance)+1.0/( (le0-Node_distance)*(le0-Node_distance) ) );
            temp_potential_energy = Spring_coefficient*exp_le0/(lmax-Node_distance);

        } else if(Node_distance < le1   &  Node_distance > lmin  )  // repulsive force
        {
            double exp_le1=exp(1.0/(Node_distance-le1));
            temp_force = ( (Spring_coefficient*exp_le1)/(Node_distance-lmin) )*( 1.0/(Node_distance-lmin)+1.0/( (Node_distance-le1)*(Node_distance-le1) ) );
            temp_potential_energy = Spring_coefficient*exp_le1/(Node_distance-lmin);
        }

        if(Node_distance>lmax || Node_distance<lmin)
        {
            if (generalParameters.Num_of_Actins == 0) {
                cout<<"k == "<<k<<"\tNode distance = "<<Node_distance<<"\tF = "<<temp_force<<endl;
                cout<<"Node A = "<<Node_A<<"\tNode B = "<<Node_B<<endl;
                cout<<"The potential between the Membrane nodes is too weak for the current temperture of the system. Or the node potential cannot sustain the applied stress. As a result, a node pair distance has exceed the allowed regien defined by the node pairwise potential. Please adjust the configuration of the springs and restart the run.\n";

                exit(EXIT_FAILURE);
            } else {
                if(temp_force>1000  || Node_distance>lmax )
                {
                    double c=800;
                    temp_force = -2.0*c*Node_distance/(lmax-lmin)+( (lmin+lmax)/(lmax-lmin) )*c;
                    temp_potential_energy=   c*Node_distance*Node_distance/(lmax-lmin)+( (lmin+lmax)/(lmin-lmax) )*c*Node_distance;
                }


                if(temp_force<-1000   ||  Node_distance<lmin )
                {
                    double c=800;
                    temp_force = -2.0*c*Node_distance/(lmax-lmin)+( (lmin+lmax)/(lmax-lmin) )*c;
                    temp_potential_energy=   c*Node_distance*Node_distance/(lmax-lmin)+( (lmin+lmax)/(lmin-lmax) )*c*Node_distance;
                }
            }
        }

        Total_Potential_Energy += temp_potential_energy;

        //damper  we have a Bug here! when we turn on the Damper, the memberane will rotate after some MD steps. so it is better not to use any Damper Coeficient other than zero.
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
void Membrane::custom(void){

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

            //damper  we have a Bug here! when we turn on the Damper, the memberane will rotate after some MD steps. so it is better not to use any Damper Coeficient other than zero.
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

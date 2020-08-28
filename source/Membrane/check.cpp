//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"


void Membrane::check(void){
    Min_node_pair_length=1000;
    Max_node_pair_length=0;
    Average_node_pair_length=0;
    
    for (int i=0; i<(Num_of_Node_Pairs- Num_of_Free_Bonds); i++) {
        double dist=0;
        dist=sqrt((Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])*(Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])+(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])*(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])+(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2])*(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2]));
        Average_node_pair_length+=dist;
        
        if (Min_node_pair_length>dist) {
            Min_node_pair_length=dist;
        }
        if (dist>Max_node_pair_length) {
            Max_node_pair_length=dist;
        }
    }
    Average_node_pair_length/=(Num_of_Node_Pairs-Num_of_Free_Bonds);
    if (!GenConst::Testmode) {
        cout<<"Node pair (bond) distances:\n";
        cout<<"Max "<<Max_node_pair_length<<"\tMin "<<Min_node_pair_length<<"\tAverage "<<Average_node_pair_length<<endl;
    }
    

    if ((Min_node_pair_length*2<Max_node_pair_length) && Bending_coefficient!=0) {
        if (!Relaxation) {
            cout<<"\nWarning:Initial node distances are not ready/optimised for triangle bending calculations. If the Membrane is not attached to an Actin or another network we strongly recommend turning on the 'Relaxation' flag or switching off the bending.\n\n";
        } else {
            cout<<"\nInitial node distances are not ready/optimised for triangle bending calculations. A few MD steps will be added to the beginning of the simulation to avoid programme break down.\n\n";
            Relaxation_Process_Model=2;
        }
    }
}

void Membrane::check_radius_update_values(void){
    if (Node_radius != New_node_radius && New_node_radius != -1)  {
        if (Begin_update_time_in_Ps == End_update_time_in_Ps) {
            cout<<"The beginning and end of the membrane node radius update time are equal. Please set different values for the paprameters.\n";
            exit(EXIT_FAILURE);
        }
        if (Begin_update_time_in_Ps > End_update_time_in_Ps) {
            cout<<"The beginning and end time of the membrane node radius update times are not in chronological order.\n";
            exit(EXIT_FAILURE);
        }
    }
}


void Membrane::calculate_mesh_properties(void){
    double temp_min=1000;
    double temp_max=0;
    double temp_avg=0;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        double dist=0;
        dist=sqrt((Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])*(Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0]) + (Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])*(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1]) + (Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2])*(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2]));
        temp_avg+=dist;
        
        if (temp_min>dist) {
            temp_min=dist;
        }
        if (dist>temp_max) {
            temp_max=dist;
        }
    }
    temp_avg/=Num_of_Node_Pairs;
//    Max_node_pair_length=temp_max;
//    Min_node_pair_length=temp_min;
    Average_node_pair_length=temp_avg;
    cout<<"Distance properties after relaxation:\n";
    cout<<"Max node distance="<<temp_max<<"\tMin node distance="<<temp_min<<endl;
    cout<<"parameters set on:\n";
    cout<<"Max_node_pair_length="<<Max_node_pair_length<<"\tMin_node_pair_length="<<Min_node_pair_length<<endl;
}

double Membrane::Average_velocity_squared(){
    double average_velocity_squard=0;
    for (int i=0; i< Num_of_Nodes; i++){
        average_velocity_squard+=sqrt(Node_Velocity[i][0]*Node_Velocity[i][0] +Node_Velocity[i][1]*Node_Velocity[i][1] + Node_Velocity[i][2]*Node_Velocity[i][2] );
    }
    average_velocity_squard= average_velocity_squard/Num_of_Nodes;
    return(average_velocity_squard);
}


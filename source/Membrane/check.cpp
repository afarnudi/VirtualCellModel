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
    
    for (int i=0; i<(Num_of_Node_Pairs); i++) {
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
    Average_node_pair_length/=(Num_of_Node_Pairs);
    if (!GenConst::Testmode) {
        cout<<"Node pair (bond) distances:\n";
        cout<<"\tMax "<<Max_node_pair_length<<"\tMin "<<Min_node_pair_length<<"\tAverage "<<Average_node_pair_length<<endl;
    }
    
    if ((Min_node_pair_length*2<Max_node_pair_length) && Bending_coefficient!=0) {
        cout<<TWWARN<<"\n!!!Warning"<<TRESET<<",Initial node distances are not ready/optimised for triangle bending calculations. You will need to make sure the system is relaxed to avoid programme break down.\n\n";

    }
}

void Membrane::check_radius_update_values(void){
    
    
    
    if (abs(Radius - New_Radius)>0.1)  {
        if (Begin_update_time_in_Ps == End_update_time_in_Ps) {
            cout<<TWWARN<<"Warning!!!"<<TRESET<<"The beginning and end of the membrane radius update time are equal. Please set different values for the paprameters.\n";
            exit(EXIT_FAILURE);
        }
        if (Begin_update_time_in_Ps > End_update_time_in_Ps) {
            cout<<TWWARN<<"Warning!!!"<<TRESET<<"The beginning and end time of the membrane radius update times are not in chronological order.\n";
            exit(EXIT_FAILURE);
        }
        update_radius_stat=true;
    } else {
        cout<<TWWARN<<"\nWarning!!!"<<TRESET<<"The initial ("<<Radius<<") and final ("<<New_Radius<<") Membrnae radii are equal."<<TWWARN<<"\nThe simulation will proceed without a radius update."<<TRESET<<endl;
        update_radius_stat=false;
    }
    if(End_update_time_in_Ps>GenConst::Simulation_Time_In_Ps) {
        cout<<TWWARN<<"Warning!!!"<<TRESET<<"The Membrane radius update end time ("<<TWARN<<End_update_time_in_Ps<<TRESET" Ps) exceeds the simulation run time ("<<TWARN<<GenConst::Simulation_Time_In_Ps<<TRESET<<" Ps)."<<endl;
        exit(0);
    }
    
    if (update_radius_stat) {
        double ratio = Radius/New_Radius;
        
        if (!GenConst::Testmode) {
            cout<<TWARN<<"Membrane radius is set to change from "<<TRESET<<Radius<<TWARN<<" to "<<TRESET<<New_Radius<<TWARN<<". The proccess will begin at "<<TRESET<<Begin_update_time_in_Ps<<TWARN<<" Ps and end at "<<TRESET<<End_update_time_in_Ps<<TWARN<<" Ps. The new node radii will scale respectivly";
            if (Node_radius_stat=="Au") {
                cout<<"."<<TRESET<<endl;
            } else{
                cout<<" from "<<TRESET<<Node_radius[0]<<TWARN<<" to "<<TRESET<<Node_radius[0]/ratio<<TWARN<<"."<<TRESET<<endl;
            }
            
        }
    }
    
}

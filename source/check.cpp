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
    for (int i=0; i<Num_of_Node_Pairs; i++) {
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
    Average_node_pair_length/=Num_of_Node_Pairs; cout<<"Max node distance="<<Max_node_pair_length<<"\tmin node distance="<<Min_node_pair_length<<"\tAverage node distance="<<Average_node_pair_length<<endl;
    if ((Min_node_pair_length*2>Max_node_pair_length) && Bending_coefficient!=0) {
        cout<<"Initial node distances are not ready/optimised for triangle bending calculations. A few MD steps will be added to the beginning of the simulation to avoid programme break down.\n";
        
    }
}

void Membrane::node_distance_correction(void){
    for(int MD_Step=0 ;MD_Step<=1000 ; MD_Step++){
        
        MD_Evolution_beginning(GenConst::MD_Time_Step);
        
        
        Elastic_Force_Calculator(0);
        
        MD_Evolution_end(GenConst::MD_Time_Step);
        
        
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
}

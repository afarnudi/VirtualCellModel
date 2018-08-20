//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::check(void){
    Min_node_pair_length=1000, Max_node_pair_length=0, Average_node_pair_length=0;
    for (int i=0; i<Membrane_num_of_Node_Pairs; i++) {
        double dist=0;
        dist=sqrt((Membrane_Node_Position[Membrane_Edges[i][0]][0]-Membrane_Node_Position[Membrane_Edges[i][1]][0])*(Membrane_Node_Position[Membrane_Edges[i][0]][0]-Membrane_Node_Position[Membrane_Edges[i][1]][0])+(Membrane_Node_Position[Membrane_Edges[i][0]][1]-Membrane_Node_Position[Membrane_Edges[i][1]][1])*(Membrane_Node_Position[Membrane_Edges[i][0]][1]-Membrane_Node_Position[Membrane_Edges[i][1]][1])+(Membrane_Node_Position[Membrane_Edges[i][0]][2]-Membrane_Node_Position[Membrane_Edges[i][1]][2])*(Membrane_Node_Position[Membrane_Edges[i][0]][2]-Membrane_Node_Position[Membrane_Edges[i][1]][2]));
        Average_node_pair_length+=dist;
        
        if (Min_node_pair_length>dist) {
            Min_node_pair_length=dist;
        }
        if (dist>Max_node_pair_length) {
            Max_node_pair_length=dist;
        }
    }
    Average_node_pair_length/=Membrane_num_of_Node_Pairs; cout<<"Max="<<Max_node_pair_length<<"\tmin="<<Min_node_pair_length<<"\tAverage="<<Average_node_pair_length<<endl;
}

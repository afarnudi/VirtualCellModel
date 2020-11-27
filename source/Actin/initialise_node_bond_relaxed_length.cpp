//
//  Membrane_node_pair_identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"


void Actin::initialise_node_bond_relaxed_length(void){
    double node_pair_distance[3];
    int node_1,node_2;
    Node_Bond_Nominal_Length_in_Nm.clear();
    Node_Bond_Nominal_Length_in_Nm.resize(Num_of_Node_Pairs);
    
    for(int i=0;i<Num_of_Node_Pairs;i++)
    {
        node_1=Node_Bond_list[i][0];
        node_2=Node_Bond_list[i][1];
        
        node_pair_distance[0]=Node_Position[node_2][0]-Node_Position[node_1][0];
        node_pair_distance[1]=Node_Position[node_2][1]-Node_Position[node_1][0];
        node_pair_distance[2]=Node_Position[node_2][2]-Node_Position[node_1][0];
        
        Node_Bond_Nominal_Length_in_Nm[i]=vector_length(node_pair_distance);
    }
}



void Actin::initialise_abp_bond_relaxed_length(void){
    
    double delta_x, delta_y, delta_z;
    int Node_A, Node_B;
    
    for(int i=0;i<Num_of_abp_Pairs;i++)
    {
        Node_A=abp_Bond_list[i][0];
        Node_B=abp_Bond_list[i][1];
        
        delta_x = Node_Position[Node_B][0] - Node_Position[Node_A][0];
        delta_y = Node_Position[Node_B][1] - Node_Position[Node_A][1];
        delta_z = Node_Position[Node_B][2] - Node_Position[Node_A][2];
        
        double a[3]={delta_x, delta_y, delta_z};
        abp_Bond_relaxed_length.push_back( vector_length(a) );
    }
}





void Actin::initialise_MT_bond_relaxed_length(void){
    
    double delta_x, delta_y, delta_z;
    int Node_A, Node_B;
    
    for(int i=0;i<Num_of_MT_Pairs;i++)
    {
        Node_A=MT_Bond_list[i][0];
        Node_B=MT_Bond_list[i][1];
        
        delta_x = Node_Position[Node_B][0] - Node_Position[Node_A][0];
        delta_y = Node_Position[Node_B][1] - Node_Position[Node_A][1];
        delta_z = Node_Position[Node_B][2] - Node_Position[Node_A][2];
        
        double a[3]={delta_x, delta_y, delta_z};
        MT_Bond_relaxed_length.push_back( vector_length(a) );
    }
}


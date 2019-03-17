//
//  Membrane_node_pair_identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"

void Actin::initialise_node_bond_relaxed_length(void){
    
    double delta_x, delta_y, delta_z;
    int Node_A, Node_B;
    
    for(int i=0;i<Num_of_Node_Pairs;i++)
    {
        Node_A=Node_Bond_list[i][0];
        Node_B=Node_Bond_list[i][1];
        
        delta_x = Node_Position[Node_B][0] - Node_Position[Node_A][0];
        delta_y = Node_Position[Node_B][1] - Node_Position[Node_A][1];
        delta_z = Node_Position[Node_B][2] - Node_Position[Node_A][2];
        
        double a[3]={delta_x, delta_y, delta_z};
        Node_Bond_relaxed_length.push_back( vector_length(a) );
    }
}


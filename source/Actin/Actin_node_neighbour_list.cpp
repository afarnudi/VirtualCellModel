//
//  Membrane_node_neighbour_list.cpp
//  Mem
//
//  Created by Ali Farnudi on 21/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"

void Actin::Node_neighbour_list_constructor(){
    Node_neighbour_list.clear();
    Node_neighbour_list.resize(Num_of_Nodes);
    Node_neighbour_list_respective_bond_index.resize(Num_of_Nodes);
    
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        int node_1 = Node_Bond_list[i][0];
        int node_2 = Node_Bond_list[i][1];
        Node_neighbour_list[ node_1 ].push_back( node_2 );
        Node_neighbour_list[ node_2 ].push_back( node_1 );
        Node_neighbour_list_respective_bond_index[ node_1 ].push_back(i);
        Node_neighbour_list_respective_bond_index[ node_2 ].push_back(i);
    }
}


//
//  Membrane_node_neighbour_list.cpp
//  Mem
//
//  Created by Ali Farnudi on 21/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::node_neighbour_list_constructor(){
    
    node_neighbour_list.resize(Membrane_num_of_Nodes);
    for (int i=0; i<Membrane_num_of_Node_Pairs; i++) {
        node_neighbour_list[Membrane_Edges[i][0]].push_back(Membrane_Edges[i][1]);
        node_neighbour_list[Membrane_Edges[i][1]].push_back(Membrane_Edges[i][0]);
    }
}

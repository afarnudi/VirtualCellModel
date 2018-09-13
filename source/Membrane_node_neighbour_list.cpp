//
//  Membrane_node_neighbour_list.cpp
//  Mem
//
//  Created by Ali Farnudi on 21/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::Node_neighbour_list_constructor(){
    
    Node_neighbour_list.resize(Num_of_Nodes);
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        Node_neighbour_list[Node_Bond_list[i][0]].push_back(Node_Bond_list[i][1]);
        Node_neighbour_list[Node_Bond_list[i][1]].push_back(Node_Bond_list[i][0]);
    }
}

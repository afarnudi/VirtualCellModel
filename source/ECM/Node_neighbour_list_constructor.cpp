//
//  ECM_Node_Pair_Identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 18/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.h"

void ECM::Node_neighbour_list_constructor(void){
    
    int Node_A, Node_B;
    Node_neighbour_list.resize(Num_of_Nodes);
    
    
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        
        Node_A=Node_Bond_list[i][0];
        Node_B=Node_Bond_list[i][1];
        
        Node_neighbour_list[Node_A].push_back(Node_B);
        Node_neighbour_list[Node_B].push_back(Node_A);
    }
}


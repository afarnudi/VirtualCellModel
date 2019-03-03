//
//  interaction.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

void Actin_Membrane_shared_Node_Identifier(Actin &act, Membrane mem){
    
    double mem_x, mem_y, mem_z, act_x, act_y, act_z;
    
    vector<int> push;
    push.resize(2);
    
    for (int mem_node_counter=0; mem_node_counter<mem.return_num_of_nodes(); mem_node_counter++) {
        for (int act_node_counter=0; act_node_counter<act.return_num_of_nodes(); act_node_counter++) {
            
            mem_x = mem.return_node_position(mem_node_counter, 0);
            mem_y = mem.return_node_position(mem_node_counter, 1);
            mem_z = mem.return_node_position(mem_node_counter, 2);
            
            act_x = act.return_node_position(mem_node_counter, 0);
            act_y = act.return_node_position(mem_node_counter, 1);
            act_z = act.return_node_position(mem_node_counter, 2);
            
            
            if( mem_x == act_x &&  mem_y == act_y &&  mem_z == act_z   )
            {
                push[0]=act_node_counter;
                push[1]=mem_node_counter;
                act.Actin_Membrane_shared_Node_list.push_back(push);
                break;
            }
        }
    }
    act.Num_of_Actin_Membrane_shared_Nodes=int(act.Actin_Membrane_shared_Node_list.size());
}

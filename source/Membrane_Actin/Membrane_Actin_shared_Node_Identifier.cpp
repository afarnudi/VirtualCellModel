//
//  interaction.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

void Actin_Membrane_shared_Node_Identifier(Actin &act, Membrane mem, int i, int j){
    
    double mem_x, mem_y, mem_z, act_x, act_y, act_z;
    
    vector<int> push;
    vector<vector<int> > Act_Mem_shared_Node_list;
    push.resize(2);
    
    for (int mem_node_counter=0; mem_node_counter<mem.get_num_of_nodes(); mem_node_counter++) {
        for (int act_node_counter=0; act_node_counter<act.get_num_of_nodes(); act_node_counter++) {
            
            mem_x = mem.get_node_position(mem_node_counter, 0);
            mem_y = mem.get_node_position(mem_node_counter, 1);
            mem_z = mem.get_node_position(mem_node_counter, 2);
            
            act_x = act.get_node_position(act_node_counter, 0);
            act_y = act.get_node_position(act_node_counter, 1);
            act_z = act.get_node_position(act_node_counter, 2);
            
            
            if( mem_x == act_x &&  mem_y == act_y &&  mem_z == act_z   )
            {
                push[0]=act_node_counter;
                push[1]=mem_node_counter;
                Act_Mem_shared_Node_list.push_back(push);
                
                break;
                
            }
        }
    }
    if (act.Actin_Membrane_shared_Node_list.size() == j) {
        act.Actin_Membrane_shared_Node_list.push_back(Act_Mem_shared_Node_list);
    }
    act.Num_of_Actin_Membrane_shared_Nodes.push_back(int(Act_Mem_shared_Node_list.size()));
    
    cout<<"\n# of shared nodes between actin " << i << " and membrane " << j << " = "<<act.Num_of_Actin_Membrane_shared_Nodes[j]<<endl;
    if (act.Num_of_Actin_Membrane_shared_Nodes[j] != mem.get_num_of_nodes()) {
        cout<<"Not all of the Membrane" << j << " node positions are compatibale with the actin " << i << " nodes. Please adjust the meshes if you wish for a fully attached membrane.\n";
    }
}

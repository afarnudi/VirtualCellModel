//
//  interaction.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

void Actin_Membrane_shared_Node_Force_calculator(Actin &act, Membrane &mem){
    
    double delta_x, delta_y, delta_z, Node_distance=0, force=0;
    int mem_node, act_node;
    
    for (int i=0 ; i< act.return_num_of_actin_membrane_shared_nodes() ; i++){
        
        act_node=act.Actin_Membrane_shared_Node_list[i][0];
        mem_node=act.Actin_Membrane_shared_Node_list[i][1];
        
        delta_x = mem.return_node_position(mem_node, 0) - act.return_node_position(act_node, 0);
        delta_y = mem.return_node_position(mem_node, 1) - act.return_node_position(act_node, 1);
        delta_z = mem.return_node_position(mem_node, 2) - act.return_node_position(act_node, 2);
        
        double a[3]={delta_x, delta_y, delta_z};
        Node_distance = vector_length(a);
        
        if (Node_distance > 0.001 || Node_distance < -0.001) {
            
            force = -GenConst::Actin_Membrane_Bond_Coefficient; //The force = -GenConst::Actin_Membrane_Bond_Coefficient*Node_distance and for the components it has to get multiplied by delta_xyz/Node_distance. We have ignored the multiplication and division to reduce computation time and round up errors.
            
            mem.add_to_force(force*delta_x, mem_node, 0);
            mem.add_to_force(force*delta_y, mem_node, 1);
            mem.add_to_force(force*delta_z, mem_node, 2);
            
            act.add_to_force(-force*delta_x, mem_node, 0);
            act.add_to_force(-force*delta_y, mem_node, 1);
            act.add_to_force(-force*delta_z, mem_node, 2);
            
        }
        
    }
    
}

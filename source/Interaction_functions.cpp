//
//  Interaction_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 22/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

double return_ecm_membrane_node_distance(Membrane mem, int mem_node, ECM ecm, int ecm_node){
    return sqrt( (mem.Node_Position[mem_node][0]-ecm.Node_Position[ecm_node][0])*(mem.Node_Position[mem_node][0]-ecm.Node_Position[ecm_node][0])+(mem.Node_Position[mem_node][1]-ecm.Node_Position[ecm_node][1])*(mem.Node_Position[mem_node][1]-ecm.Node_Position[ecm_node][1])+(mem.Node_Position[mem_node][2]-ecm.Node_Position[ecm_node][2])*(mem.Node_Position[mem_node][2]-ecm.Node_Position[ecm_node][2]) );
}

double return_triangle_membrane_distance(Membrane mem, int mem_node, ECM ecm, int tri_index, double tri_com[3]){
    int node_A=ecm.ECM_triangle_list[tri_index][0];
    int node_B=ecm.ECM_triangle_list[tri_index][1];
    int node_C=ecm.ECM_triangle_list[tri_index][2];
    
    tri_com[0]=(ecm.Node_Position[node_A][0]+ecm.Node_Position[node_B][0]+ecm.Node_Position[node_C][0])/3.0;
    tri_com[1]=(ecm.Node_Position[node_A][1]+ecm.Node_Position[node_B][1]+ecm.Node_Position[node_C][1])/3.0;
    tri_com[2]=(ecm.Node_Position[node_A][2]+ecm.Node_Position[node_B][2]+ecm.Node_Position[node_C][2])/3.0;
    
    return sqrt( (mem.Node_Position[mem_node][0]-tri_com[0])*(mem.Node_Position[mem_node][0]-tri_com[0])+(mem.Node_Position[mem_node][1]-tri_com[1])*(mem.Node_Position[mem_node][1]-tri_com[1])+(mem.Node_Position[mem_node][2]-tri_com[2])*(mem.Node_Position[mem_node][2]-tri_com[2]) );
}


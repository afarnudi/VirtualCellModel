//
//  Interaction_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 22/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

void Membrane_ECM_shared_node_force (ECM &ecm, Membrane &mem){
    
    int mem_nodes=mem.return_num_of_nodes();
    double force=0, temp_potential_energy=0;
    double delta_x=0, delta_y=0, delta_z=0, Node_distance=0;
    double sigma = 0.8;
    double epsilon = 4 * GenConst::K * GenConst::MD_T * mem.return_ECM_interaction_strength();
    
    for (int i=0; i<mem_nodes; i++) {
        if (mem.ECM_Node_neighbour_list[i].size() != 0) {
            int ecm_Node = mem.ECM_Node_neighbour_list[i][0];
            
            delta_x = mem.return_node_position(i, 0) - ecm.return_node_position(ecm_Node, 0);
            delta_y = mem.return_node_position(i, 1) - ecm.return_node_position(ecm_Node, 1);
            delta_z = mem.return_node_position(i, 2) - ecm.return_node_position(ecm_Node, 2);
            
            double a[3]={delta_x, delta_y, delta_z};
            Node_distance = vector_length(a);
            
            double r_1 = (Node_distance)/sigma;
            double r_3 = r_1*r_1*r_1;
            double r_5 = r_3*r_1*r_1;
            
            force = 2*epsilon*( 2.0/r_5 - 1.0/r_3 )/sigma;
            temp_potential_energy = epsilon * ( r_1*1.0/r_5 - r_1*1.0/r_3  );
            
            force/=Node_distance;
            
            mem.add_to_force(force*delta_x, i, 0);
            mem.add_to_force(force*delta_y, i, 1);
            mem.add_to_force(force*delta_z, i, 2);
//            cout<<"got one\n";
//            ecm.add_to_force(-force*delta_x, ecm_Node, 0);
//            ecm.add_to_force(-force*delta_y, ecm_Node, 1);
//            ecm.add_to_force(-force*delta_z, ecm_Node, 2);
        }
    }
}

//double return_ecm_membrane_node_distance(Membrane mem, int mem_node, ECM ecm, int ecm_node){
//    return sqrt( (mem.Node_Position[mem_node][0]-ecm.Node_Position[ecm_node][0])*(mem.Node_Position[mem_node][0]-ecm.Node_Position[ecm_node][0])+(mem.Node_Position[mem_node][1]-ecm.Node_Position[ecm_node][1])*(mem.Node_Position[mem_node][1]-ecm.Node_Position[ecm_node][1])+(mem.Node_Position[mem_node][2]-ecm.Node_Position[ecm_node][2])*(mem.Node_Position[mem_node][2]-ecm.Node_Position[ecm_node][2]) );
//}
//
//double return_triangle_membrane_distance(Membrane mem, int mem_node, ECM ecm, int tri_index, double tri_com[3]){
//    int node_A=ecm.Triangle_List[tri_index][0];
//    int node_B=ecm.Triangle_List[tri_index][1];
//    int node_C=ecm.Triangle_List[tri_index][2];
//    
//    tri_com[0]=(ecm.Node_Position[node_A][0]+ecm.Node_Position[node_B][0]+ecm.Node_Position[node_C][0])/3.0;
//    tri_com[1]=(ecm.Node_Position[node_A][1]+ecm.Node_Position[node_B][1]+ecm.Node_Position[node_C][1])/3.0;
//    tri_com[2]=(ecm.Node_Position[node_A][2]+ecm.Node_Position[node_B][2]+ecm.Node_Position[node_C][2])/3.0;
//    
//    return sqrt( (mem.Node_Position[mem_node][0]-tri_com[0])*(mem.Node_Position[mem_node][0]-tri_com[0])+(mem.Node_Position[mem_node][1]-tri_com[1])*(mem.Node_Position[mem_node][1]-tri_com[1])+(mem.Node_Position[mem_node][2]-tri_com[2])*(mem.Node_Position[mem_node][2]-tri_com[2]) );
//}


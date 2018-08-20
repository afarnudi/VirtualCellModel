//
//  interaction.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"
double return_ecm_membrane_node_distance(Membrane mem, int mem_node, ECM ecm, int ecm_node){
    return sqrt( (mem.Membrane_Node_Position[mem_node][0]-ecm.ECM_Node_Position[ecm_node][0])*(mem.Membrane_Node_Position[mem_node][0]-ecm.ECM_Node_Position[ecm_node][0])+(mem.Membrane_Node_Position[mem_node][1]-ecm.ECM_Node_Position[ecm_node][1])*(mem.Membrane_Node_Position[mem_node][1]-ecm.ECM_Node_Position[ecm_node][1])+(mem.Membrane_Node_Position[mem_node][2]-ecm.ECM_Node_Position[ecm_node][2])*(mem.Membrane_Node_Position[mem_node][2]-ecm.ECM_Node_Position[ecm_node][2]) );
}

void interaction_1(int MD_Step, Membrane membrane, ECM ecm, vector<int> ECM_membrane_neighbour_list){
    //Build neighbour list
    if (MD_Step% 100 ==0 ) {
        ECM_membrane_neighbour_list.resize(ecm.return_num_of_nodes(),-1);
        for (int ecm_index=0; ecm_index<ecm.return_num_of_nodes(); ecm_index++) {
            
            if (ecm_index%100==0) {
                cout<<"ecm_index="<<ecm_index<<endl;
            }
            
            vector<int> temp_neighbour_list_index;
            vector<double> temp_neighbour_list_distance;
            
            double normal_prob_value=0;
            for (int membrane_index=0; membrane_index<membrane.return_num_of_nodes(); membrane_index++) {
                double temp_ecm_membrane_distance=return_ecm_membrane_node_distance(membrane, membrane_index, ecm, ecm_index);
                
                if (temp_ecm_membrane_distance<ecm.return_interaction_range()) {
                    temp_neighbour_list_index.push_back(membrane_index);
                    temp_neighbour_list_distance.push_back(temp_ecm_membrane_distance);
                    normal_prob_value+=temp_ecm_membrane_distance;
                }
            }
            double random_num=((double) rand() / (RAND_MAX));
            double prob=0;
            for (int i=0; i<temp_neighbour_list_index.size(); i++) {
                prob+=temp_neighbour_list_distance[i]/normal_prob_value;
                if (random_num<prob) {
                    ECM_membrane_neighbour_list[ecm_index]=temp_neighbour_list_index[i];
                    break;
                }
            }
            
        }
    }
    
    //Force Implimintation
    
    
    for (int i=0; i<ecm.return_num_of_nodes(); i++) {
        if (ECM_membrane_neighbour_list[i]!=-1) {
            int mem_index=ECM_membrane_neighbour_list[i];
            double ecm_mem_node_distance=return_ecm_membrane_node_distance(membrane, mem_index, ecm, i);
            double delta_x=membrane.Membrane_Node_Position[mem_index][0]-ecm.ECM_Node_Position[i][0];
            double delta_y=membrane.Membrane_Node_Position[mem_index][1]-ecm.ECM_Node_Position[i][1];
            double delta_z=membrane.Membrane_Node_Position[mem_index][2]-ecm.ECM_Node_Position[i][2];
            double reduced_radius=1.5*ecm.return_epsilon()/ecm_mem_node_distance;
            double reduced_radius_squared=reduced_radius*reduced_radius;
            double reduced_radius_cubed=reduced_radius_squared*reduced_radius;
            double reduced_radius_quint=reduced_radius_cubed*reduced_radius_squared;
            double force_magnitude=4*ecm.return_sigma()*(reduced_radius_quint- reduced_radius_cubed)/(1.5*ecm.return_epsilon()*ecm_mem_node_distance);
            
            ecm.ECM_Node_Force[i][0]+=force_magnitude*delta_x;
            ecm.ECM_Node_Force[i][1]+=force_magnitude*delta_y;
            ecm.ECM_Node_Force[i][2]+=force_magnitude*delta_z;
            
            membrane.Membrane_Node_Force[mem_index][0] += -force_magnitude*delta_x;
            membrane.Membrane_Node_Force[mem_index][1] += -force_magnitude*delta_y;
            membrane.Membrane_Node_Force[mem_index][2] += -force_magnitude*delta_z;
            
        }
    }
}


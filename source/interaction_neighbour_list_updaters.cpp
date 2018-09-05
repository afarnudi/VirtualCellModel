//
//  interaction_neighbour_list_updaters.cpp
//  Mem
//
//  Created by Ali Farnudi on 22/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

void update_neighbour_list(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag){
    
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
                costume_interaction_flag=true;
                //                cout<<"got one\n";
                break;
            }
        }
        
    }
    
    
}

void update_neighbour_list_2(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list){
    
    vector<int> neighbour_pool;
    for (int i=0; i<ECM_membrane_neighbour_list.size(); i++) {
        if (ECM_membrane_neighbour_list[i]!=-1) {
            neighbour_pool.push_back(ECM_membrane_neighbour_list[i]);
        }
    }
    neighbour_pool_constructor(neighbour_pool, membrane);
    
    double range=ecm.return_interaction_range();
    
    for (int ecm_index=0; ecm_index<ecm.return_num_of_nodes(); ecm_index++) {
        
        if (ecm_index%100==0) {
            cout<<"ecm_index="<<ecm_index<<endl;
        }
        
        vector<int> temp_neighbour_list_index;
        vector<double> temp_neighbour_list_distance;
        
        double normal_prob_value=0;
        
        
        for (int membrane_index=0; membrane_index<neighbour_pool.size(); membrane_index++) {
            double temp_ecm_membrane_distance=return_ecm_membrane_node_distance(membrane, neighbour_pool[membrane_index], ecm, ecm_index);
            
            if (temp_ecm_membrane_distance<range) {
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
                
                //                cout<<"got one\n";
                break;
            }
        }
        
    }
    
    
}

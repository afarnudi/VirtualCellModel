//
//  neighbour_pool_constructor.cpp
//  Mem
//
//  Created by Ali Farnudi on 21/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"
//#include "General_constants.h"

void Membrane_ECM_neighbour_finder (ECM ecm, Membrane &mem){
    bool initiate = true;
    int mem_nodes = mem.return_num_of_nodes();
    int ecm_nodes = ecm.return_num_of_nodes();
    for (int i=0; i<mem_nodes; i++) {
        if (mem.ECM_Node_neighbour_list[i].size() != 0) {
            initiate=false;
        }
    }
    
    if (initiate) {
        initialise_ecm_mem_neighbour_list(ecm, mem);
    } else {
//        update_ecm_mem_neighbour_list(ecm, mem);
    }
    
    
    
    
}

void initialise_ecm_mem_neighbour_list (ECM ecm, Membrane &mem){
    int mem_nodes = mem.return_num_of_nodes();
    int ecm_nodes = ecm.return_num_of_nodes();
    double delta_x, delta_y, delta_z, distance;
    double cut_off = mem.return_ECM_interaction_cut_off();
    vector<vector<double> >  dist_list;
    vector<vector<int> >  inde_list;
    dist_list.resize(mem_nodes);
    inde_list.resize(mem_nodes);
    
    for (int i=0; i<mem_nodes; i++) {
//        cout<<"i = "<<i<<endl;
        
        map<double, int> dist;
        map<double, int>::iterator it, it_4, it_end;
        
        for (int j=0; j<ecm_nodes; j++) {
            
            delta_x = mem.return_node_position(i, 0) - ecm.return_node_position(j, 0);
            delta_y = mem.return_node_position(i, 1) - ecm.return_node_position(j, 1);
            delta_z = mem.return_node_position(i, 2) - ecm.return_node_position(j, 2);
            
            double a[3]={delta_x, delta_y, delta_z};
            distance = vector_length(a);
            if (distance <= cut_off) {
                dist[distance] = j;
//                cout<<"dist = "<<distance<<"\tj = "<<j<<endl;
            }
        }
        if (dist.size() > 1) {
            it=dist.begin();
            it_end=dist.end();
            it_end--;
            it_4=it; it_4++; it_4++; it_4++;
            
            while ( it != it_4 ) {
                dist_list[i].push_back( it->first);
                inde_list[i].push_back( it->second);
                it++;
                if (it == it_end) {
                    break;
                }
            }
        } else if (dist.size() == 1){
            it=dist.begin();
            dist_list[i].push_back( it->first);
            inde_list[i].push_back( it->second);
        }
        
    }
//    cout<<"here\n";
    prune_list(mem_nodes, dist_list, inde_list);
    
    add_nodes_to_neighbour_list(mem, dist_list, inde_list);
    
}

void prune_list(int mem_nodes, vector<vector<double> >  &dist_list, vector<vector<int> >  &inde_list){
    
    for (int i=0; i<mem_nodes-1; i++) {
        for (int j=0; j<inde_list[i].size(); j++) {
            
            for (int k=i+1; k<mem_nodes; k++) {
                for (int l=0; l<inde_list[k].size(); l++) {
                    
                    if (inde_list[i][j] == inde_list[k][l]) {
                        
                        if (dist_list[i][j] < dist_list[k][l]) {
                            dist_list[k].erase(dist_list[k].begin()+l);
                            inde_list[k].erase(inde_list[k].begin()+l);
                            break;
                        } else {
                            dist_list[i].erase(dist_list[i].begin()+j);
                            inde_list[i].erase(inde_list[i].begin()+j);
                            j--;
                            break;
                        }
                    }
                    
                }
            }
            
        }
    }
//    for (int i=0; i<mem_nodes; i++) {
//        if (dist_list[i].size() != 0) {
//            cout<< "i= "<<i<<" : dist= "<<dist_list[i][0]<<" -> "<<inde_list[i][0]<<endl;
//        }
//    }
    
}

void add_nodes_to_neighbour_list (Membrane &mem, vector<vector<double> >  dist_list, vector<vector<int> >  inde_list){
    
    int mem_nodes=mem.return_num_of_nodes();
//    double force=0, temp_potential_energy=0;
    
    for (int i=0; i<mem_nodes; i++) {
        if (dist_list[i].size() != 0) {
            mem.ECM_Node_neighbour_list[i].resize(1, inde_list[i][0]);
    
            
        }
    }
}

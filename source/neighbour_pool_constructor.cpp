//
//  neighbour_pool_constructor.cpp
//  Mem
//
//  Created by Ali Farnudi on 21/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//
#include <algorithm>
#include "interaction.hpp"
//#include "General_constants.h"

void Membrane_ECM_neighbour_finder (ECM &ecm, Membrane &mem){
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
        update_ecm_mem_neighbour_list(ecm, mem);
    }
    
}

void initialise_ecm_mem_neighbour_list (ECM &ecm, Membrane &mem){
    int mem_nodes = mem.return_num_of_nodes();
    int ecm_nodes = ecm.return_num_of_nodes();
    
    mem.ECM_Node_neighbour_list.clear();
    mem.ECM_Node_neighbour_list.resize(mem_nodes);
    
    double delta_x, delta_y, delta_z, distance;
    double cut_off = mem.return_ECM_interaction_cut_off();
//    vector<vector<double> >  dist_list;
//    vector<vector<int> >  inde_list;
//
    vector<vector<pair<double, int> > > neighbour_pairs;
//    dist_list.resize(mem_nodes);
//    inde_list.resize(mem_nodes);
    
    for (int i=0; i<mem_nodes; i++) {
//        cout<<"i = "<<i<<endl;
        
//        map<double, int> dist;
//        map<double, int>::iterator it, it_4, it_end;
        neighbour_pairs.resize(i+1);
        for (int j=0; j<ecm_nodes; j++) {
            
            delta_x = mem.return_node_position(i, 0) - ecm.return_node_position(j, 0);
            delta_y = mem.return_node_position(i, 1) - ecm.return_node_position(j, 1);
            delta_z = mem.return_node_position(i, 2) - ecm.return_node_position(j, 2);
//            cout<<delta_x<<"\t"<<delta_y<<"\t"<<delta_z<<endl;
            double a[3]={delta_x, delta_y, delta_z};
            distance = vector_length(a);
//            cout<<cut_off<<endl;
            if (distance < cut_off) {
//                cout<<distance<<endl;
                neighbour_pairs[i].push_back(make_pair(distance, j));
//                cout<<"i="<<i<<" : "<<distance<<"\t-> "<<j<<endl;
//                dist[distance] = j;
            }
        }
        
        if (neighbour_pairs[i].size() > 1) {
            sort(neighbour_pairs[i].begin(), neighbour_pairs[i].end());
            if (neighbour_pairs[i].size() > 4) {
                int size=neighbour_pairs[i].size();
                for (int k=0; k< neighbour_pairs[i].size()-4; k++) {
                    neighbour_pairs[i].erase(neighbour_pairs[i].begin()+4);
                }
            }
        }
        
    }
    
    prune_list(mem_nodes, neighbour_pairs);
    
    add_nodes_to_neighbour_list(mem, neighbour_pairs);
    
}

void prune_list(int mem_nodes, vector<vector<pair<double, int> > > neighbour_pairs){
    
    for (int i=0; i<mem_nodes-1; i++) {
        for (int j=0; j<neighbour_pairs[i].size(); j++) {
            
            for (int k=i+1; k<mem_nodes; k++) {
                for (int l=0; l<neighbour_pairs[k].size(); l++) {
                    
                    if (neighbour_pairs[i][j].second == neighbour_pairs[k][l].second) {
                        
                        if (neighbour_pairs[i][j].first < neighbour_pairs[k][l].second) {
                            neighbour_pairs[k].erase(neighbour_pairs[k].begin()+l);
                            break;
                        } else {
                            neighbour_pairs[i].erase(neighbour_pairs[i].begin()+j);
                            j--;
                            break;
                        }
                    }
                    
                }
            }
            
        }
    }

}

void add_nodes_to_neighbour_list (Membrane &mem, vector<vector<pair<double, int> > > neighbour_pairs){
    
    int mem_nodes=mem.return_num_of_nodes();
//    double force=0, temp_potential_energy=0;
    
    for (int i=0; i<mem_nodes; i++) {
        if (neighbour_pairs[i].size() != 0) {
            mem.ECM_Node_neighbour_list[i].resize(1, neighbour_pairs[i][0].second);
        }
    }
}

void update_ecm_mem_neighbour_list (ECM &ecm, Membrane &mem){
    initialise_ecm_mem_neighbour_list(ecm, mem);
//    int mem_nodes = mem.return_num_of_nodes();
//    int ecm_nodes = ecm.return_num_of_nodes();
//
//    mem.ECM_Node_neighbour_list.clear();
//    mem.ECM_Node_neighbour_list.resize(mem_nodes);
//
//    double delta_x, delta_y, delta_z, distance;
//    double cut_off = mem.return_ECM_interaction_cut_off();
//    vector<vector<double> >  dist_list;
//    vector<vector<int> >  inde_list;
//    dist_list.resize(mem_nodes);
//    inde_list.resize(mem_nodes);
//
//    for (int i=0; i<mem_nodes; i++) {
////                cout<<"i = "<<i<<endl;
//
//        map<double, int> dist;
//        map<double, int>::iterator it, it_4, it_end;
//
//        for (int j=0; j<ecm_nodes; j++) {
//
//            delta_x = mem.return_node_position(i, 0) - ecm.return_node_position(j, 0);
//            delta_y = mem.return_node_position(i, 1) - ecm.return_node_position(j, 1);
//            delta_z = mem.return_node_position(i, 2) - ecm.return_node_position(j, 2);
//
//            double a[3]={delta_x, delta_y, delta_z};
//            distance = vector_length(a);
//            if (distance <= cut_off) {
//                dist[distance] = j;
////                cout<<"i="<<i<<": "<<distance<<"\t-> "<<j<<endl;
//            }
//        }
//
//        if (dist.size() == 1){
//            it=dist.begin();
//            dist_list[i].push_back( it->first);
//            inde_list[i].push_back( it->second);
//        } else if (dist.size() == 2){
//            it=dist.begin();
//            dist_list[i].push_back( it->first);
//            inde_list[i].push_back( it->second);
//            it++;
//            dist_list[i].push_back( it->first);
//            inde_list[i].push_back( it->second);
//        } else if (dist.size() == 3) {
//            it=dist.begin();
//            dist_list[i].push_back( it->first);
//            inde_list[i].push_back( it->second);
//            it++;
//            dist_list[i].push_back( it->first);
//            inde_list[i].push_back( it->second);
//            it++;
//            dist_list[i].push_back( it->first);
//            inde_list[i].push_back( it->second);
//        } else if (dist.size() >= 4) {
//            it=dist.begin();
//            it_4=it; it_4++; it_4++; it_4++;
//
//            dist_list[i].push_back( it->first);
//            inde_list[i].push_back( it->second);
////            cout<<it->first<<endl;
//
//            while ( it != it_4 ) {
//                it++;
//                dist_list[i].push_back( it->first);
//                inde_list[i].push_back( it->second);
////                cout<<it->first<<endl;
//
//            }
//        }
//
//    }
////        cout<<"here\n";
//    prune_list(mem_nodes, dist_list, inde_list);
//
//    add_nodes_to_neighbour_list(mem, dist_list, inde_list);
}


void Vesicle_particle_neighbour_finder (Membrane &particle, Membrane &vesicle){
    bool initiate = true;
    int vesicle_nodes = vesicle.return_num_of_nodes();
    int particle_nodes = particle.return_num_of_nodes();
    for (int i=0; i<particle_nodes; i++) {
        if (particle.Vesicle_Node_neighbour_list[i].size() != 0) {
            initiate=false;
        }
    }
    
    if (initiate) {
        initialise_vesicle_particle_neighbour_list(particle, vesicle);
    } else {
        update_particle_vesicle_neighbour_list(particle, vesicle);
    }
    
}


void initialise_vesicle_particle_neighbour_list (Membrane &particle, Membrane &vesicle){
    int vesicle_nodes = vesicle.return_num_of_nodes();
    int particle_nodes = particle.return_num_of_nodes();
    
    particle.Vesicle_Node_neighbour_list.clear();
    particle.Vesicle_Node_neighbour_list.resize(particle_nodes);
    
    double delta_x, delta_y, delta_z, distance;
    double cut_off = 10.0;//particle.Average_node_pair_length;
//    vector<vector<double> >  dist_list;
//    vector<vector<int> >  inde_list;
//
    vector<vector<pair<double, int> > > neighbour_pairs;
//    dist_list.resize(mem_nodes);
//    inde_list.resize(mem_nodes);
    
    for (int i=0; i<particle_nodes; i++) {
//        cout<<"i = "<<i<<endl;
        
//        map<double, int> dist;
//        map<double, int>::iterator it, it_4, it_end;
        neighbour_pairs.resize(i+1);
        for (int j=0; j<vesicle_nodes; j++) {
            
            delta_x = particle.return_node_position(i, 0) - vesicle.return_node_position(j, 0);
            delta_y = particle.return_node_position(i, 1) - vesicle.return_node_position(j, 1);
            delta_z = particle.return_node_position(i, 2) - vesicle.return_node_position(j, 2);
//            cout<<delta_x<<"\t"<<delta_y<<"\t"<<delta_z<<endl;
            double a[3]={delta_x, delta_y, delta_z};
            distance = vector_length(a);
//            cout<<cut_off<<endl;
            if (distance < cut_off) {
//                cout<<distance<<endl;
                neighbour_pairs[i].push_back(make_pair(distance, j));
//                cout<<"i="<<i<<" : "<<distance<<"\t-> "<<j<<endl;
//                dist[distance] = j;
            }
        }
        
        if (neighbour_pairs[i].size() > 1) {
            sort(neighbour_pairs[i].begin(), neighbour_pairs[i].end());
            if (neighbour_pairs[i].size() > 4) {
                int size=neighbour_pairs[i].size();
                for (int k=0; k< neighbour_pairs[i].size()-4; k++) {
                    neighbour_pairs[i].erase(neighbour_pairs[i].begin()+4);
                }
            }
        }
        
    }
    
    prune_list(particle_nodes, neighbour_pairs);
    
    add_nodes_to_particle_neighbour_list(particle, neighbour_pairs);
    
}


void add_nodes_to_particle_neighbour_list (Membrane &particle, vector<vector<pair<double, int> > > neighbour_pairs){
    
    int particle_nodes=particle.return_num_of_nodes();
//    double force=0, temp_potential_energy=0;
    
    for (int i=0; i<particle_nodes; i++) {
        if (neighbour_pairs[i].size() != 0) {
            particle.Vesicle_Node_neighbour_list[i].resize(1, neighbour_pairs[i][0].second);
        }
    }
}

void update_particle_vesicle_neighbour_list (Membrane &particle, Membrane &vesicle){
initialise_vesicle_particle_neighbour_list (particle,vesicle);
//..
}





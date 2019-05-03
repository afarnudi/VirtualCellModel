//
//  interaction.hpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#ifndef interaction_hpp
#define interaction_hpp

#include <stdio.h>
#include "Membrane.h"
#include "ECM.h"
#include "Chromatin.h"
#include "Actin.h"

//void interaction_1(int MD_Step, Membrane &membrane, ECM &ecm, vector<int> &membrane_ECM_neighbour_list, bool &costume_interaction_flag);
//void interaction_2(int MD_Step, Membrane &membrane, ECM &ecm, vector<int> &membrane_ECM_neighbour_list, bool &costume_interaction_flag);
//void interaction_4(int MD_Step, Membrane &membrane, ECM &ecm, vector<int> &membrane_ECM_neighbour_list, bool &costume_interaction_flag);
//
//void interaction_3(Membrane &membrane);
//
//void neighbour_pool_constructor (vector<int> &neighbour_pool, Membrane membrane);
//void update_neighbour_list(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag);
//void update_neighbour_list_2(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list);
//void Node_ecm_Barrier(Membrane &membrane, ECM &ecm, vector<int> ECM_membrane_neighbour_list);
//
//
//void triangle_update_neighbour_list(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag);
//void triangle_update_neighbour_list_2(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list);
//
//void triangle_update_neighbour_list_new(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag);
//void triangle_update_neighbour_list_new_2(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list);

double return_ecm_membrane_node_distance(Membrane mem, int mem_node, ECM ecm, int ecm_node);
double return_triangle_membrane_distance(Membrane mem, int mem_node, ECM ecm, int tri_index, double tri_com[3]);
bool barrier_2(Membrane mem, int mem_index);


//Chromatin-Membrane
void Chromatin_Membrane_neighbour_finder(Chromatin &chromo, Membrane Mem);
void Chromatin_Membrane_hard_sphere(Chromatin &chromo, Membrane &Mem);
//void Chromatin_Membrane_triangle_collision(Chromatin chromo, Membrane Mem);

//Actin-Membrane
void Actin_Membrane_shared_Node_Identifier(Actin &actin, Membrane Mem, int j);
void Actin_Membrane_shared_Node_Force_calculator(Actin &actin, Membrane &Mem, int j);

//ECM-Membrane
void Membrane_ECM_neighbour_finder (ECM &ecm, Membrane &mem);
void initialise_ecm_mem_neighbour_list (ECM &ecm, Membrane &mem);
void update_ecm_mem_neighbour_list (ECM &ecm, Membrane &mem);
void prune_list(int mem_nodes, vector<vector<pair<double, int> > > neighbour_pairs);
void add_nodes_to_neighbour_list (Membrane &mem, vector<vector<pair<double, int> > > neighbour_pairs);
void Membrane_ECM_shared_node_force (ECM &ecm, Membrane &mem);

#endif /* interaction_hpp */

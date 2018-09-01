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
#include "ECM.hpp"

void interaction_1(int MD_Step, Membrane &membrane, ECM &ecm, vector<int> &membrane_ECM_neighbour_list, bool &costume_interaction_flag);
void interaction_2(int MD_Step, Membrane &membrane, ECM &ecm, vector<int> &membrane_ECM_neighbour_list, bool &costume_interaction_flag);
void interaction_4(int MD_Step, Membrane &membrane, ECM &ecm, vector<int> &membrane_ECM_neighbour_list, bool &costume_interaction_flag);

void interaction_3(Membrane &membrane);

void neighbour_pool_constructor (vector<int> &neighbour_pool, Membrane membrane);
void update_neighbour_list(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag);
void update_neighbour_list_2(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list);
void Node_ecm_Barrier(Membrane &membrane, ECM &ecm, vector<int> ECM_membrane_neighbour_list);


void triangle_update_neighbour_list(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag);
void triangle_update_neighbour_list_2(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list);

void triangle_update_neighbour_list_new(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list, bool &costume_interaction_flag);
void triangle_update_neighbour_list_new_2(Membrane membrane, ECM ecm, vector<int> &ECM_membrane_neighbour_list);

double return_ecm_membrane_node_distance(Membrane mem, int mem_node, ECM ecm, int ecm_node);
double return_triangle_membrane_distance(Membrane mem, int mem_node, ECM ecm, int tri_index, double tri_com[3]);
bool barrier_2(Membrane mem, int mem_index);
#endif /* interaction_hpp */

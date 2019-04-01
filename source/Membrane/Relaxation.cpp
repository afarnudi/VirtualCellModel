//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::Relax_1(void){
   // node_distance_correction();
    calculate_mesh_properties();
}
void Membrane::Relax_2(void){
    node_distance_correction();
    calculate_mesh_properties();
}


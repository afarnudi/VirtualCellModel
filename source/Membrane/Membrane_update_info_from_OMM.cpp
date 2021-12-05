//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"
#include "General_constants.h"
#include "General_functions.hpp"

using std::ofstream;



void Membrane::update_info_from_omm(MyAtomInfo atoms[], int atom_count){
       
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Position[i][0] = atoms[atom_count + i].posInNm[0];
        Node_Position[i][1] = atoms[atom_count + i].posInNm[1];
        Node_Position[i][2] = atoms[atom_count + i].posInNm[2];
        
        Node_Velocity[i][0] = atoms[atom_count + i].velocityInNmperPs[0];
        Node_Velocity[i][1] = atoms[atom_count + i].velocityInNmperPs[1];
        Node_Velocity[i][2] = atoms[atom_count + i].velocityInNmperPs[2];
    }
}

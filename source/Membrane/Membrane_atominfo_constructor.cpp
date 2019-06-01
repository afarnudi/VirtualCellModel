//
//  Membrane_spring_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"


void Membrane::my_atom_info_constructor(void){
    my_atom_info* atoms = new my_atom_info[Num_of_Nodes];
    for (int i=0; i<Num_of_Nodes; i++) {
        //         pdb     mass   charge  vdwRad vdwEnergy   gbsaRad gbsaScale                  initPos
        atoms[i]={label.c_str(), Node_Mass, 0,      0,      0,         0,       0,    Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
    }
}

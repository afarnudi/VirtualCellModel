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


//Actin-Membrane
void Actin_Membrane_shared_Node_Identifier(Actin &actin, Membrane Mem, int i, int j);
void Membrnae_Membrane_shared_Node_Identifier(Membrane &mem0, Membrane mem1, int atom_count_0, int atom_count_1);
#endif /* interaction_hpp */

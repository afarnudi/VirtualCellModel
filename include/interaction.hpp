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

#endif /* interaction_hpp */

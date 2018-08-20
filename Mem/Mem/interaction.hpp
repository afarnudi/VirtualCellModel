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

void interaction_1(int MD_Step, Membrane membrane, ECM ecm, vector<int> membrane_ECM_neighbour_list);


#endif /* interaction_hpp */

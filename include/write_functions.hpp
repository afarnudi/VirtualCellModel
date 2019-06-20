//
//  write_functions.hpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#ifndef write_functions_hpp
#define write_functions_hpp

#include <stdio.h>
#include "ECM.h"
#include "Membrane.h"
#include <iomanip>
#include <limits>


void Results (ECM ecm, string label, char* buffer);
void Results (Membrane membrane, string label, char* buffer);
void ParticleGeneratingReport (char* buffer, Membrane particle );
void MembraneGeneratingReport (char* buffer, Membrane membrane );
void EcmGeneratingReport (char* buffer, ECM ecm );
void checkingForce (Membrane membrane, int MD_Step, char* buffer);

#endif /* write_functions_hpp */

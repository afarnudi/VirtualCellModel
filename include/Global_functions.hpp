//
//  General_functions.hpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef Global_functions_hpp
#define Global_functions_hpp

#include <string>


#include "General_constants.h"
#include "Chromatin.h"
#include "Membrane.h"
#include "OpenMM_structs.h"
#include "Interaction_table.hpp"

/// \file

/// \brief  Set Temperature function
/// \param  MD_step An integer argument indicating the MD time step.
/// \param  temperature A double argument setting the system temperature
/// \param  buffer A double argument setting the target system temperature
void set_temperature(int MD_step, double temperature, int buffer);

using std::string;
using std::vector;

void collect_data(MyOpenMMData*,
                  MyAtomInfo atoms[],
                  NonBondInteractionMap intertable,
                  vector<Membrane>  &mems,
                  double timeInPs);

void print_statistics(int num_of_atoms,
                      int num_of_bonds,
                      int num_of_dihedrals,
                      vector<Membrane> &mems,
                      vector<Chromatin> &chromos);
#endif /* Global_functions_hpp */

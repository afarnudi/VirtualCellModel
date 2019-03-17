//
//  General_functions.hpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef Global_functions_hpp
#define Global_functions_hpp
#include "General_constants.h"

/// \file

/// \brief  Set Temperature function
/// \param  MD_step An integer argument indicating the MD time step.
/// \param  temperature A double argument setting the system temperature
/// \param  buffer A double argument setting the target system temperature
void set_temperature(int MD_step, double temperature, int buffer);

#endif /* Global_functions_hpp */

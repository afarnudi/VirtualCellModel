//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#include "Global_functions.hpp"
#include "General_constants.h"
#include <math.h>

///    \fn Temperature
///     \brief a brief Temperature
///
///    a more detailed Temperature descriptions.
///

/** @brief Dummy class used for illustration purposes. Doing something with it.
 
 Detailed description follows here.
 @author X. XYZ, DESY
 @date March 2008
 */

///
/// This function will linearly change the temperature from "Buffer_temperature" to "MD_T" set in the General Constants in ?? MD steps.
///
void set_temperature(int MD_step, double temperature, int buffer){
    if (MD_step < buffer) {
        double slope=temperature/double(buffer);
        GenConst::Buffer_temperature=MD_step*slope;
    } else {
        GenConst::Buffer_temperature=GenConst::MD_T;
    }
}

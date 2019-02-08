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


void set_temperature(int MD_step, double temperature, int buffer){
    if (MD_step < buffer) {
        double slope=temperature/double(buffer);
        GenConst::Buffer_temperature=MD_step*slope;
    } else {
        GenConst::Buffer_temperature=GenConst::MD_T;
    }
}

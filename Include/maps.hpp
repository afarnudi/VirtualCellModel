//
//  General_constants.h
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 26/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef maps_hpp
#define maps_hpp
#include <map>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

void read_general_parameters(map<string, double> general_param_map, string input_file_name, vector<string> &membrane_list);

#endif /* maps_hpp */

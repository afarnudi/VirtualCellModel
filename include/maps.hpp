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
#include <string>
#include <sstream>
#include <iterator>

using namespace std;

void read_general_parameters(string input_file_name, vector<string> &membrane_config_list, vector<string> &chromatin_config_list, vector<string> &actin_config_list, vector<string> &ecm_config_list, vector<string> &pointparticle_config_list);
void set_parameter(map<string, double> &general_param_map, string param_name, double param_value);

void read_interaction_map(vector<vector<int> > &inter_map);

#endif /* maps_hpp */

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

/**The interaction map specifies the class instances that are allowed to interact with one another and the nature of their interaction. The map is a text file that lists the class member instances of the enviroment followed by the interaction specifier that is represented with an integer. example:
 * For a simple cell that contains an outer layer of membrane, an actin network, a nucleus membrane, 1 chromatin, and an ECM substrait the interaction will be:
 * mem0 1
 * mem1 0   1
 * act0 1   1   1
 * ecm0 1   0   0   1
 * chr0 0   1   0   0   1
 *
 * Here we assume that mem0 is the outer membrane. The class instance indecies are set during runtime in the order in which the respective configuration file directory is written in the general configuration file. The list of labels on the left (mem0, mem1, etc) is ignored by the programme and its function is for the user to keep track of the columens. It should be noted that the programme expects to come across a word in each line (which is ignored) so the user must not delete the labels altogether. But the actual label written in the interaction map is up to the user as long as it is declared in a single word, for example 'abcdef1234ghi and not 'abc 12 def'.
 * The programme sets the interaction between the class instances in the following order: Membranes, Actins, ECMs, Chromatins, Point Particles. The order in which these interactions are specified in the map is important.
 */
void read_interaction_map(vector<vector<int> > &inter_map);

void configfile_generator(int status);

#endif /* maps_hpp */

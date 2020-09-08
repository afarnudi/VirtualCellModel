//
//  OpenMM_structs.h
//  Membrae
//
//  Created by Ali Farnudi on 10/06/2019.
//  Copyright © 2019 Ali Farnudi. All rights reserved.
//

#ifndef Configfile_hpp
#define Configfile_hpp

#include <string>
#include <vector>
#include <map>
#include "ConfigfileStructs.hpp"

using namespace std;

void configfile_generator(int status);

vector<string> split_and_check_for_comments(string line);
int checkflag(string line, map<string, int> FLAG, map<string, vector<string> > &FLAGlines, int flag);
string getmapkey(map<string, int> FLAG, int value);


void get_class_numbers(map<string, vector<string> > config_lines);
int count_class_numbers(vector<string> lines);
vector<vector<string> > sort_class_configs(vector<string> lines);
void check_class_label(string classname, string label, map<string, int> label_map);
vector<vector<string> > check_and_apply_inheritance(vector<vector<string> > configs, map<string, int> labels);

map<string, vector<string> > read_configfile(string configfilename);
void parse_genconfig_parameters(vector<string> lines);
void assign_general_parameters(GeneralParameters &defaults);

vector<vector<int> > parse_interactiontable_parameters(vector<string> lines);



#endif /* OpenMM_structs_h */


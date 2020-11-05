//
//  OpenMM_structs.h
//  Membrae
//
//  Created by Ali Farnudi on 10/06/2019.
//  Copyright © 2019 Ali Farnudi. All rights reserved.
//

#ifndef Interaction_table_hpp
#define Interaction_table_hpp

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <bitset>

#include "General_constants.h"
#include "Configfile.hpp"

using namespace std;

class NonBondInteractionMap{
    bool   status = false;
    int table_size=0;
    int rows =1;
    int force_sum=0;
    int ForceGroupCount=0;
    
    
    map<string, int>         inter_config_type;
    map<int, string>         inter_config_type_reverse;
    map<int, string>         inter_type_name;
    vector<vector<int> >     inter_table;
    vector<vector<bool> >    force_report;
    vector<vector<string> >  force_label;
    vector<vector<bool> >    use_max_radius;
    
    map<int, string>         forceGroupLabel;
    map<int, int>            forceGroupIndex;
    
    bool is_interaction(string word);
    bool needs_report(string word);
    bool optimise_radius(string word);
    
public:
    NonBondInteractionMap(vector<string> lines);
    
    int parse_interacton_type(string word);
    
    string get_interacton_type(int row,int col){
        return inter_type_name[inter_table[row][col]];
    }
    void check_force_consistency();
    
    void set_force_label(int row, int col,string new_label){
        force_label[row][col]=new_label;
    }
    
    string get_force_label(int row, int col){
        return force_label[row][col];
    }
    
    bool get_report_status(int row, int col){
        return force_report[row][col];
    }
    
    int setForceGroup(int row, int col);
    int setForceGroup(int row, int col, int chromotype);
    
    int get_ForceGroupCount(){
        return ForceGroupCount;
    }
    
    string get_ForceGroupLabel(int index){
        return forceGroupLabel[index];
    }
    
    int get_ForceGroup(int index);
    
    bool get_radius_optimisation_status(int row, int col);
};


#endif /* OpenMM_structs_h */


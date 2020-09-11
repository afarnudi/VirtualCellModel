//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
#include <iostream>
#include <fstream>

#include "General_constants.h"
#include "Arg_pars.hpp"
#include "cxxopts.hpp"
#include "Configfile.hpp"

using namespace std;
vector<vector<string> > check_and_apply_inheritance(vector<vector<string> > configs, map<string, int> labels){
    for (auto &conflines: configs) {
        vector<string> split = split_and_check_for_comments(conflines[0]);
        if (split.size()>3) {
            string daughter = split[1];
            string mother = split[3];
            check_class_label(split[0],daughter, labels);
            check_class_label(split[0],mother, labels);
                
            vector<string> copy = configs[labels[daughter]];
            configs[labels[daughter]].clear();
            configs[labels[daughter]]= configs[labels[mother]];
            configs[labels[daughter]][0]=copy[0];
            for (int i=1; i<copy.size();i++){
                vector<string> copysplit = split_and_check_for_comments(copy[i]);
                for (int j=0; j<configs[labels[daughter]].size(); j++) {
                    vector<string> daughtersplit = split_and_check_for_comments(configs[labels[daughter]][j]);
                    if (daughtersplit[0]==copysplit[0]) {
                        configs[labels[daughter]].erase(configs[labels[daughter]].begin()+j);
                    }
                }
            }
            
            for (int i=1; i<copy.size();i++){
                configs[labels[daughter]].push_back(copy[i]);
            }
        }
    }
    return configs;
}


void get_class_numbers(map<string, vector<string> > config_lines){
    GenConst::Num_of_Membranes=0;
    GenConst::Num_of_Actins=0;
    GenConst::Num_of_ECMs=0;
    GenConst::Num_of_Chromatins=0;
    
    if (config_lines["-Membrane"].size()!=0) {
        GenConst::Num_of_Membranes=count_class_numbers(config_lines["-Membrane"]);
    }
    if (config_lines["-Membrane"].size()!=0) {
        GenConst::Num_of_Actins=count_class_numbers(config_lines["-Actin"]);
    }
    if (config_lines["-Membrane"].size()!=0) {
        GenConst::Num_of_ECMs=count_class_numbers(config_lines["-ECM"]);
    }
    if (config_lines["-Membrane"].size()!=0) {
        GenConst::Num_of_Chromatins=count_class_numbers(config_lines["-Chromatin"]);
    }
}
int count_class_numbers(vector<string> config_lines){
    int class_count=0;
    for (auto & element : config_lines) {
        auto split = split_and_check_for_comments(element);
        if (split[0][0]=='-') {
            class_count++;
        }
    }
    return class_count;
}


void check_class_label(string classname, string label, map<string, int> label_map){
    if (label_map.count(label)==0) {
        string errorMessage = TWARN;
        errorMessage+="No "+classname+" with label "+label+" referenced in the configuration file.";
        errorMessage+=TRESET;
        throw std::runtime_error(errorMessage);
    }
    
}

int checkflag(string line, map<string, int> FLAG, map<string, vector<string> > &FLAGlines, int flag){
    vector<string> split = split_and_check_for_comments(line);
    if (split[0][0]=='-') {
        map<string, int>::iterator it;
        it = FLAG.find(split[0]);
        if (it != FLAG.end()){
            flag =it->second;
            return flag;
            
        } else {
            return -1;
        }
    }
    return flag;
}

string getmapkey(map<string, int> FLAG, int value){
    for (auto &i : FLAG) {
       if (i.second == value) {
           return i.first;
       }
    }
    return "";
}

vector<string> split_and_check_for_comments(string line){
    char comment='#';
    vector<string> empty;
    if(!line.empty()){
        istringstream iss(line);
        vector<string> words(istream_iterator<string>{iss}, istream_iterator<string>());
        
        
        
        if (words[0][0] == comment) {
            words.clear();
            return words;
        }else{
            for (int i=0; i<words.size(); i++) {
                if (words[i][0]==comment) {
                    words.erase(words.begin()+i,words.end());
                    return words;
                }
            }
        }
        return words;
    }
    return empty;
}


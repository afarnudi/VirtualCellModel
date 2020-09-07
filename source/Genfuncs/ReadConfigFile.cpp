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

map<string, vector<string> > read_configfile(string configfilename){
    FLAGindex flagindex;
    auto FLAG = flagindex.FLAG;
    
    ifstream readconfig(configfilename.c_str());
    vector<string> generalLines;
    map<string, vector<string> > FLAGlines;
    
    int flag =-1;
    
    if (readconfig.is_open()) {
        string line;
        while(getline(readconfig, line)){
            
            if(line.empty()){
                continue;
            }
            
            flag = checkflag(line, FLAG, FLAGlines, flag);
            string key = getmapkey(FLAG, flag);
            switch (flag) {
                case 0:
                    FLAGlines[key].push_back(line);
                    break;
                case 1:
                    FLAGlines[key].push_back(line);
                    break;
                case 2:
                    FLAGlines[key].push_back(line);
                    break;
                case 3:
                    FLAGlines[key].push_back(line);
                    break;
                case 4:
                    FLAGlines[key].push_back(line);
                    break;
                    
                default:
                    continue;
            }
        }
    } else {
        string errorMessage = "Write Error: Could not create '"+configfilename+"'";
        throw std::runtime_error(errorMessage);
    }
    return FLAGlines;
}

void parse_genconfig_parameters(vector<string> lines){
    GeneralParameters values;
    //first line is the key
    lines[0].erase(lines[0].begin(),lines[0].end());
    for (int i=0; i<lines.size(); i++) {
        if(lines[i].size()!=0){
            vector<string> split = split_and_check_for_comments(lines[i]);
            if (split.size()!=0) {
                map<string, vector<string> >::iterator it;
                it = values.GenParams.find(split[0]);
                if (it != values.GenParams.end()) {
                    lines[i].erase(lines[i].begin(),lines[i].begin()+int(split[0].size()) );
                    it->second[0] = lines[i];
                }
            }
        }
        
    }
    assign_general_parameters(values);
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
    
    istringstream iss(line);
    vector<string> split(istream_iterator<string>{iss}, istream_iterator<string>());
    
    
    
    if (split[0][0] == comment) {
        split.clear();
        return split;
    }else{
        for (int i=0; i<split.size(); i++) {
            if (split[i][0]==comment) {
                split.erase(split.begin()+i,split.end());
                return split;
            }
        }
    }
    
    return split;
}
void assign_general_parameters(GeneralParameters &defaults){
    for (auto const& it : defaults.GenParams)
    {
        vector<string> split = split_and_check_for_comments(it.second[0]);
        
        if (it.first == "SimulationTimeInPs") {
            GenConst::Simulation_Time_In_Ps=stod(split[0]);
        } else if (it.first == "ReportIntervalInFs") {
            GenConst::Report_Interval_In_Fs=stod(split[0]);
        } else if (it.first == "StepSizeInFs") {
            GenConst::Step_Size_In_Fs=stod(split[0]);
        } else if (it.first == "MC_step") {
            GenConst::MC_step=stoi(split[0]);
        } else if (it.first == "Mem_fluidity") {
            GenConst::Mem_fluidity=stoi(split[0]);
        } else if (it.first == "SimulationBoxLength") {
            GenConst::Lbox=stod(split[0]);
            if (int(GenConst::Lbox) == 0) {
                GenConst::Periodic_box=true;
            } else {
                GenConst::Periodic_box=false;
            }
        } else if (it.first == "Integrator") {
            int type=-1;
            if (split[0][0]=='V'){
                type=0;
            } else if (split[0][0]=='B'){
                type=1;
            } else if (split[0][0]=='L'){
                type=2;
            }
            GenConst::Integrator_type=type;
        } else if (it.first == "FrictionInPs") {
            GenConst::frictionInPs=stod(split[0]);
        } else if (it.first == "Temperature") {
            GenConst::temperature=stod(split[0]);
        } else if (it.first == "ReportEnergy") {
            bool stat;
            if (split[0]=="true") {
                stat=true;
            } else if (split[0]=="false"){
                stat=false;
            } else {
                string errorMessage = "Configfile parameter parse error: Could not parse  '"+split[0]+"'. use 'true' or 'false'.";
                throw std::runtime_error(errorMessage);
            }
            GenConst::WantEnergy=stat;
        } else if (it.first == "WriteVelocitiesForces") {
            bool stat;
            if (split[0]=="true") {
                stat=true;
            } else if (split[0]=="false"){
                stat=false;
            } else {
                string errorMessage = "Configfile parameter parse error: Could not parse  '"+split[0]+"'. use 'true' or 'false'.";
                throw std::runtime_error(errorMessage);
            }
            GenConst::WantForce=stat;
            GenConst::WriteVelocitiesandForces=stat;
        } else if (it.first == "CMMotionRemover") {
            GenConst::CMMotionRemoverStep=stoi(split[0]);
            if (stoi(split[0])>0) {
                GenConst::CMMotionRemover = true;
            }
        } else if (it.first == "FrictionInPs") {
            GenConst::frictionInPs=stod(split[0]);
        } else if (it.first == "MCBarostatPressure") {
            GenConst::MCBarostatPressure=stod(split[0]);
        } else if (it.first == "MCBarostatTemperature") {
            GenConst::MCBarostatTemperature=stod(split[0]);
        } else if (it.first == "MCBarostatFrequency") {
            GenConst::MCBarostatFrequency=stod(split[0]);
        } else if (it.first == "MCStep") {
            GenConst::MC_step=stoi(split[0]);
        } else if (it.first == "MemFluidity") {
            GenConst::Mem_fluidity=stoi(split[0]);
        } else if (it.first == "PeriodicBoxVector0") {
            GenConst::PeriodicBoxVector0.push_back(stod(split[0]));
            GenConst::PeriodicBoxVector0.push_back(stod(split[1]));
            GenConst::PeriodicBoxVector0.push_back(stod(split[2]));
        } else if (it.first == "PeriodicBoxVector1") {
            GenConst::PeriodicBoxVector1.push_back(stod(split[0]));
            GenConst::PeriodicBoxVector1.push_back(stod(split[1]));
            GenConst::PeriodicBoxVector1.push_back(stod(split[2]));
        } else if (it.first == "PeriodicBoxVector2") {
            GenConst::PeriodicBoxVector2.push_back(stod(split[0]));
            GenConst::PeriodicBoxVector2.push_back(stod(split[1]));
            GenConst::PeriodicBoxVector2.push_back(stod(split[2]));
        } else if (it.first == "OutputFileName") {
            GenConst::trajectory_file_name = split[0];
        }
        
        
    }
}



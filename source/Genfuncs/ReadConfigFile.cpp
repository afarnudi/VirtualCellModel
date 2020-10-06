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
#include "Membrane.h"
#include "Actin.h"
#include "ECM.h"
#include "Chromatin.h"

using namespace std;

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
                case 5:
                    FLAGlines[key].push_back(line);
                    break;
                    
                default:
                    continue;
            }
        }
    } else {
        string errorMessage = TWARN;
        errorMessage +="Read Error: Could not read '"+configfilename+"'";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    return FLAGlines;
}

vector<vector<string>> sort_class_configs(vector<string> lines){
    vector<vector<string>> class_configs;
    int count = count_class_numbers(lines);
    int class_order = -1;
    if (count>1) {
        map<string, int> class_labels;
        vector<vector<string>> temp_configs;
        for (auto &line:lines){
            vector<string> split = split_and_check_for_comments(line, "Class config sorter: ");
            if (split.size()>0) {
                if (split[0][0]=='-') {
                    class_order++;
                    if (split.size()>1) {
                        class_labels[split[1]]=class_order;
                    } else {
                        class_labels[to_string(class_order)]=class_order;
                    }
                    
                    temp_configs.resize(class_order+1);
                    temp_configs[class_order].push_back(line);
                } else {
                    
                    temp_configs[class_order].push_back(line);
                }
            }
        }
        class_configs=check_and_apply_inheritance(temp_configs,class_labels);
    } else {
        class_configs.push_back(lines);
    }
    
    return class_configs;
}


void parse_genconfig_parameters(vector<string> lines){
    GeneralParameters values;
    //first line is the key
    lines.erase(lines.begin()+0);
    for (int i=0; i<lines.size(); i++) {
        if(lines[i].size()!=0){
            vector<string> split = split_and_check_for_comments(lines[i], "General config parser: ");
            if (split.size()!=0) {
                map<string, vector<string> >::iterator it;
                it = values.GenParams.find(split[0]);
                if (it != values.GenParams.end()) {
                    lines[i].erase(lines[i].begin(),lines[i].begin()+int(split[0].size()) );
                    it->second[0] = lines[i];
                } else {
                    cout<<TWARN<<"Note: \""<<TFILE<<split[0]<<TWARN<<"\" is not a GeneralParameter."<<TRESET<<endl;
                    cout<<"If you wish to edit the configfile, exit. If not, press any key to continue."<<endl;
                    getchar();
                    cout<<TRESET;
                }
            }
        }
        
    }
    GenConst::Membrane_label="mem";
    GenConst::Actin_label="act";
    GenConst::ECM_label="ecm";
    GenConst::Chromatin_label="chr";
    GenConst::Checkpoint_path = "Results/Resumes/OpenMM/";
    
    assign_general_parameters(values);
    
}

void assign_general_parameters(GeneralParameters &defaults){
    for (auto const& it : defaults.GenParams)
    {
        vector<string> split = split_and_check_for_comments(it.second[0], "General parameters: "+it.first);
        
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
        } else if (it.first == "WriteBondsPDB") {
            bool stat;
            if (split[0]=="true") {
                stat=true;
            } else if (split[0]=="false"){
                stat=false;
            } else {
                string errorMessage = TWARN;
                errorMessage+="Configfile parameter parse error: Could not parse  '"+split[0]+"'. use 'true' or 'false'.";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
            GenConst::write_bonds_to_PDB=stat;
        } else if (it.first == "ReportEnergy") {
            bool stat;
            if (split[0]=="true") {
                stat=true;
            } else if (split[0]=="false"){
                stat=false;
            } else {
                string errorMessage = TWARN;
                errorMessage+="Configfile parameter parse error: Could not parse  '"+split[0]+"'. use 'true' or 'false'.";
                errorMessage+= TRESET;
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
        } else if (it.first == "ProjectName") {
            GenConst::ProjectName = split[0];
        }
        
        
    }
}


void generate_report(map<string, vector<string> > config_lines, string Reporname){
    ofstream write_report;
    write_report.open(Reporname.c_str());
    
    //First write everything under -GeneralParameters in the config file
    for (auto &i : config_lines["-GeneralParameters"]) {
        write_report<<i<<endl;
    }
    
    write_report<<"#Default values that where used:"<<endl;
    //Now fill in the defaults missing from the configfile
    GeneralParameters values;
    for (auto &param: values.GenParams){
        bool add =true;
        for (auto  &cparam: config_lines["-GeneralParameters"]) {
            vector<string> split = split_and_check_for_comments(cparam, "General parameters Reporter: Default values: ");
            if (split[0] == param.first) {
                add =false;
                break;
            }
        }
        if (add) {
            write_report<<param.first<<" "<<param.second[0]<<endl;
        }
    }
    
    if (GenConst::Num_of_Membranes!=0) {
        write_report<<endl;
        for (auto &i : config_lines["-Membrane"]) {
            write_report<<i<<endl;
        }
        
        write_report<<"#Default values that where used:"<<endl;
        //Now fill in the defaults missing from the configfile
        Membrane mem;
        for (auto &param: mem.get_map()){
            bool add =true;
            for (auto  &cparam: config_lines["-Membrane"]) {
                vector<string> split = split_and_check_for_comments(cparam, "Membrane Reporter: Default values: ");
                if (split[0] == param.first) {
                    add =false;
                    break;
                }
            }
            if (add) {
                write_report<<param.first<<" "<<param.second[0]<<endl;
            }
        }
    }
    if (GenConst::Num_of_Actins!=0) {
        write_report<<endl;
        for (auto &i : config_lines["-Actin"]) {
            write_report<<i<<endl;
        }
        
        write_report<<"#Default values that where used:"<<endl;
        //Now fill in the defaults missing from the configfile
        Actin act;
        for (auto &param: act.get_map()){
            bool add =true;
            for (auto  &cparam: config_lines["-Actin"]) {
                vector<string> split = split_and_check_for_comments(cparam, "Actin Reporter: Default values: ");
                if (split[0] == param.first) {
                    add =false;
                    break;
                }
            }
            if (add) {
                write_report<<param.first<<" "<<param.second[0]<<endl;
            }
        }
    }
    
    if (GenConst::Num_of_ECMs!=0) {
        write_report<<endl;
        for (auto &i : config_lines["-ECM"]) {
            write_report<<i<<endl;
        }
        
        write_report<<"#Default values that where used:"<<endl;
        //Now fill in the defaults missing from the configfile
        ECM ecm;
        for (auto &param: ecm.get_map()){
            bool add =true;
            for (auto  &cparam: config_lines["-ECM"]) {
                vector<string> split = split_and_check_for_comments(cparam, "ECM Reporter: Default values: ");
                if (split[0] == param.first) {
                    add =false;
                    break;
                }
            }
            if (add) {
                write_report<<param.first<<" "<<param.second[0]<<endl;
            }
        }
    }
    
    if (GenConst::Num_of_Chromatins!=0) {
        write_report<<endl;
        for (auto &i : config_lines["-Chromatin"]) {
            write_report<<i<<endl;
        }
        
        write_report<<"#Default values that where used:"<<endl;
        //Now fill in the defaults missing from the configfile
        Chromatin chromo;
        for (auto &param: chromo.get_map()){
            bool add =true;
            for (auto  &cparam: config_lines["-Chromatin"]) {
                vector<string> split = split_and_check_for_comments(cparam, "Chromatin Reporter: Default values: ");
                if (split[0] == param.first) {
                    add =false;
                    break;
                }
            }
            if (add) {
                write_report<<param.first<<" "<<param.second[0]<<endl;
            }
        }
    }
    for (auto &i : config_lines["-InteractionTable"]) {
        write_report<<"AA "<<i<<endl;
    }
    
}

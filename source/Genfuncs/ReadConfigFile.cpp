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
            vector<string> split = split_and_check_for_comments(line, "reader: ");
            if(line.empty()){
                continue;
            } else if (split.size()==0){
                continue;
            }
            
            flag = checkflag(line, FLAG, FLAGlines, flag);
            string key = getmapkey(FLAG, flag);
//            cout<<line<<endl;
//            cout<<flag<<endl;
//            cout<<key<<endl;
            if (flag != -1) {
                FLAGlines[key].push_back(line);
            } else {
                string errorMessage = TWARN;
                errorMessage +="Read Error: Config Parser: No such class as '"+line+"'";
                errorMessage +=TRESET;
                throw std::runtime_error(errorMessage);
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
    //first line is the key
    lines.erase(lines.begin()+0);
    for (int i=0; i<lines.size(); i++) {
        if(lines[i].size()!=0){
            vector<string> split = split_and_check_for_comments(lines[i], "General config parser: ");
            if (split.size()!=0) {
                map<string, vector<string> >::iterator it;
                it = generalParameters.GenParams.find(split[0]);
                if (it != generalParameters.GenParams.end()) {
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
    
    GenConst::Checkpoint_path = "Results/Resumes/OpenMM/";
    
    
    assign_general_parameters();
    general_parameters_consistency_check();
}

void general_parameters_consistency_check(void){
    bool MCAnisoBarostatScalestat=false;
    for (int i=0; i<3; i++) {
        if (generalParameters.MCAnisoBarostatPressure[i]!=0) {
            generalParameters.MCAnisoBarostatOn=true;
        }
        if (generalParameters.MCAnisoBarostatScaleXYZ[i]!=false) {
            MCAnisoBarostatScalestat=true;
        }
    }
    if (generalParameters.MCAnisoBarostatOn && !MCAnisoBarostatScalestat) {
        string errorMessage = TWARN;
        errorMessage +="General parameters consistency check: MCAnisoBarostat Non of the axes are allowed to scale. At least one axis has to be allowed to scale. Please edit the configurations and try again.\n";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    if (generalParameters.MCAnisoBarostatOn && generalParameters.MCBarostatPressure!=0) {
        string errorMessage = TWARN;
        errorMessage +="You cannot switch on the MCAnisoBarostat and MCBarostat simultaneously. Please edit the configurations and try again.\n";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
}

void assign_general_parameters(void){
//    generalParameters.BoltzmannKJpermolkelvin = 0.008314459920816468;
    
    for (auto const& it : generalParameters.GenParams)
    {
        vector<string> split = split_and_check_for_comments(it.second[0], "General parameters: "+it.first);
        
        if (it.first == "SimulationTimeInPs") {
            generalParameters.Simulation_Time_In_Ps=stod(split[0]);
        } else if (it.first == "ReportIntervalInFs") {
            generalParameters.Report_Interval_In_Fs=stod(split[0]);
        } else if (it.first == "StepSizeInFs") {
            generalParameters.Step_Size_In_Fs=stod(split[0]);
        } else if (it.first == "MC_step") {
            GenConst::MC_step=stoi(split[0]);
        } else if (it.first == "Mem_fluidity") {
            GenConst::Mem_fluidity=stoi(split[0]);
        } else if (it.first == "SimulationBoxLength") {
            generalParameters.Simulation_box_length=stod(split[0]);
            if (int(generalParameters.Simulation_box_length) == 0) {
                generalParameters.Periodic_condtion_status=false;
            } else {
                generalParameters.Periodic_condtion_status=true;
            }
        } else if (it.first == "Integrator") {
            generalParameters.Integrator_type="Verlet";
            if (split[0][0]=='B'){
                generalParameters.Integrator_type="Brownian";
            } else if (split[0][0]=='L'){
                generalParameters.Integrator_type="Langevin";
            } else if (split[0][0]=='C'){
                generalParameters.Integrator_type="Custom";
                if (split.size()>1) {
                    generalParameters.customtemperature=stod(split[1]);
                }
            }
        } else if (it.first == "FrictionInPs") {
            generalParameters.frictionInPs=stod(split[0]);
        } else if (it.first == "Temperature") {
            generalParameters.temperature=stod(split[0]);
            if (generalParameters.customtemperature<0) {
                generalParameters.customtemperature=generalParameters.temperature;
            }
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
            generalParameters.WantEnergy=stat;
        } else if (it.first == "WantVelocity") {
            bool stat;
            if (split[0]=="true") {
                stat=true;
            } else if (split[0]=="false"){
                stat=false;
            } else {
                string errorMessage = "Configfile parameter parse error: Could not parse  '"+split[0]+"'. use 'true' or 'false'.";
                throw std::runtime_error(errorMessage);
            }
            generalParameters.WantVelocity=stat;
        } else if (it.first == "Minimise") {
            bool stat;
            if (split[0]=="true") {
                stat=true;
            } else if (split[0]=="false"){
                stat=false;
            } else {
                string errorMessage = "Configfile parameter parse error: Could not parse  '"+split[0]+"'. use 'true' or 'false'.";
                throw std::runtime_error(errorMessage);
            }
            generalParameters.Minimise=stat;
        } else if (it.first == "MinimiseTolerance") {
            generalParameters.MinimiseTolerance=stod(split[0]);
        } else if (it.first == "MinimiseMaxIterations") {
            generalParameters.MinimiseMaxIterations=stoi(split[0]);
        } else if (it.first == "CMMotionRemover") {
            generalParameters.CMMotionRemoverStep=stoi(split[0]);
            if (stoi(split[0])>0) {
                generalParameters.CMMotionRemover = true;
            }
        } else if (it.first == "MCBarostatPressure") {
            generalParameters.MCBarostatPressure=stod(split[0]);
        } else if (it.first == "MCBarostatTemperature") {
            generalParameters.MCBarostatTemperature=stod(split[0]);
        } else if (it.first == "MCBarostatFrequency") {
            generalParameters.MCBarostatFrequency=stod(split[0]);
        } else if (it.first == "MCAnisoBarostatPressure") {
            generalParameters.MCAnisoBarostatPressure.push_back(stod(split[0]));
            generalParameters.MCAnisoBarostatPressure.push_back(stod(split[1]));
            generalParameters.MCAnisoBarostatPressure.push_back(stod(split[2]));
        } else if (it.first == "MCAnisoBarostatTemperature") {
            generalParameters.MCAnisoBarostatTemperature=stod(split[0]);
        } else if (it.first == "MCAnisoBarostatScaleXYZ") {
            for (int i=0; i<3; i++) {
                bool stat;
                if (split[i]=="true") {
                    generalParameters.MCAnisoBarostatScaleXYZ.push_back(true);
                } else if (split[i]=="false"){
                    generalParameters.MCAnisoBarostatScaleXYZ.push_back(false);
                } else {
                    string errorMessage = "Configfile parameter parse error: MCAnisoBarostatPressure: Could not parse  '"+split[i]+"'. use 'true' or 'false'.";
                    throw std::runtime_error(errorMessage);
                }
            }
        } else if (it.first == "MCAnisoBarostatFrequency") {
            generalParameters.MCAnisoBarostatFrequency=stod(split[0]);
        } else if (it.first == "MCStep") {
            GenConst::MC_step=stoi(split[0]);
        } else if (it.first == "MemFluidity") {
            GenConst::Mem_fluidity=stoi(split[0]);
        } else if (it.first == "PeriodicBoxVector0") {
            generalParameters.PeriodicBoxVector0.push_back(stod(split[0]));
            generalParameters.PeriodicBoxVector0.push_back(stod(split[1]));
            generalParameters.PeriodicBoxVector0.push_back(stod(split[2]));
        } else if (it.first == "PeriodicBoxVector1") {
            generalParameters.PeriodicBoxVector1.push_back(stod(split[0]));
            generalParameters.PeriodicBoxVector1.push_back(stod(split[1]));
            generalParameters.PeriodicBoxVector1.push_back(stod(split[2]));
        } else if (it.first == "PeriodicBoxVector2") {
            generalParameters.PeriodicBoxVector2.push_back(stod(split[0]));
            generalParameters.PeriodicBoxVector2.push_back(stod(split[1]));
            generalParameters.PeriodicBoxVector2.push_back(stod(split[2]));
        } else if (it.first == "ProjectName") {
            generalParameters.ProjectName = split[0];
        } else if (it.first == "MinimumForceDecleration") {
            
            if (split[0]=="true") {
                generalParameters.MinimumForceDecleration = true;
            } else if (split[0]=="false"){
                generalParameters.MinimumForceDecleration = false;
            } else {
                string errorMessage = "Configfile parameter parse error: MinimumForceDecleration: Could not parse  '"+split[0]+"'. use 'true' or 'false'.";
                throw std::runtime_error(errorMessage);
            }
        
        }
        
        
    }
}


void generate_report(map<string, vector<string> > config_lines){
    string Reportname = generalParameters.trajectory_file_name+"_report.txt";
    ofstream write_report;
    write_report.open(Reportname.c_str());
    
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
    
    if (generalParameters.Num_of_Membranes!=0) {
        write_report<<endl;
        vector<vector<string> > membrane_configs = sort_class_configs(config_lines["-Membrane"]);
        for (int ind=0; ind<membrane_configs.size(); ind++) {
            for (auto &i : membrane_configs[ind]) {
                write_report<<i<<endl;
            }
            write_report<<"#Default values that where used:"<<endl;
            //Now fill in the defaults missing from the configfile
            Membrane mem;
            for (auto &param: mem.get_map()){
                bool add =true;
                for (auto  &cparam: membrane_configs[ind]) {
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
            write_report<<endl;
        }
    }
    if (generalParameters.Num_of_Actins!=0) {
        write_report<<endl;
        vector<vector<string> > actin_configs = sort_class_configs(config_lines["-Actin"]);
        for (int ind=0; ind<actin_configs.size(); ind++) {
            for (auto &i : actin_configs[ind]) {
                write_report<<i<<endl;
            }
            write_report<<"#Default values that where used:"<<endl;
            //Now fill in the defaults missing from the configfile
            Actin act;
            for (auto &param: act.get_map()){
                bool add =true;
                for (auto  &cparam: actin_configs[ind]) {
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
            write_report<<endl;
        }
    }
    
    if (generalParameters.Num_of_ECMs!=0) {
        write_report<<endl;
        vector<vector<string> > ecm_configs = sort_class_configs(config_lines["-ECM"]);
        for (int ind=0; ind<ecm_configs.size(); ind++) {
            for (auto &i : ecm_configs[ind]) {
                write_report<<i<<endl;
            }
            write_report<<"#Default values that where used:"<<endl;
            //Now fill in the defaults missing from the configfile
            ECM ecm;
            for (auto &param: ecm.get_map()){
                bool add =true;
                for (auto  &cparam: ecm_configs[ind]) {
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
            write_report<<endl;
        }
    }
    
    if (generalParameters.Num_of_Chromatins!=0) {
        write_report<<endl;
        vector<vector<string> > chromo_configs = sort_class_configs(config_lines["-Chromatin"]);
        for (int ind=0; ind<chromo_configs.size(); ind++) {
            for (auto &i : chromo_configs[ind]) {
                write_report<<i<<endl;
            }
            write_report<<"#Default values that where used:"<<endl;
            //Now fill in the defaults missing from the configfile
            Chromatin chromo;
            for (auto &param: chromo.get_map()){
                bool add =true;
                for (auto  &cparam: chromo_configs[ind]) {
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
            write_report<<endl;
        }
        
    }
    write_report<<endl;
    for (auto &i : config_lines["-InteractionTable"]) {
        write_report<<i<<endl;
    }
    
}

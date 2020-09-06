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
#include "maps.hpp"

using namespace std;

void read_configfile(string configfilename);
vector<string> split_and_check_for_comments(string line);
int checkflag(string line, map<string, int> FLAG, map<string, vector<vector<string> > > &FLAGlines, int flag);
string getmapkey(map<string, int> FLAG, int value);
void parse_genconfig_parameters(vector<string> lines);
void assign_parameters(GeneralParameters &defaults);

ArgStruct_VCM cxxparser_vcm(int argc, char **argv){
    ArgStruct_VCM args;
    string programme_name =argv[0];
    try
    {
        
        cxxopts::Options options(programme_name,
                                 "\n Welcome to VCM.\n"
                                 "To generate a configuration file, use the configuration generator. If you already hava set your configurations you can run the programme or use the 'c' option to input the configuration file."
                                 );
        options
        .positional_help("")
        .show_positional_help();
        
        options
        //                .allow_unrecognised_options()
        .add_options()
        ("h,help", "Print help.")
        ("g", "Configuration generator. 0: generate the template. 1: Generate a customised configuration.", cxxopts::value<int>(),"int")
        ("c",
        "Path to the configuration file.", cxxopts::value<std::string>())
        ;
        
        options.parse_positional({"configfile"});
        
        auto result = options.parse(argc, argv);
        
        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            exit(0);
        }
        
        if (result.count("g"))
        {
            configfile_generator( result["g"].as<int>() );
            exit(0);
        }
        if (result.count("c"))
        {
            read_configfile(result["c"].as<string>() );
            exit(0);
        }
    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(6);
    }
    
    return args;
}

void read_configfile(string configfilename){
    ifstream readconfig(configfilename.c_str());
    
    vector<string> generalLines;
    map<string, int> FLAG;
    FLAG["-GeneralParameters"]=0;
    FLAG["-Membrane"]=1;
    FLAG["-Actin"]=2;
    FLAG["-ECM"]=3;
    FLAG["-Chromatin"]=4;
    
    map<string, vector<vector<string> > > FLAGlines;
    for (auto &i : FLAGlines) {
        i.second.resize(1);
    }
    
    int flag =-1;
    
    if (readconfig.is_open()) {
        
        string line;
        
        
        while(getline(readconfig, line)){
            
            
            if(line.empty()){
                continue;
            }
            flag = checkflag(line, FLAG, FLAGlines, flag);
            string key = getmapkey(FLAG, flag);
            int last = int(FLAGlines[key].size() )-1;
            switch (flag) {
                case 0:
                    
                    FLAGlines[key][last].push_back(line);
                    break;
                    
                default:
                    continue;
            }
        }
//        for (int i=0; i<FLAGlines["-GeneralParameters"][0].size(); i++) {
//            cout<<FLAGlines["-GeneralParameters"][0][i]<<endl;
//        }
        parse_genconfig_parameters(FLAGlines["-GeneralParameters"][0]);
    } else {
        string errorMessage = "Write Error: Could not create '"+configfilename+"'";
        throw std::runtime_error(errorMessage);
    }
}

void parse_genconfig_parameters(vector<string> lines){
    GeneralParameters values;
    //first line is the key
    lines[0].erase(lines[0].begin(),lines[0].end());
    for (int i=0; i<lines.size(); i++) {
        if(lines[i].size()!=0){
            vector<string> split = split_and_check_for_comments(lines[i]);
            if (split.size()!=0) {
                unordered_map<string, vector<string> >::iterator it;
                it = values.GenParams.find(split[0]);
                if (it != values.GenParams.end()) {
                    lines[i].erase(lines[i].begin(),lines[i].begin()+int(split[0].size()) );
                    it->second[0] = lines[i];
//                    cout<<it->first<<" "<<it->second[0]<<endl;
                }
            }
        }
        
    }
    cout<<"here\n";
    assign_parameters(values);
}
int checkflag(string line, map<string, int> FLAG, map<string, vector<vector<string> > > &FLAGlines, int flag){
    
    vector<string> split = split_and_check_for_comments(line);
    if (split[0][0]=='-') {
        map<string, int>::iterator it;
        it = FLAG.find(split[0]);
        if (it != FLAG.end()){
            flag =it->second;
            string key = getmapkey(FLAG, flag);
            FLAGlines[key].resize(FLAGlines[key].size()+1);
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
void assign_parameters(GeneralParameters &defaults){
    for (auto const& it : defaults.GenParams)
    {
        vector<string> split = split_and_check_for_comments(it.second[0]);
        if (it.first == "SimulationTimeInPs") {
            GenConst::Simulation_Time_In_Ps=stod(split[0]);
            cout<<"GenConst::Simulation_Time_In_Ps "<<GenConst::Simulation_Time_In_Ps<<endl;
        } else if (it.first == "ReportIntervalInFs") {
            GenConst::Report_Interval_In_Fs=stod(split[0]);
            cout<<"GenConst::Report_Interval_In_Fs "<<GenConst::Report_Interval_In_Fs<<endl;
        } else if (it.first == "StepSizeInFs") {
            GenConst::Step_Size_In_Fs=stod(split[0]);
            cout<<"GenConst::Step_Size_In_Fs "<<GenConst::Step_Size_In_Fs<<endl;
        } else if (it.first == "MC_step") {
            GenConst::MC_step=stoi(split[0]);
            cout<<"GenConst::MC_step "<<GenConst::MC_step<<endl;
        } else if (it.first == "Mem_fluidity") {
            GenConst::Mem_fluidity=stoi(split[0]);
            cout<<"GenConst::Mem_fluidity "<<GenConst::Mem_fluidity<<endl;
        } else if (it.first == "SimulationBoxLength") {
            GenConst::Lbox=stod(split[0]);
            if (int(GenConst::Lbox) == 0) {
                GenConst::Periodic_box=true;
            } else {
                GenConst::Periodic_box=false;
            }
            cout<<"GenConst::Lbox "<<GenConst::Lbox<<endl;
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
            cout<<"GenConst::Integrator_type "<<GenConst::Integrator_type<<endl;
        } else if (it.first == "FrictionInPs") {
            GenConst::frictionInPs=stod(split[0]);
            cout<<"GenConst::frictionInPs "<<GenConst::frictionInPs<<endl;
        } else if (it.first == "Temperature") {
            GenConst::temperature=stod(split[0]);
            cout<<"GenConst::temperature "<<GenConst::temperature<<endl;
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
            cout<<"GenConst::WantEnergy "<<GenConst::WantEnergy<<endl;
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
            cout<<"GenConst::WantForce "<<GenConst::WantForce<<endl;
            cout<<"GenConst::WriteVelocitiesandForces "<<GenConst::WriteVelocitiesandForces<<endl;
        } else if (it.first == "CMMotionRemover") {
            GenConst::CMMotionRemoverStep=stoi(split[0]);
            if (stoi(split[0])>0) {
                GenConst::CMMotionRemover = true;
            }
            cout<<"GenConst::CMMotionRemover "<<GenConst::CMMotionRemover<<endl;
            cout<<"GenConst::CMMotionRemoverStep "<<GenConst::CMMotionRemoverStep<<endl;
        } else if (it.first == "FrictionInPs") {
            GenConst::frictionInPs=stod(split[0]);
            cout<<"GenConst::frictionInPs "<<GenConst::frictionInPs<<endl;
        }
        
    }
}



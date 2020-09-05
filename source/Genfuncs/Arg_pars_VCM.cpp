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

void configfile_parser(string configfilename);
vector<string> check_for_comments(string line);

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
            configfile_parser(result["c"].as<string>() );
            exit(0);
        }
    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(6);
    }
    
    return args;
}

void configfile_parser(string configfilename){
    ifstream readconfig(configfilename.c_str());
    if (readconfig.is_open()) {
        
        string line;
        int line_num=0;
        
        while(getline(readconfig, line)){
            line_num++;
            
            if(line.empty()){
                continue;
            }
            
            vector<string> split = check_for_comments(line);
            for (int i=0; i<split.size(); i++) {
                cout<<split[i]<<" ";
            }
            cout<<endl;
        }
    } else {
        string errorMessage = "Write Error: Could not create '"+configfilename+"'";
        throw std::runtime_error(errorMessage);
    }
}

vector<string> check_for_comments(string line){
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

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
#include "OpenMM_funcs.hpp"

using namespace std;






ArgStruct_VCM cxxparser_vcm(int argc, char **argv){
    string programme_name =argv[0];
    ArgStruct_VCM userinputs;
    userinputs.configfilename = "None";
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
        ("a,availablePlatforms",
        "generate flags for available platforms.")
        ("platformID", "ID of platform to be used for the simulation. If you want a list of available platforms, use the avalablePlatforms flag.", cxxopts::value<int>(),"int")
        ("platformDeviceID", "ID of platform device to be used for the simulation. If you want a list of available platforms, use the avalablePlatforms flag.", cxxopts::value<int>(),"int")
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
            userinputs.configfilename =result["c"].as<string>() ;
        }
        if (result.count("a"))
        {
            print_platform_info();
            exit(0);
        }
        if (result.count("platformID"))
        {
            userinputs.platforminfo.platform_id =result["platformID"].as<int>() ;
            userinputs.platforminput=true;
        }
        if (result.count("platformDeviceID"))
        {
            userinputs.platforminfo.platform_device_id =result["platformDeviceID"].as<int>() ;
        }
    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(6);
    }
    
    return userinputs;
}





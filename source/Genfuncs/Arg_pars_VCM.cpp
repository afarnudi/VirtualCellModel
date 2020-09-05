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
    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(6);
    }
    
    return args;
}

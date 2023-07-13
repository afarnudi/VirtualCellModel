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
        "Path to the configuration file.", cxxopts::value<std::string>(), "Path/to/the/configuration/file")
        ("r,resume", "'Directory path' of the interupted simulation. The simulation can only resume on the same machine it was originally running on (for a consist run).", cxxopts::value<std::string>(),"Path/to/the/directory")
        ("a,availablePlatforms",
        "generate flags for available platforms.")
        ("platformID", "ID of platform to be used for the simulation. If you want a list of available platforms, use the avalablePlatforms flag.", cxxopts::value<int>(),"int")
        ("platformDeviceID", "ID of platform device to be used for the simulation. If you want a list of available platforms, use the avalablePlatforms flag.", cxxopts::value<int>(),"int")
        ("ommPluginPath", "Path to OpenMM's plugins. Usually located at '/lib/plugins' in openmm's installation path", cxxopts::value<std::string>(), "Path/to/the/openmm/lib/plugins")
//        ("use-voronoi", "with each step multiply forces with the voronoi area of the membrane. Default false", cxxopts::value<bool>()->default_value("false"))
        ("write-at-end", "Write all outputs at the end of simulation. Warning resume will not be supported. Default false", cxxopts::value<bool>()->default_value("false"))
//        ("s,seed", "OpenMM Documentation: \"Set the random number seed. The precise meaning of this parameter is undefined, and is left up to each Platform to interpret in an appropriate way. It is guaranteed that if two simulations are run with different random number seeds, the sequence of random forces will be different. On the other hand, no guarantees are made about the behavior of simulations that use the same seed. In particular, Platforms are permitted to use non-deterministic algorithms which produce different results on successive runs, even if those runs were initialized identically.\"", cxxopts::value<std::string>(),"Default value, 0, randomises the seed")
        
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
        if (result.count("r"))
        {
            userinputs.resumePath =result["r"].as<string>() ;
            userinputs.configfilename = find_resume_config(userinputs.resumePath, generalParameters.Checkpoint_path, generalParameters.Checkpoint_platformName);
            generalParameters.Resume=true;
        }
//        if (result.count("s"))
//        {
//            generalParameters.Seed=result["s"].as<int>();
//        }
        if (result.count("platformID"))
        {
            userinputs.platforminfo.platform_id =result["platformID"].as<int>() ;
            userinputs.platforminput=true;
        }
        if (result.count("platformDeviceID"))
        {
            userinputs.platforminfo.platform_device_id =result["platformDeviceID"].as<int>() ;
        }
//        if (result.count("use-voronoi"))
//        {
//            userinputs.use_voronoi = true;
//        }
        if (result.count("ommPluginPath"))
        {
            userinputs.platforminfo.platformPluginPath =result["ommPluginPath"].as<string>() ;
            userinputs.platformPluginInput=true;
        }
        if (result.count("write-at-end"))
        {
            userinputs.write_at_end = true;
        }
        if (result.count("a"))
        {
            print_platform_info(userinputs);
            exit(0);
        }
    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(6);
    }
    
    return userinputs;
}





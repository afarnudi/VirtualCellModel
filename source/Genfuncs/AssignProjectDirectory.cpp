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


#include <boost/filesystem.hpp>

void assign_project_directories(char* buffer){
    int instanceID=0;
    string trajectory = "Results";
    boost::filesystem::path pTrajectory(trajectory);
    if(!boost::filesystem::exists(pTrajectory)){
        boost::filesystem::create_directory(pTrajectory);
    }
    trajectory +="/"+generalParameters.ProjectName;
    pTrajectory=trajectory;
    if(!boost::filesystem::exists(pTrajectory)){
        boost::filesystem::create_directory(pTrajectory);
    }
    trajectory += "/";
    trajectory += buffer;
    
    string IDtoString = to_string(instanceID);
    string ID = string(3 - IDtoString.length(), '0') + IDtoString;
    pTrajectory=trajectory+"_"+ID;
    while(boost::filesystem::exists(pTrajectory)){
        instanceID++;
        IDtoString = to_string(instanceID);
        ID = string(3 - IDtoString.length(), '0') + IDtoString;
        pTrajectory=trajectory+"_"+ID;
    }
    boost::filesystem::create_directory(pTrajectory);
    
    trajectory+="_"+ID;
    trajectory+="/";
    trajectory+= buffer;
    trajectory+="_"+ID;
    cout<<trajectory<<endl;
    generalParameters.trajectory_file_name=trajectory;
}

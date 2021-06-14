//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
#include <iostream>
#include <fstream>


#include "General_functions.hpp"

using namespace std;


#include <boost/filesystem.hpp>

void assign_project_directories(char* buffer,
                                string projectName,
                                string &trajPath,
                                string &buffPath
                                ){
    int instanceID=0;
    string trajectory = "Results";
    boost::filesystem::path pTrajectory(trajectory);
    if(!boost::filesystem::exists(pTrajectory)){
        boost::filesystem::create_directories(pTrajectory);
    }
    trajectory +="/" + projectName;
    pTrajectory=trajectory;
    if(!boost::filesystem::exists(pTrajectory)){
        boost::filesystem::create_directories(pTrajectory);
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
    boost::filesystem::create_directories(pTrajectory);
    
    trajectory+="_"+ID;
    trajectory+="/";
    pTrajectory=trajectory+"buffs";
    if(!boost::filesystem::exists(pTrajectory)){
        boost::filesystem::create_directories(pTrajectory);
    }
    buffPath = trajectory + "buffs/" + buffer + "_" + ID;
    trajectory+= buffer;
    trajectory+="_"+ID;
    cout<<trajectory<<endl;
    trajPath=trajectory;
}

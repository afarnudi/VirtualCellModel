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
#include "OpenMM_structs.h"
#include <boost/filesystem.hpp>

using namespace std;

void loadCheckpoint(MyOpenMMData* omm, string ckeckpoint_name, string ckeckpointBackup_name, bool &usingBackupCheckpoint){
    cout<<"Loading checkpoint from: "<<ckeckpoint_name<<endl;
    try {
        std::filebuf rfb;
        rfb.open (ckeckpoint_name.c_str(),std::ios::in);
        std::istream rcheckpoint(&rfb);
        omm->context->loadCheckpoint(rcheckpoint);
    } catch (const std::exception& e) {
        try {
            usingBackupCheckpoint=true;
            std::filebuf rfb;
            string backupcheckpoint = ckeckpointBackup_name + "Backup";
            rfb.open (backupcheckpoint.c_str(),std::ios::in);
            std::istream rcheckpoint(&rfb);
            omm->context->loadCheckpoint(rcheckpoint);
        } catch (const std::exception& e) {
            string errorMessage = TWARN;
            errorMessage+="Loading Checkpoint: Both the checkpoint and the backup are curropt. This simulation cannot be resumed.\n";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
    }
}



void saveCheckpoint(MyOpenMMData* omm, string ckeckpoint_name, string ckeckpointBackup_name, bool &usingBackupCheckpoint){
    string ckeckpoint_name_backup = ckeckpointBackup_name + "Backup";
    if (!usingBackupCheckpoint) {
        boost::filesystem::copy_file(ckeckpoint_name, ckeckpoint_name_backup, boost::filesystem::copy_option::overwrite_if_exists);
    } else {
        usingBackupCheckpoint = false;
    }
    
        
    
    std::filebuf wfb;
    wfb.open (ckeckpoint_name.c_str(),std::ios::out);
    std::ostream wcheckpoint(&wfb);
    omm->context->createCheckpoint(wcheckpoint);
    wfb.close();
}

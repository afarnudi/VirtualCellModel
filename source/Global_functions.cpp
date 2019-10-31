//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#include "Global_functions.hpp"
#include "General_constants.h"
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>

///    \fn Temperature
///     \brief a brief Temperature
///
///    a more detailed Temperature descriptions.
///

/** @brief Dummy class used for illustration purposes. Doing something with it.
 
 Detailed description follows here.
 @author X. XYZ, DESY
 @date March 2008
 */

///
/// This function will linearly change the temperature from "Buffer_temperature" to "MD_T" set in the General Constants in ?? MD steps.
///
void set_temperature(int MD_step, double temperature, int buffer){
    if (MD_step < buffer) {
        double slope=temperature/double(buffer);
        GenConst::Buffer_temperature=MD_step*slope;
    } else {
        GenConst::Buffer_temperature=GenConst::MD_T;
    }
}

using std::string;
using std::endl;


void write_velocities(string buffer){
    
    
    string traj_file_name="Results/"+GenConst::trajectory_file_name+buffer+"_vels.txt";
    std::ofstream write_vels;
    write_vels.open(traj_file_name.c_str());
    
    for (int t=0; t<GenConst::velocity_save[0].size(); t++) {
        for (int ind=0; ind<3; ind++) {
            write_vels<<GenConst::vel_times[t]<<"\t";
            for (int n=0; n<GenConst::velocity_save.size(); n++) {
                if (n==GenConst::velocity_save.size()-1) {
                    write_vels<<GenConst::velocity_save[n][t][ind]<<endl;
                } else {
                    write_vels<<GenConst::velocity_save[n][t][ind]<<"\t";
                }
            }
        }
    }
}

void write_CM(string buffer, vector<Chromatin> chromos){
    string cm_file_name="Results/CMs/"+GenConst::trajectory_file_name+buffer;
    
    for (int ch=0; ch< chromos.size(); ch++) {
        string write_name = cm_file_name + std::to_string(ch) + ".txt";
        std::ofstream write_cms;
        write_cms.open(write_name.c_str());
        int N = chromos[ch].get_num_of_nodes();
        for (int x=0; x<N; x++) {
            for (int y=0; y<N; y++) {
                if (y != N-1) {
                    write_cms<<chromos[ch].get_cm(x, y)<<"";
                } else {
                    write_cms<<chromos[ch].get_cm(x, y)<<"\n";
                }
            }
        }
    }
    
}

using std::vector;

void collect_data(MyAtomInfo atoms[],
                  double timeInPs){
    GenConst::velocity_save.resize(6);
    vector<double> temp_vel(3,0);
    vector<int> atom_list ={ 0, 20, 75, 150, 300, 320, 450};
    int counter = 0;
    for (int n : atom_list) {
        for (int ind=0; ind<3; ind++) {
            temp_vel[ind] = atoms[n].velocityInAngperPs[ind]  * OpenMM::AngstromsPerNm;
        }
        GenConst::velocity_save[counter].push_back(temp_vel);
        counter++;
    }
    GenConst::vel_times.push_back(timeInPs);
}


void analysis(string buffer,
              vector<Chromatin> chromos){
    
    write_velocities(buffer);
    write_CM(buffer, chromos);
    
}

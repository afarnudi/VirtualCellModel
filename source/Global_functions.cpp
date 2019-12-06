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

const int EndOfList=-1;
using std::vector;


void write_data(MyAtomInfo atoms[], string buffer);
//void write_CM(string buffer, vector<Chromatin> chromos);
//double calculate_pressure(vector<Membrane> mems);

void collect_data(MyAtomInfo atoms[],
                  string buffer,
                  vector<Chromatin> chromos,
                  vector<Membrane> mems,
                  double timeInPs){
    
    GenConst::data_colection_times.push_back(timeInPs);
    write_data(atoms, buffer);
}


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

void write_data(MyAtomInfo atoms[],
                string buffer){
    
    
    string traj_file_name="Results/"+GenConst::trajectory_file_name+buffer+"_vels_forces.txt";
    std::ofstream wdata;
    wdata.open(traj_file_name.c_str(), std::ios::app);
    
    wdata<<"time: "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1]<<"\tv_x_InAngperPs v_y_InAngperPs v_z_InAngperPs f_x_inN f_y_inN f_z_inN Pressure_in_N/An^2\n";
    for (int t=0; atoms[t].type != EndOfList; t++) {
        wdata<<t<<"\t"<<atoms[t].velocityInAngperPs[0] << "\t" << atoms[t].velocityInAngperPs[1] << "\t" << atoms[t].velocityInAngperPs[2];
        if (GenConst::WantForce) {
            wdata<<"\t"<<atoms[t].force[0] << "\t" << atoms[t].force[1] << "\t" << atoms[t].force[2];
        }
        wdata<<"\n";
    }
}

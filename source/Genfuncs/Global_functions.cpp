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

const int EndOfList=-1;
using std::vector;


//void write_CM(string buffer, vector<Chromatin> chromos);
//double calculate_pressure(vector<Membrane> mems);

void write_data(   MyAtomInfo atoms[]);

void collect_data(MyAtomInfo atoms[],
                  vector<Chromatin> &chromos,
                  vector<Membrane>  &mems,
                  double timeInPs){
    
    GenConst::data_colection_times.push_back(timeInPs);
    
    for (int i=0; i<mems.size(); i++) {
        if (mems[i].get_GeometricProps_flag()) {
            mems[i].write_geometrics();
        }
    }
    
    write_data(atoms);
}


using std::string;
using std::endl;

void write_data(MyAtomInfo atoms[]){
    
    if (GenConst::WantVelocity) {
        string vel_file_name=GenConst::trajectory_file_name+"_vels.txt";
        std::ofstream wvel;
        wvel.open(vel_file_name.c_str(), std::ios::app);
        
        wvel<<"AtomIndex\tVx\tVy\tVz (Nm/Ps)\n";
        wvel<<"timeInPS "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1]<<endl;
        for (int t=0; atoms[t].type != EndOfList; t++) {
            wvel<<t<<"\t"<<atoms[t].velocityInNmperPs[0] << "\t" << atoms[t].velocityInNmperPs[1] << "\t" << atoms[t].velocityInNmperPs[2];
            wvel<<"\n";
        }
    }
}

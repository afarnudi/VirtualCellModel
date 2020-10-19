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
using namespace std;


//void write_CM(string buffer, vector<Chromatin> chromos);
//double calculate_pressure(vector<Membrane> mems);

void write_data(MyOpenMMData* omm,
                MyAtomInfo atoms[],
                NonBondInteractionMap intertable);

void collect_data(MyOpenMMData* omm,
                  MyAtomInfo atoms[],
                  NonBondInteractionMap intertable,
                  vector<Membrane>  &mems,
                  double timeInPs){
    
    GenConst::data_colection_times.push_back(timeInPs);
    
    for (int i=0; i<mems.size(); i++) {
        if (mems[i].get_GeometricProps_flag()) {
            mems[i].write_geometrics();
        }
    }
    
    write_data(omm, atoms, intertable);
}

using OpenMM::Vec3;

void write_data(MyOpenMMData* omm,
                MyAtomInfo atoms[],
                NonBondInteractionMap intertable){
    
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
    if (GenConst::WantForce) {
        int infoMask = 0;
        infoMask += OpenMM::State::Forces;     // for pot. energy (expensive)
        
        const OpenMM::State state = omm->context->getState(infoMask,GenConst::Periodic_box);
        const std::vector<Vec3>& Forces = state.getForces();
        
        ofstream write_total_force;
        string totalforcepath =GenConst::trajectory_file_name+"_Net_Force.txt";
        write_total_force.open(totalforcepath.c_str(),std::ios_base::app);
        write_total_force<<"AtomIndex\tFx\tFy\tFz (KJ/Nm.mol)"<<endl;
        write_total_force<<"timeInPS "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1]<<endl;
        
        for (int t=0; t < (int)Forces.size(); ++t){
            write_total_force<<t<<"\t"<<Forces[t][0] << "\t" << Forces[t][1] << "\t" << Forces[t][2];
            write_total_force<<"\n";
        }
        
        for (int i=1; i<intertable.get_ForceGroupCount()+1; i++) {
            
            const OpenMM::State statei = omm->context->getState(infoMask,GenConst::Periodic_box, intertable.get_ForceGroup(i));
            const std::vector<Vec3>& Forcesi = statei.getForces();
            
            
            ofstream write_force;
            string forcepath =GenConst::trajectory_file_name+"_"+intertable.get_ForceGroupLabel(i)+".txt";
            write_force.open(forcepath.c_str(),std::ios_base::app);
            write_force<<"AtomIndex\tFx\tFy\tFz (KJ/Nm.mol)"<<endl;
            write_force<<"timeInPS "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1]<<endl;
            for (int t=0; t < (int)Forcesi.size(); ++t){
                write_force<<t<<"\t"<<Forcesi[t][0] << "\t" << Forcesi[t][1] << "\t" << Forcesi[t][2];
                write_force<<"\n";
            }
        }
        
        
    }
    
    
}

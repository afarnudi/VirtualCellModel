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
    
    int atom_count=0;
    for (int i=0; i<mems.size(); i++) {
        if (mems[i].get_GeometricProps_flag()) {
            mems[i].update_info_from_omm(atoms, atom_count);
            mems[i].write_geometrics();
        }
        atom_count+=mems[i].get_num_of_nodes();
    }
    
    write_data(omm, atoms, intertable);
}

using OpenMM::Vec3;

void write_data(MyOpenMMData* omm,
                MyAtomInfo atoms[],
                NonBondInteractionMap intertable){
    
//    if (generalParameters.WantVelocity) {
//        string vel_file_name=generalParameters.trajectory_file_name+"_vels.txt";
//        std::ofstream wvel;
//        wvel.open(vel_file_name.c_str(), std::ios::app);
//        
//        wvel<<"AtomIndex\tVx\tVy\tVz (Nm/Ps)\n";
//        wvel<<"timeInPS "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1]<<endl;
//        for (int t=0; atoms[t].type != EndOfList; t++) {
//            wvel<<t<<"\t"<<atoms[t].velocityInNmperPs[0] << "\t" << atoms[t].velocityInNmperPs[1] << "\t" << atoms[t].velocityInNmperPs[2];
//            wvel<<"\n";
//        }
//    }
    
    
    int infoMask = 0;
    infoMask += OpenMM::State::Forces;
    infoMask += OpenMM::State::Energy;      // for pot. energy (expensive)
    infoMask += OpenMM::State::Velocities;  // for kinetic energy (cheapm)
    const OpenMM::State state = omm->context->getState(infoMask,generalParameters.Periodic_condtion_status);
    const std::vector<Vec3>& Forces = state.getForces();
    double potential_energyInKJ =state.getPotentialEnergy();
    double kinetic_energyInKJ =state.getKineticEnergy();
    double potential_energyInKJ_i;
    double timinPs = state.getTime();
    
    
    string kineticEnergyPath =generalParameters.trajectory_file_name+"_Net_kinetic_energy.txt";
    ofstream write_kineticEnergy;
    write_kineticEnergy.open(kineticEnergyPath.c_str(),std::ios_base::app);
    write_kineticEnergy<<timinPs<<" "<<kinetic_energyInKJ<<setprecision(16)<<endl;
    if (isnan(kinetic_energyInKJ)) {
        string errorMessage = TWARN;
        errorMessage +="Output: Total energy of the system is nan.\n";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    if (kinetic_energyInKJ > 1.5*(3/2)*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature*1002) {
        string errorMessage = TWARN;
        errorMessage +="Output: Total energy of the system is larger than 1.5 times the total kinetic energy.\n";
        errorMessage +="time= "+to_string(timinPs)+" E_k = "+to_string(kinetic_energyInKJ)+".\n";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
//    string energyPath =generalParameters.trajectory_file_name+"_Net_energy.txt";
//    ofstream write_energy;
//    write_energy.open(energyPath.c_str(),std::ios_base::app);
//    write_energy<<potential_energyInKJ<<setprecision(16)<<endl;
    if (isnan(potential_energyInKJ)) {
        string errorMessage = TWARN;
        errorMessage +="Output: Total energy of the system is nan.\n";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    
    
    
    if (generalParameters.WantForce) {
        int infoMask = 0;
        infoMask += OpenMM::State::Forces;
        infoMask += OpenMM::State::Energy;      // for pot. energy (expensive)
        infoMask += OpenMM::State::Velocities;  // for kinetic energy (cheapm)
        const OpenMM::State state = omm->context->getState(infoMask,generalParameters.Periodic_condtion_status);
        const std::vector<Vec3>& Forces = state.getForces();
        double potential_energyInKJ =state.getPotentialEnergy();
        double kinetic_energyInKJ =state.getKineticEnergy();
        double potential_energyInKJ_i;
        
        
        string energyPath =generalParameters.trajectory_file_name+"_Net_energy.txt";
        ofstream write_energy;
        write_energy.open(energyPath.c_str(),std::ios_base::app);
        write_energy<<potential_energyInKJ<<setprecision(16)<<endl;
        if (isnan(potential_energyInKJ)) {
            string errorMessage = TWARN;
            errorMessage +="Output: Total energy of the system is nan.\n";
            errorMessage +=TRESET;
            throw std::runtime_error(errorMessage);
        }
        
        
        
        
        ofstream write_total_force;
        string totalforcepath =generalParameters.trajectory_file_name+"_Net_Force.txt";
        write_total_force.open(totalforcepath.c_str(),std::ios_base::app);
        write_total_force<<"AtomIndex\tFx\tFy\tFz (KJ/Nm.mol)"<<endl;
        write_total_force<<"timeInPS "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1];
        write_total_force<<" potential_energy_inKJpermol ";
        write_total_force<<potential_energyInKJ<<setprecision(16)<<endl;
        
        for (int t=0; t < (int)Forces.size(); ++t){
            write_total_force<<t<<"\t"<<Forces[t][0] << "\t" << Forces[t][1] << "\t" << Forces[t][2];
            write_total_force<<"\n";
        }
//        cout<<generalParameters.force_group_count<<endl;
        for (int i=1; i<intertable.get_ForceGroupCount()+1; i++) {
//            cout<<intertable.get_ForceGroup(i)<<endl;
            const OpenMM::State statei = omm->context->getState(infoMask,generalParameters.Periodic_condtion_status, intertable.get_ForceGroup(i));
            const std::vector<Vec3>& Forcesi = statei.getForces();
            potential_energyInKJ_i =statei.getPotentialEnergy();
            
            ofstream write_force;
            string forcepath =generalParameters.trajectory_file_name+"_"+intertable.get_ForceGroupLabel(i)+".txt";
            cout<<forcepath<<endl;
            write_force.open(forcepath.c_str(),std::ios_base::app);
            write_force<<"AtomIndex\tFx\tFy\tFz (KJ/Nm.mol)"<<endl;
            write_force<<"timeInPS "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1];
            write_force<<" potential_energy_inKJpermol ";
            write_force<<potential_energyInKJ_i<<setprecision(16)<<endl;
            for (int t=0; t < (int)Forcesi.size(); ++t){
                write_force<<t<<"\t"<<Forcesi[t][0] << "\t" << Forcesi[t][1] << "\t" << Forcesi[t][2];
                write_force<<"\n";
            }
        }
        
        for (int i=0; i<generalParameters.force_group_count; i++) {
            int forcegroup = 1 << (i+1);
//            cout<<forcegroup<<endl;
            const OpenMM::State statei = omm->context->getState(infoMask,generalParameters.Periodic_condtion_status, forcegroup);
            const std::vector<Vec3>& Forcesi = statei.getForces();
//            cout<<Forcesi[0][0]<<endl;
            potential_energyInKJ_i =statei.getPotentialEnergy();
            
            ofstream write_force;
            
            string forcepath =generalParameters.trajectory_file_name+"_"+generalParameters.force_group_label[i]+".txt";
//            cout<<forcepath<<endl;
            string energyPath =generalParameters.trajectory_file_name+"_"+generalParameters.force_group_label[i]+"energy.txt";
            ofstream write_energy;
            write_energy.open(energyPath.c_str(),std::ios_base::app);
            write_energy<<potential_energyInKJ_i<<setprecision(16)<<endl;
            write_force.open(forcepath.c_str(),std::ios_base::app);
            write_force<<"AtomIndex\tFx\tFy\tFz (KJ/Nm.mol)"<<endl;
            write_force<<"timeInPS "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1];
            write_force<<" potential_energy_inKJpermol ";
            write_force<<potential_energyInKJ_i<<setprecision(16)<<endl;
            for (int t=0; t < (int)Forcesi.size(); ++t){
                write_force<<t<<"\t"<<Forcesi[t][0] << "\t" << Forcesi[t][1] << "\t" << Forcesi[t][2];
                write_force<<"\n";
            }
        }
        
        
        
        
        
    }
    if (generalParameters.Integrator_type=="Custom") {
        int forcegroup = 1 << 31;
        int infoMask = 0;
        infoMask += OpenMM::State::Forces;     // for pot. energy (expensive)
        
        const OpenMM::State statei = omm->context->getState(infoMask,generalParameters.Periodic_condtion_status, forcegroup);
        const std::vector<Vec3>& Forcesi = statei.getForces();
        
        ofstream write_force31;
        string force31path =generalParameters.trajectory_file_name+"_f31.txt";
        write_force31.open(force31path.c_str(),std::ios_base::app);
        write_force31<<"AtomIndex\tFx\tFy\tFz (KJ/Nm.mol)"<<endl;
        write_force31<<"timeInPS "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1]<<endl;
        for (int t=0; t < (int)Forcesi.size(); ++t){
            write_force31<<t<<"\t"<<Forcesi[t][0] << "\t" << Forcesi[t][1] << "\t" << Forcesi[t][2];
            write_force31<<"\n";
        }
    }
    
    
}

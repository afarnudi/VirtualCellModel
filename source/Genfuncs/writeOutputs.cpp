//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
#include <iostream>
#include "write_functions.hpp"


void writeOutputs(int                atom_count,
                  int                frameNum,
                  const MyAtomInfo   all_atoms[],
                  double             time,
                  double             energyInKJ,
                  double             potential_energyInKJ,
                  bool               usingBackupCheckpoint
                  ){
    
    if (generalParameters.simulation_initial_energy<0) {
        generalParameters.simulation_initial_energy=energyInKJ;
    }
    
//    double kinetic_energyInKJ = energyInKJ-potential_energyInKJ;
//    string kineticEnergyPath =generalParameters.trajectory_file_name+"_Net_kinetic_energy.txt";
//    ofstream write_kineticEnergy;
//    write_kineticEnergy.open(kineticEnergyPath.c_str(),std::ios_base::app);
//    write_kineticEnergy<<time<<" "<<kinetic_energyInKJ<<setprecision(16)<<endl;
//    if (isnan(kinetic_energyInKJ)) {
//        string errorMessage = TWARN;
//        errorMessage +="Output: Total energy of the system is nan.\n";
//        errorMessage +=TRESET;
//        throw std::runtime_error(errorMessage);
//    }
//    double kinetic_energy_threshold =1.5*(3./2.)*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature*1002;
    
//    if (kinetic_energyInKJ > kinetic_energy_threshold) {
//        string errorMessage = TWARN;
//        errorMessage +="Output: Total energy of the system is larger than 1.5 times the total kinetic energy.\n";
//        errorMessage +="time= "+to_string(time)+" E_k = "+to_string(kinetic_energyInKJ)+" E_k max = "+to_string(kinetic_energy_threshold)+".\n";
//        errorMessage +=TRESET;
//        throw std::runtime_error(errorMessage);
//    }
//    cout<<1.5*generalParameters.simulation_initial_energy<<endl;exit(0);
    
    
    
    
    
//    if (energyInKJ > 1.5*generalParameters.simulation_initial_energy) {
//        string errorMessage = TWARN;
//        errorMessage +="Output: Total energy of the system is larger than 1.5 times the total kinetic energy.\n";
//        errorMessage +="time= "+to_string(time)+" E_T = "+to_string(energyInKJ)+" E_k max = "+to_string(1.5*generalParameters.simulation_initial_energy)+".\n";
//        errorMessage +=TRESET;
//        throw std::runtime_error(errorMessage);
//    }
    
    if (isnan(energyInKJ)) {
        string errorMessage = TWARN;
        errorMessage +="Output: Total energy of the system is nan.\n";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    string energyPath =generalParameters.trajectory_file_name+"_Net_energy.txt";
    ofstream write_energy;
    write_energy.open(energyPath.c_str(),std::ios_base::app);
    write_energy<<time<<" "<<energyInKJ<<setprecision(16)<<endl;
    
    if (generalParameters.WantPDB) {
        myWritePDBFrame(frameNum, time, energyInKJ, potential_energyInKJ, all_atoms, !usingBackupCheckpoint);
    }
    if (generalParameters.WantXYZ) {
        writeXYZFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ, !usingBackupCheckpoint);
    }
    if (generalParameters.WantVelocity) {
        writeVelocitiesFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ, !usingBackupCheckpoint);
    }
    if (generalParameters.WantXYZbin) {
        writeXYZbinFrame(all_atoms, generalParameters.precision, !usingBackupCheckpoint);
    }
    if (generalParameters.WantVelocityBin) {
        writeVELbinFrame(all_atoms, generalParameters.precision, !usingBackupCheckpoint);
    }
    if (generalParameters.WantTPKBin) {
        writeTPKbinFrame(time, energyInKJ, potential_energyInKJ, generalParameters.precision, !usingBackupCheckpoint);
    }
    
}


void appendOutputs(int                atom_count,
                   int                frameNum,
                   const MyAtomInfo   all_atoms[],
                   double             time,
                   double             energyInKJ,
                   double             potential_energyInKJ,
                   xyzStashInfo      &xyzStash
                   ){
    
    //    if (generalParameters.WantPDB) {
    //        myWritePDBFrame(frameNum, time, energyInKJ, potential_energyInKJ, all_atoms, !usingBackupCheckpoint);
    //    }
    
    if (generalParameters.WantXYZ) {
        xyzStash.xyzframes += stashXYZFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ);
    }
    if (generalParameters.WantVelocity) {
        xyzStash.velframes += stashVelocitiesFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ);
    }
    if (generalParameters.WantXYZbin) {
        xyzStash.xyzbinframes += stashXYZbinFrame(all_atoms, generalParameters.precision);
    }
    if (generalParameters.WantVelocityBin) {
        xyzStash.velbinframes += stashVELbinFrame(all_atoms, generalParameters.precision);
    }
    if (generalParameters.WantTPKBin) {
        xyzStash.tpkbinframes += stashTPKbinFrame(time, energyInKJ, potential_energyInKJ, generalParameters.precision);
    }
    
}

void flushOutputs(xyzStashInfo      xyzStash,
                  int                atom_count,
                                    int                frameNum,
                                    const MyAtomInfo   all_atoms[],
                                    double             time,
                                    double             energyInKJ,
                                    double             potential_energyInKJ){
    
    //    if (generalParameters.WantPDB) {
    //        myWritePDBFrame(frameNum, time, energyInKJ, potential_energyInKJ, all_atoms, !usingBackupCheckpoint);
    //    }
    
    if (generalParameters.WantXYZ) {
        string readline;
        string buff_name= generalParameters.buffer_file_name+".xyz_buff";
        string traj_name= generalParameters.trajectory_file_name+".xyz";
        
        
        
        ifstream readxyzb(buff_name.c_str());
        ofstream writexyz(traj_name.c_str(), ios_base::app);
        if (!readxyzb.is_open()) {
            string errorMessage = TWARN;
            errorMessage+="No trajectory buffer available to read. If this is an immature simulation please run it from the beginning.\n";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        
        
        while (getline(readxyzb,readline)) {
            writexyz<<readline<<"\n";
        }
        readxyzb.close();
        writexyz<<xyzStash.xyzframes;
        
        writexyz.close();
        xyzStash.xyzframes = "";
    }
    if (generalParameters.WantVelocity) {
        string readline;
        string buff_name= generalParameters.buffer_file_name+".vel_buff";
        string traj_name= generalParameters.trajectory_file_name+".vel";
        
        
        
        ifstream readxyzb(buff_name.c_str());
        ofstream writexyz(traj_name.c_str(), ios_base::app);
        if (!readxyzb.is_open()) {
            string errorMessage = TWARN;
            errorMessage+="No trajectory buffer available to read. If this is an immature simulation please run it from the beginning.\n";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        
        
        while (getline(readxyzb,readline)) {
            writexyz<<readline<<"\n";
        }
        readxyzb.close();
        writexyz<<xyzStash.velframes;
        
        writexyz.close();
        xyzStash.velframes = "";
    }
    if (generalParameters.WantXYZbin) {
        
        string readline;
        string buff_name= generalParameters.buffer_file_name+".bin_xyz_"+generalParameters.precision+"_buff";
        string traj_name= generalParameters.trajectory_file_name+".bin_xyz_"+generalParameters.precision;
        
        ofstream writexyz(traj_name.c_str(), std::ios::app );
        ifstream readxyzb(buff_name.c_str() );
        
        if (!readxyzb.is_open()) {
            string errorMessage = TWARN;
            errorMessage+="No trajectory buffer available to read. If this is an immature simulation please run it from the beginning.\n";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        
        writexyz<<readxyzb.rdbuf();
        writexyz<<xyzStash.xyzbinframes;
        
        writexyz.close();
        xyzStash.xyzbinframes="";
    }
    if (generalParameters.WantVelocityBin) {
        
        string readline;
        string buff_name= generalParameters.buffer_file_name+".bin_vel_"+generalParameters.precision+"_buff";
        string traj_name= generalParameters.trajectory_file_name+".bin_vel_"+generalParameters.precision;
        
        ofstream writexyz(traj_name.c_str(), std::ios::app );
        ifstream readxyzb(buff_name.c_str() );
        
        if (!readxyzb.is_open()) {
            string errorMessage = TWARN;
            errorMessage+="No trajectory buffer available to read. If this is an immature simulation please run it from the beginning.\n";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        
        writexyz<<readxyzb.rdbuf();
        writexyz<<xyzStash.velbinframes;
        
        
        writexyz.close();
        xyzStash.velbinframes="";
    }
    if (generalParameters.WantTPKBin) {
        
        string readline;
        string buff_name= generalParameters.buffer_file_name+".bin_tpk_"+generalParameters.precision+"_buff";
        string traj_name= generalParameters.trajectory_file_name+".bin_tpk_"+generalParameters.precision;
        
        ofstream writetpk(traj_name.c_str(), std::ios::app );
        ifstream readtpkb(buff_name.c_str() );
        
        if (!readtpkb.is_open()) {
            string errorMessage = TWARN;
            errorMessage+="No tpk buffer available to read. If this is an immature simulation please run it from the beginning.\n";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        
        writetpk<<readtpkb.rdbuf();
        writetpk<<xyzStash.tpkbinframes;
        
        writetpk.close();
        
        
    }
    
}

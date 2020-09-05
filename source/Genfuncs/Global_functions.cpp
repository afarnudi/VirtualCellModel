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


void write_data(MyAtomInfo atoms[],
                string buffer,
                double volume,
                double area,
                double bending_energy);
//void write_CM(string buffer, vector<Chromatin> chromos);
//double calculate_pressure(vector<Membrane> mems);

void write_data2(   MyAtomInfo atoms[],
                 string buffer,
                 double volume,
                 double area,
                 double bending_energy,
                 vector<double> voronoi_area);

void collect_data(MyAtomInfo atoms[],
                  string buffer,
                  vector<Chromatin> chromos,
                  vector<Membrane> mems,
                  double timeInPs){
    
    GenConst::data_colection_times.push_back(timeInPs);
    
    double mem_volume = 0;
    double mem_surface = 0;
    double bending_energy = 0;
    vector <double> voronoi_area;
    //    double stretch_energy = 0;
    if (mems.size()!=0) {
        mems[0].calculate_volume_and_surface_area();
        if (GenConst::Wantvoronoi) {
            voronoi_area = mems[0].return_voronoi_node_area();
        }
        
        mem_volume  = mems[0].return_volume();
        mem_surface = mems[0].get_surface_area();
        bending_energy = mems[0].calculate_bending_energy();
    } else {
        voronoi_area.resize(1,0);
    }
    //    write_data(atoms, buffer, mem_volume, mem_surface, bending_energy);
    write_data2(atoms, buffer, mem_volume, mem_surface, bending_energy, voronoi_area);
}


using std::string;
using std::endl;

void write_data(MyAtomInfo atoms[],
                string buffer,
                double volume,
                double area,
                double bending_energy){
    
    
    string traj_file_name="Results/"+GenConst::trajectory_file_name+buffer+"_vels_forces.txt";
    std::ofstream wdata;
    wdata.open(traj_file_name.c_str(), std::ios::app);
    
    wdata<<"time: "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1]<<"\tvx, vy, vz(Nm/Ps) fx, fy, fz(KJ/Nm) Volume="<<volume<<" (Nm^3) Area="<<area<<" (Nm^2) bending_energy="<<bending_energy<<" (KJ)\n";
    for (int t=0; atoms[t].type != EndOfList; t++) {
        wdata<<t<<"\t"<<atoms[t].velocityInNmperPs[0] << "\t" << atoms[t].velocityInNmperPs[1] << "\t" << atoms[t].velocityInNmperPs[2];
        if (GenConst::WantForce) {
            wdata<<"\t"<<atoms[t].force[0] << "\t" << atoms[t].force[1] << "\t" << atoms[t].force[2];
        }
        wdata<<"\n";
    }
}

void write_data2(MyAtomInfo atoms[],
                 string buffer,
                 double volume,
                 double area,
                 double bending_energy,
                 vector<double> voronoi_area){
    
    
    string traj_file_name="Results/"+GenConst::trajectory_file_name+buffer+"_properties.txt";
    std::ofstream wdata;
    wdata.open(traj_file_name.c_str(), std::ios::app);
    
    //header
    
    wdata<<"vx, vy, vz(Nm/Ps)";
    if (GenConst::WantForce) {
        wdata<<" fx, fy, fz(KJ/Nm)";
    }
    if (GenConst::Wantvoronoi){
        wdata<<" voronoi_area(Nm^2)";
    }
    wdata<<endl;
    
    wdata<<"time "<<GenConst::data_colection_times[GenConst::data_colection_times.size()-1]<<" Volume_(Nm^3) "<<volume<<"  Area_(Nm^2) "<<area<<" bending_energy_(KJ/mol) "<<bending_energy;
    if (GenConst::MCBarostatFrequency!=0) {
        wdata<<" LboxX_(Nm) "<<GenConst::Lboxdims[0][0]<<" LboxY_(Nm) "<<GenConst::Lboxdims[1][1]<<" LboxZ_(Nm) "<<GenConst::Lboxdims[2][2];
    }
    wdata<<"\n";
    for (int t=0; atoms[t].type != EndOfList; t++) {
        wdata<<t<<"\t"<<atoms[t].velocityInNmperPs[0] << "\t" << atoms[t].velocityInNmperPs[1] << "\t" << atoms[t].velocityInNmperPs[2];
        if (GenConst::WantForce) {
            wdata<<"\t"<<atoms[t].force[0] << "\t" << atoms[t].force[1] << "\t" << atoms[t].force[2];
        }
        if (GenConst::Wantvoronoi) {
            wdata<<"\t"<<voronoi_area[t];
        }
        wdata<<"\n";
    }
}

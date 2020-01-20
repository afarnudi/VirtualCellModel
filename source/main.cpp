/** @file doxygen_example.cpp
 @author Lastname:Firstname:A00123456:cscxxxxx
 @version Revision 1.1
 @brief Illustrates doxygen-style comments for documenting a C++
 program file and the functions in that file.
 @details If you want to add any further detailed description of
 what is in the file, then place it here (after the first statement)
 and it will appear in the detailed description section of the HTML
 output description for the file.
 @date Monday, September 19, 2011
 */

/// \file

#include <stdio.h>
#include <ctime>
#include <sstream>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <random>
#include <string>
#include <math.h>
#include <chrono>

#include "Membrane.h"
#include "Chromatin.h"
#include "Actin.h"
#include "ECM.h"

#include "General_functions.hpp"
#include "write_functions.hpp"
#include "interaction.hpp"
#include "maps.hpp"
#include "Global_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_class_functions.h"


/** -----------------------------------------------------------------------------
 *                           OpenMM-USING CODE
 * -----------------------------------------------------------------------------
 * The OpenMM API is visible only at this point and below. Normally this would
 * be in a separate compilation module; we're including it here for simplicity.
 * -----------------------------------------------------------------------------
 */

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
#pragma warning(disable:4996)   // sprintf is unsafe
#endif

#include "OpenMM.h"


namespace GenConst {
    int MD_num_of_steps;
    double Simulation_Time_In_Ps;
    int MD_traj_save_step;
    double Report_Interval_In_Fs;
    double Step_Size_In_Fs;
    double MD_T;
    double K;
    int MD_thrmo_step;
    int MC_step;
    double Mem_fluidity;
    double Lbox;
    bool Periodic_condtion_status;
    int Num_of_Membranes;
    int Num_of_Chromatins;
    int Num_of_Actins;
    int Num_of_ECMs;
    int Num_of_pointparticles;
    string trajectory_file_name;;
    double Buffer_temperature;
    double Bussi_tau;
    double Actin_Membrane_Bond_Coefficient;
    bool Interaction_map;
    string Interaction_map_file_name;
    bool Excluded_volume_interaction;
    bool OpenMM;
    double sigma_LJ_12_6;
    double epsilon_LJ_12_6;
    string Membrane_label;
    string Actin_label;
    string Chromatin_label;
    string ECM_label;
    int Integrator_type;
    double frictionInPs;
    double temperature;
    bool CreateCheckpoint;
    bool Load_from_checkpoint;
    std::string Checkpoint_path;
    std::string Checkpoint_file_name;


    bool   write_bonds_to_PDB;
    bool   WantEnergy;
    bool   WantForce;
    bool   WriteVelocitiesandForces;
    bool   CMMotionRemover;
    int    CMMotionRemoverStep;



//    std::vector<std::vector<std::vector<double> > > data;
    std::vector<double> data_colection_times;
}



const int EndOfList=-1;


int main(int argc, char **argv)
{

    
    // get the current time.
    time_t t = time(0);
    auto chrono_clock_start = chrono::steady_clock::now();
    auto chrono_sys_clock_start = chrono::system_clock::now();
    
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_time_%H_%M",now);
    
    string general_file_name="general-config.txt";
    cout<<"\nHi!\nPlease enter the path (relative to the binary file) + name of the config file:\nexample:\t../../myconfigfile.txt\n\nPath to configuration file: ";
    cin>>general_file_name;
    clock_t tStart = clock();//Time the programme
    vector<string> membrane_config_list;
    vector<string> chromatin_config_list;
    vector<string> actin_config_list;
    vector<string> ecm_config_list;
    vector<string> pointparticle_config_list;
    
    read_general_parameters(general_file_name, membrane_config_list, chromatin_config_list, actin_config_list, ecm_config_list, pointparticle_config_list);
    
    vector<vector<int> > interaction_map;
    read_interaction_map(interaction_map);
    
    ofstream Trajectory;
    ofstream bondenergycheck1;
    ofstream bondenergycheck2;
    ofstream bondlengthcheck;
    ofstream bondlengthcheck2;
    bondenergycheck1.open("memBond", ios::app);
    bondenergycheck2.open("openmmBond", ios::app);
    bondlengthcheck.open("L_bondlength", ios::app);
    bondlengthcheck2.open("R_bondlength", ios::app);
    string traj_file_name="Results/"+GenConst::trajectory_file_name+buffer+".xyz";
    string ckeckpoint_name=GenConst::Checkpoint_path+GenConst::trajectory_file_name+buffer;
    
    
    vector<Membrane> Membranes;
    vector<std::set<int> > membrane_set;
    
    vector<Actin> Actins;
    vector<std::set<int> > actin_set;
    
    vector<ECM> ECMs;
    vector<std::set<int> > ecm_set;
    
    vector<Chromatin> Chromatins;
//    vector<std::set<int> > chromatin_set;
    vector<vector<std::set<int> > > chromatin_set;
    
    vector<point_particle> pointparticles;
    
    bool Include_Membrane  = false;
    bool Include_Chromatin = false;
    bool Include_Actin     = false;
    bool Include_ECM       = false;
    bool Include_pointparticle = false;
    
    int num_of_atoms=0;
    int num_of_bonds=0;
    int num_of_dihedrals=0;
    
    if (!GenConst::Load_from_checkpoint) {
        if (GenConst::Num_of_Membranes!=0) {
            Include_Membrane = true;
            
            Membranes.resize(GenConst::Num_of_Membranes);
            membrane_set.resize(GenConst::Num_of_Membranes);
            for (int i=0; i<GenConst::Num_of_Membranes; i++) {
                string label=GenConst::Membrane_label+to_string(i);
                Membranes[i].set_label(label);
                Membranes[i].set_file_time(buffer);
                Membranes[i].set_index(i);
                Membranes[i].import_config(membrane_config_list[i]);
                Membranes[i].generate_report();
            }
        }
        
        
        if (GenConst::Num_of_Actins!=0) {
            Include_Actin = true;
            
            Actins.resize(GenConst::Num_of_Actins);
            actin_set.resize(GenConst::Num_of_Actins);
            for (int i=0; i<GenConst::Num_of_Actins; i++) {
                string label=GenConst::Actin_label+to_string(i);
                Actins[i].set_label(label);
                Actins[i].set_file_time(buffer);
                Actins[i].set_index(i);
                Actins[i].import_config(actin_config_list[i]);
                Actins[i].generate_report();
            }
        }
        
        if (GenConst::Num_of_ECMs!=0){
            Include_ECM=true;
            
            ECMs.resize(GenConst::Num_of_ECMs);
            ecm_set.resize(GenConst::Num_of_ECMs);
            for (int i=0; i<GenConst::Num_of_ECMs; i++) {
                string label=GenConst::ECM_label+to_string(i);
                ECMs[i].set_label(label);
                ECMs[i].set_file_time(buffer);
                ECMs[i].set_index(i);
                ECMs[i].import_config(ecm_config_list[i]);
                ECMs[i].generate_report();
            }
            
        }
        
        
        if (GenConst::Num_of_Chromatins!=0) {
            Include_Chromatin = true;
            
            Chromatins.resize(GenConst::Num_of_Chromatins);
            chromatin_set.resize(GenConst::Num_of_Chromatins);
            for (int i=0; i<GenConst::Num_of_Chromatins; i++) {
                string label=GenConst::Chromatin_label+to_string(i);
                Chromatins[i].set_label(label);
                Chromatins[i].set_file_time(buffer);
                Chromatins[i].set_index(i);
                Chromatins[i].import_config(chromatin_config_list[i]);
                chromatin_set[i].resize(Chromatins[i].get_num_of_node_types() );
                Chromatins[i].generate_report();
            }
        }
        
        if (GenConst::Num_of_pointparticles!=0){
            Include_pointparticle=true;
            pointparticles.resize(GenConst::Num_of_pointparticles);
            for (int i=0; i<GenConst::Num_of_pointparticles; i++) {
                pointparticles[i].set_file_time(buffer);
                pointparticles[i].set_index(i);
                pointparticles[i].import_config(pointparticle_config_list[i]);
                pointparticles[i].generate_report();
            }
            
        }
        
        //Used in the old code (no OpenMM engine)
        if (!GenConst::OpenMM) {
            if (Include_Membrane && Include_ECM) {
                for (int i=0; i<GenConst::Num_of_Membranes; i++) {
                    for (int j=0; j<GenConst::Num_of_ECMs; j++) {
                        Membrane_ECM_neighbour_finder(ECMs[j], Membranes[i]);
                    }
                }
            }
            if (Include_Membrane && Membranes.size()>1) {
                for (int i=1; i<GenConst::Num_of_Membranes; i++) {
                    Vesicle_particle_neighbour_finder (Membranes[i], Membranes[0]);
                }
            }
        }
        
        
        
        
        
        
        if (Include_Membrane) {
            for (int i=0; i<Membranes.size(); i++) {
                num_of_atoms        += Membranes[i].get_num_of_nodes();
                num_of_bonds        += Membranes[i].get_num_of_node_pairs();
                num_of_dihedrals    += Membranes[i].get_num_of_triangle_pairs();
            }
        }
        
        if (Include_Actin) {
            for (int i=0; i<Actins.size(); i++) {
                num_of_atoms        += Actins[i].get_num_of_nodes();
                num_of_bonds        += Actins[i].get_num_of_node_pairs();
            }
        }
        if (Include_ECM) {
            for (int i=0; i<ECMs.size(); i++) {
                num_of_atoms += ECMs[i].get_num_of_nodes();
                num_of_bonds += ECMs[i].get_num_of_node_pairs();
            }
        }
        if (Include_Chromatin) {
            for (int i=0; i<Chromatins.size(); i++) {
                num_of_atoms    += Chromatins[i].get_num_of_nodes();
                num_of_bonds    += Chromatins[i].get_num_of_nodes()-1;
            }
        }
        
        
        if (Include_pointparticle){
            num_of_atoms+=GenConst::Num_of_pointparticles;
        }
        
        if (Include_Membrane){
            if (Include_Actin){
                for (int i=0; i<GenConst::Num_of_Actins; i++) {
                    for (int j=0; j<GenConst::Num_of_Membranes; j++) {
                        Actin_Membrane_shared_Node_Identifier(Actins[i], Membranes[j] , i, j);
                        num_of_bonds        += Actins[i].return_num_of_actin_membrane_shared_nodes(j);
                        if (Membranes[j].get_relax_with_actin_flag()) {
                            Membranes[j].Relax_1();
                        }
                    }
                    
                } //for (int i=0; i<GenConst::Num_of_Actins; i++)
            } else { //if (Include_Actin){
                for (int i=0; i<Membranes.size(); i++){
                    Membranes[i].Relax_1();
                }// End of for (int i=0; i<Membranes.size(); i++)
            }//end else
        } // End of if (Include_Membrane)
    }
    
    float progressp=0;
    int   progress=0;
    //openmm**
    if (GenConst::OpenMM) {
        cout<<"\nBeginnig the OpenMM section:\n";
        std::string   platformName;
        int atom_count=0;
        int bond_count=0;
        int dihe_count=0;
        
        int mem_atom_count=0;
        //int act_atom_count=0;
        
        //The +1 is for the last member of the list that is set to -1 to indicate the end of list.
        MyAtomInfo* all_atoms    = new MyAtomInfo[num_of_atoms+1];
        Bonds* all_bonds         = new Bonds[num_of_bonds+1];
        Dihedrals* all_dihedrals = new Dihedrals[num_of_dihedrals+1];
        
        all_atoms[num_of_atoms].type         =EndOfList;
        all_bonds[num_of_bonds].type         =EndOfList;
        all_dihedrals[num_of_dihedrals].type =EndOfList;
        
        if (Include_Membrane) {
            OpenMM_membrane_info_relay(Membranes,
                                       membrane_set,
                                       all_atoms,
                                       all_bonds,
                                       all_dihedrals,
                                       atom_count,
                                       bond_count,
                                       dihe_count);
        }
        
        mem_atom_count = atom_count;
        
        
        if (Include_Actin) {
            OpenMM_Actin_info_relay(Actins,
                                    actin_set,
                                    all_atoms,
                                    all_bonds,
                                    all_dihedrals,
                                    atom_count,
                                    bond_count,
                                    dihe_count);
        }
        
        
        if (Include_Membrane  && Include_Actin) {
            OpenMM_ActMem_info_relay(Actins,
                                     Membranes,
                                    all_bonds,
                                    mem_atom_count,
                                    bond_count);

        }
        
        if (Include_ECM) {
            OpenMM_ECM_info_relay(ECMs,
                                  ecm_set,
                                  all_atoms,
                                  all_bonds,
                                  all_dihedrals,
                                  atom_count,
                                  bond_count,
                                  dihe_count);
        }
        if (Include_Chromatin) {
            OpenMM_Chromatin_info_relay(Chromatins,
                                        chromatin_set,
                                        all_atoms,
                                        all_bonds,
                                        all_dihedrals,
                                        atom_count,
                                        bond_count,
                                        dihe_count);
        }
       
        
        
        //autocorrelation calculations:
//        GenConst::velocity_save.resize(6);
        
        
        // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
        // usage and runtime errors are caught and reported.
        
        try {
            MyOpenMMData* omm = new MyOpenMMData();
            TimeDependantData* time_dependant_data = new TimeDependantData();
            if (!GenConst::Load_from_checkpoint) {
                omm = myInitializeOpenMM(all_atoms, GenConst::Step_Size_In_Fs, platformName, time_dependant_data, all_bonds, all_dihedrals, membrane_set, actin_set, ecm_set, chromatin_set, interaction_map);
            } else {
                std::filebuf rfb;
                string checkpoint_load_name = GenConst::Checkpoint_path + GenConst::Checkpoint_file_name;
                checkpoint_load_name = "Results/Resumes/OpenMM/chromo2019_11_17_time_11_57";
                rfb.open (checkpoint_load_name.c_str(),std::ios::in);
                if (rfb.is_open()) {
                    cout<<"Loading checkpoint from: "<<checkpoint_load_name<<endl;
                } else {
                    cout<<"Checkpoint loadding is enabled. Cannpt find the checkpoint file. Please check the checkpoint path and filename in the general config file."<<endl;
                    exit(EXIT_FAILURE);
                }
                std::istream rcheckpoint(&rfb);
                omm->context->loadCheckpoint(rcheckpoint);
                
                //wrok in progress.
                //Need to retrive all information from the checkpoint and relay them to the respective classes.
            }
            // Run the simulation:
            //  (1) Write the first line of the PDB file and the initial configuration.
            //  (2) Run silently entirely within OpenMM between reporting intervals.
            //  (3) Write a PDB frame when the time comes.
            printf("REMARK  Using OpenMM platform %s\n", platformName.c_str());
            std::filebuf wfb;
            wfb.open (ckeckpoint_name.c_str(),std::ios::out);
            std::ostream wcheckpoint(&wfb);
            
            //Time the programme
            tStart = clock();
            chrono_clock_start = chrono::steady_clock::now();
            chrono_sys_clock_start = chrono::system_clock::now();
            
            std::string traj_name="Results/"+GenConst::trajectory_file_name+buffer+".pdb";
            
            const int NumSilentSteps = (int)(GenConst::Report_Interval_In_Fs / GenConst::Step_Size_In_Fs + 0.5);
            
            int total_step_num = 0;
            double last_update_time=0;
            bool expanding = false;
            bool set_spring = false;
            
            for (int frame=1; ; ++frame) {

                double time, energyInKJ, potential_energyInKJ;
                
                myGetOpenMMState(omm, time, energyInKJ, potential_energyInKJ, all_atoms);
                if(GenConst::WriteVelocitiesandForces){
                    collect_data(all_atoms, buffer, Chromatins, Membranes, time);
                }
                myWritePDBFrame(frame, time, energyInKJ, potential_energyInKJ, all_atoms, all_bonds, traj_name);
                
                //Begin: Exporting congiguration of classes for simulation resume.
                Export_classes_for_resume(Membranes, Actins, ECMs, Chromatins, time, all_atoms);
                
                if (GenConst::CreateCheckpoint) {
                    omm->context->createCheckpoint(wcheckpoint);
                    //End: Exporting congiguration of classes for simulation resume.
                }
                
                
                
                if (time >= GenConst::Simulation_Time_In_Ps)
                    break;
                
                myStepWithOpenMM(omm,time_dependant_data, all_atoms, NumSilentSteps, total_step_num);
                
                if (100*time/GenConst::Simulation_Time_In_Ps>progressp){
//                    printf("[ %2.1f ]\t time: %4.1f Ps [out of %4.1f Ps]    \r",progressp, time, GenConst::Simulation_Time_In_Ps);
                    cout<<"[ "<<int(progressp*10)/10.0<<"% ]   \t time: "<<time<<" Ps [out of "<<GenConst::Simulation_Time_In_Ps<<" Ps]    \r" << std::flush;
                    progressp+=0.1;
                    progress++;
                }
                
                if (check_for_membrane_update(Membranes, time, last_update_time)) {
                    updateOpenMMforces(Membranes, Chromatins, omm, time, all_atoms, all_bonds, membrane_set, interaction_map);
                }
                
                if (time>100 && set_spring) {
                    int atom1, atom2 ;
                    double length, stiffness;
                    
                    for(int i=0; i<Membranes[0].get_num_of_node_pairs() ; i++)
                    {
                        omm->harmonic->getBondParameters(i, atom1, atom2, length, stiffness);
                        stiffness = 0.5 * OpenMM::KJPerKcal
                                        * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
                        
                        omm->harmonic->setBondParameters(i, atom1, atom2, length, stiffness);
                    }
                    omm->harmonic->updateParametersInContext(*omm->context);
                    set_spring=false;
                }
                
                if (time>800 && expanding) {
                    expand(Chromatins, omm);
                    expanding=false;
                }
                //the monte_carlo part

                 //if(progress==0 or progress==25 or progress==50 or progress==75){ 
                if ((progress%5==0 or progress==0) and GenConst::MC_step !=0){
                 
                    //Membranes[0].check_the_flip(omm, all_bonds , all_dihedrals);
                    

                    Monte_Carlo_Reinitialize(omm, all_bonds , all_dihedrals, Membranes[0], all_atoms);
                }
                
             
            }
            
            cout<<"[ 100% ]\t time: "<<GenConst::Simulation_Time_In_Ps<<"Ps\n";
            
            
            print_wall_clock_time((double)((clock() - tStart)/CLOCKS_PER_SEC));
            print_real_time(chrono_clock_start, chrono::steady_clock::now());
            print_system_time(chrono_sys_clock_start, chrono::system_clock::now());
            
            
            
            
            
            // Clean up OpenMM data structures.
            myTerminateOpenMM(omm,time_dependant_data);
            
            cout<<"\nDone!"<<endl;
            return 0; // Normal return from main.
        }
        
        // Catch and report usage and runtime errors detected by OpenMM and fail.
        catch(const std::exception& e) {
            printf("EXCEPTION: %s\n", e.what());
            return 1;
        }
    }
    
    
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    progress=0;
    
    cout<<"\nBeginnig the MD\nProgress:\n";
    for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_steps ; MD_Step++){
        //        cout<<Membranes[0].return_node_position(0, 0);
        
        //Thermostat step first step
        if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {
            if (Include_Membrane) {
                for (int i=0; i<Membranes.size(); i++) {
                    Membranes[i].Thermostat_Bussi(GenConst::Buffer_temperature);
                }
            }
            if (Include_Actin) {
                for (int i=0; i<Actins.size(); i++) {
                    //                    Actins[i].Thermostat_Bussi(GenConst::MD_T);
                }
            }
            if (Include_Chromatin) {
                for (int i=0; i<Chromatins.size(); i++) {
                    Chromatins[i].Thermostat_Bussi(GenConst::MD_T*0.01);
                }
            }
            if (Include_ECM) {
                for (int i=0; i<ECMs.size(); i++) {
                    //                    ECMs[i].Thermostat_Bussi(GenConst::MD_T);
                }
            }
        }
        
        
        //Velocity Verlet first step
        if (Include_Membrane)
        {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].MD_Evolution_beginning(GenConst::Step_Size_In_Fs);
            }
        }
        if (Include_Chromatin)
        {
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatins[i].MD_Evolution_beginning(GenConst::Step_Size_In_Fs);
            }
        }
        if (Include_Actin)
        {
            for (int i=0; i<Actins.size(); i++) {
                Actins[i].MD_Evolution_beginning(GenConst::Step_Size_In_Fs);
            }
        }
        if (Include_ECM)
        {
            for (int i=0; i<ECMs.size(); i++) {
                ECMs[i].MD_Evolution_beginning(GenConst::Step_Size_In_Fs);
            }
        }
        
        if (Include_pointparticle)
        {
            for (int i=0; i<pointparticles.size(); i++) {
                if (pointparticles[i].on_or_off_MD_evolution){
                    pointparticles[i].MD_Evolution_beginning(GenConst::Step_Size_In_Fs);
                }
            }
        }
        
        //force implamentation
        if (Include_Membrane)
        {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].Elastic_Force_Calculator(0);
            }
        }
        if (Include_Chromatin)
        {
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatins[i].Force_Calculator_2();
            }
        }
        if (Include_Actin)
        {
            for (int i=0; i<Actins.size(); i++) {
                Actins[i].Elastic_Force_Calculator();
            }
        }
        if (Include_ECM)
        {
            for (int i=0; i<ECMs.size(); i++) {
                //                ECMs[i].Elastic_Force_Calculator();
            }
        }
        
        //Shared Forces
        if (Include_Chromatin && Include_Membrane) {
            if (MD_Step%2000==0) {
                for (int i=0; i<Chromatins.size(); i++) {
                    Chromatin_Membrane_neighbour_finder(Chromatins[i], Membranes[i]);
                    Chromatin_Membrane_hard_sphere(Chromatins[i], Membranes[i]);
                }
            }
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatin_Membrane_hard_sphere(Chromatins[i], Membranes[i]);
            }
            
            
        }
        
        if (Include_Membrane && Include_Actin) {
            for (int i=0; i<Actins.size(); i++) {
                for (int j=0; j<Membranes.size(); j++) {
                    Actin_Membrane_shared_Node_Force_calculator(Actins[i], Membranes[j], j);
                }
            }
        }
        
        
        if (Include_Membrane && Membranes.size()>1) {
            for (int i=1; i<Membranes.size(); i++) {
                particle_vesicle_shared_node_force (Membranes[i] , Membranes[0]);
                if (MD_Step%2000==0) {
                    update_particle_vesicle_neighbour_list (Membranes[i] , Membranes[0]);
                }
            }
        }
        
        
        if (Include_Membrane && Include_ECM) {
            for (int i=0; i<Membranes.size(); i++) {
                for (int j=0; j<ECMs.size(); j++) {
                    Membrane_ECM_shared_node_force (ECMs[j], Membranes[i]);
                    if (MD_Step%2000==0) {
                        update_ecm_mem_neighbour_list (ECMs[j], Membranes[i]);
                    }
                }
            }
        }
        
        
        if (Include_Membrane && Include_pointparticle) {
            for (int i=0; i<Membranes.size(); i++) {
                for (int j=0; j<pointparticles.size(); j++){
                    Vesicle_pointparticle_neighbour_finder (pointparticles[j], Membranes[i]);
                    pointparticle_vesicle_shared_node_force (pointparticles[j], Membranes[i]);
                }
            }
        }
        
        if(Include_pointparticle && pointparticles.size()>1){
            for (int i=0; i<pointparticles.size(); i++){
                for (int j=0; j<pointparticles.size(); j++){
                    if(i != j){
                        pointparticle_pointparticle_interaction(pointparticles[i],pointparticles[j]);
                    }
                }
            }
        }
        
        //Velocity Verlet second step
        if (Include_Membrane) {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].MD_Evolution_end(GenConst::Step_Size_In_Fs);
            }
        }
        if (Include_Chromatin) {
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatins[i].MD_Evolution_end(GenConst::Step_Size_In_Fs);
            }
        }
        if (Include_Actin) {
            for (int i=0; i<Actins.size(); i++) {
                Actins[i].MD_Evolution_end(GenConst::Step_Size_In_Fs);
            }
        }
        if (Include_ECM) {
            for (int i=0; i<ECMs.size(); i++) {
                ECMs[i].MD_Evolution_end(GenConst::Step_Size_In_Fs);
            }
        }
        
        if (Include_pointparticle)
        {
            for (int i=0; i<pointparticles.size(); i++) {
                if (pointparticles[i].on_or_off_MD_evolution){
                    pointparticles[i].MD_Evolution_end(GenConst::Step_Size_In_Fs);
                    
                }
            }
        }
        
        
        //Thermostat second step
        if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {
            if (Include_Membrane) {
                for (int i=0; i<Membranes.size(); i++) {
                    Membranes[i].Thermostat_Bussi(GenConst::Buffer_temperature);
                }
            }
            if (Include_Actin) {
                for (int i=0; i<Actins.size(); i++) {
                    //                    Actins[i].Thermostat_Bussi(GenConst::MD_T);
                }
            }
            if (Include_Chromatin) {
                for (int i=0; i<Chromatins.size(); i++) {
                    Chromatins[i].Thermostat_Bussi(GenConst::MD_T*0.01);
                }
            }
            if (Include_ECM) {
                for (int i=0; i<ECMs.size(); i++) {
                    //                    ECMs[i].Thermostat_Bussi(GenConst::MD_T);
                }
            }
        }
        
        //saving Results
        if (MD_Step%GenConst::MD_traj_save_step == 0)
        {
            Trajectory << num_of_atoms<<endl;
            Trajectory << " nodes  "<<endl;
            
            
            if (Include_Membrane) {
                for (int i=0; i<Membranes.size(); i++) {
                    //                    string label="Membrane_"+to_string(i);
                    Membranes[i].write_traj(traj_file_name);
                    Membranes[i].export_for_resume(MD_Step);
                }
            }
            
            if (Include_Chromatin) {
                for (int i=0; i<Chromatins.size(); i++) {
                    string label="Chromatin_"+to_string(i);
                    Chromatins[i].write_traj(traj_file_name, label);
                    Chromatins[i].export_for_resume(MD_Step);
                }
            }
            if (Include_Actin) {
                for (int i=0; i<Actins.size(); i++) {
                    string label="Actin_"+to_string(i);
                    Actins[i].write_traj(traj_file_name, label);
                    //                    Actins[i].export_for_resume(MD_Step);
                }
            }
            if (Include_ECM) {
                for (int i=0; i<ECMs.size(); i++) {
                    string label="ECM_"+to_string(i);
                    ECMs[i].write_traj(traj_file_name, label);
                    //                    Actins[i].export_for_resume(MD_Step);
                }
            }
            
            if (Include_pointparticle) {
                for (int i=0; i<pointparticles.size(); i++) {
                    
                    pointparticles[i].write_traj(traj_file_name);
                    
                }
            }
        }// End of if (MD_Step%100==0)
        
        
        if (int(100*MD_Step/GenConst::MD_num_of_steps)>progress){
            cout<<"[ "<<progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;
            progress+=5;
        }
        
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    cout<<"[ 100% ]\t step: "<<GenConst::MD_num_of_steps<<"\n";
    cout<<"\nDone!"<<endl;
    printf("Time taken: %.2f Minutes\n", (double)((clock() - tStart)/CLOCKS_PER_SEC)/60.0);
    return 0;
}


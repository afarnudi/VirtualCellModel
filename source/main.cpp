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
#include "Global_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_class_functions.h"
#include "Arg_pars.hpp"
#include "Configfile.hpp"
#include "Interaction_table.hpp"

//#include "Tests.hpp"

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
PotentialModelIndex potential;
double Simulation_Time_In_Ps;
double Report_Interval_In_Fs;
double Step_Size_In_Fs;
int    MC_step;
int    Mem_fluidity;
bool   Periodic_box;
double Lbox;
double Simulation_box_length;
bool   Periodic_condtion_status;
int    Num_of_Membranes;
int    Num_of_Chromatins;
int    Num_of_Actins;
int    Num_of_ECMs;
string trajectory_file_name;
string ProjectName;
string force_file_name;
double Buffer_temperature; //***********OLDCODE
double Bussi_tau;
double Actin_Membrane_Bond_Coefficient;
bool   Excluded_volume_interaction;
double sigma_LJ_12_6;
double epsilon_LJ_12_6;
string Membrane_label;
string Actin_label;
string Chromatin_label;
string ECM_label;
int    Integrator_type;
double frictionInPs;
double temperature;
bool   CreateCheckpoint;
bool   Load_from_checkpoint;
string Checkpoint_path;
string Checkpoint_file_name;
bool   ChromatinVirtualSites;

bool   WantEnergy;
bool   WantForce;
bool   WantVelocity;
bool   CMMotionRemover;
int    CMMotionRemoverStep;
bool   Wantvoronoi;
bool   Testmode;


double MCBarostatPressure;
double MCBarostatTemperature;
int    MCBarostatFrequency;

std::vector<double> PeriodicBoxVector0;
std::vector<double> PeriodicBoxVector1;
std::vector<double> PeriodicBoxVector2;


//    std::vector<std::vector<std::vector<double> > > data;
std::vector<double> data_colection_times;
std::vector<std::vector<double> > Lboxdims;

std::string hardwareReport;
}



const int EndOfList=-1;


int main(int argc, char **argv)
{
    cout<<TRESET;
    // get the current time.
    time_t t = time(0);
    auto chrono_clock_start = chrono::steady_clock::now();
    auto chrono_sys_clock_start = chrono::system_clock::now();
    
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_time_%H_%M",now);
    
    string configfilename;
    string general_file_name="General_param_map.txt";
    
    configfilename = cxxparser_vcm(argc, argv);
    
    if (configfilename == "None") {
        cout<<TBOLD<<"\nHi!\nPlease enter the path + name of the configuration file. If you do not have a configuration file, run \""<<argv[0]<<" -h\" for more options.\n"<<TRESET<<"Example:\t../../myconfigfile.txt\n\nPath to configuration file: ";
        cout<<TFILE;
        cin>>configfilename;
        cout<<TRESET;

    }
    
    
    clock_t tStart = clock();//Time the programme

    map<string, vector<string> > config_lines =read_configfile(configfilename);
    get_class_numbers(config_lines);
    parse_genconfig_parameters(config_lines["-GeneralParameters"]);

    NonBondInteractionMap interaction_map(config_lines["-InteractionTable"]);
    
    string ckeckpoint_name=GenConst::trajectory_file_name+"_Checkpoint";
    
    vector<Membrane> Membranes;
    vector<std::set<int> >  membrane_set;
    vector<vector<string> > membrane_configs =sort_class_configs(config_lines["-Membrane"]);
    
    vector<Actin> Actins;
    vector<std::set<int> > actin_set;
    vector<vector<string> > actin_configs =sort_class_configs(config_lines["-Actin"]);

    vector<ECM> ECMs;
    vector<std::set<int> > ecm_set;
    vector<vector<string> > ecm_configs =sort_class_configs(config_lines["-ECM"]);
    
    vector<Chromatin> Chromatins;
    vector<vector<std::set<int> > > chromatin_set;
    vector<vector<string> > chromatin_configs =sort_class_configs(config_lines["-Chromatin"]);
    
    int num_of_atoms=0;
    int num_of_bonds=0;
    int num_of_dihedrals=0;
    
    if (!GenConst::Load_from_checkpoint) {
        if (GenConst::Num_of_Membranes!=0) {
            Membranes.resize(GenConst::Num_of_Membranes);
            membrane_set.resize(GenConst::Num_of_Membranes);
            for (int i=0; i<GenConst::Num_of_Membranes; i++) {
                string label=GenConst::Membrane_label+to_string(i);
                Membranes[i].set_label(label);
                Membranes[i].set_file_time(buffer);
                Membranes[i].set_index(i);
                try{
                    Membranes[i].import_config(membrane_configs[i]);
                }
                catch(const std::exception& e) {
                    printf("EXCEPTION: %s\n", e.what());
                    return 0;
                }
            }
        }
        
        if (GenConst::Num_of_Actins!=0) {
            Actins.resize(GenConst::Num_of_Actins);
            actin_set.resize(GenConst::Num_of_Actins);
            for (int i=0; i<GenConst::Num_of_Actins; i++) {
                string label=GenConst::Actin_label+to_string(i);
                Actins[i].set_label(label);
                Actins[i].set_file_time(buffer);
                Actins[i].set_index(i);
                try {
                    Actins[i].import_config(actin_configs[i]);
                }
                catch(const std::exception& e) {
                    printf("EXCEPTION: %s\n", e.what());
                    return 0;
                }
            }
        }
        
        if (GenConst::Num_of_ECMs!=0){
            ECMs.resize(GenConst::Num_of_ECMs);
            ecm_set.resize(GenConst::Num_of_ECMs);
            for (int i=0; i<GenConst::Num_of_ECMs; i++) {
                string label=GenConst::ECM_label+to_string(i);
                ECMs[i].set_label(label);
                ECMs[i].set_file_time(buffer);
                ECMs[i].set_index(i);
                try {
                    ECMs[i].import_config(ecm_configs[i]);
                }
                catch(const std::exception& e) {
                    printf("EXCEPTION: %s\n", e.what());
                    return 0;
                }
            }
            
        }
        
        
        if (GenConst::Num_of_Chromatins!=0) {
            Chromatins.resize(GenConst::Num_of_Chromatins);
            chromatin_set.resize(GenConst::Num_of_Chromatins);
            for (int i=0; i<GenConst::Num_of_Chromatins; i++) {
                string label=GenConst::Chromatin_label+to_string(i);
                Chromatins[i].set_label(label);
                Chromatins[i].set_file_time(buffer);
                Chromatins[i].set_index(i);
                try {
                Chromatins[i].import_config(chromatin_configs[i]);
                }
                catch(const std::exception& e) {
                    printf("EXCEPTION: %s\n", e.what());
                    return 0;
                }
                chromatin_set[i].resize(Chromatins[i].get_num_of_node_types() );
            }
        }
        
        
        if (GenConst::Num_of_Membranes!=0) {
            
            for (int i=0; i<Membranes.size(); i++) {
                num_of_atoms        += Membranes[i].get_num_of_nodes();
                num_of_bonds        += Membranes[i].get_num_of_node_pairs();
                num_of_dihedrals    += Membranes[i].get_num_of_triangle_pairs();
            }
        }
        
        if (GenConst::Num_of_Actins!=0) {
            for (int i=0; i<Actins.size(); i++) {
                num_of_atoms        += Actins[i].get_num_of_nodes();
                num_of_bonds        += 4*Actins[i].get_num_of_node_pairs() + 4*Actins[i].get_num_of_abp_pairs() + 4*Actins[i].get_num_of_MT_pairs();
            }
        }
        if (GenConst::Num_of_ECMs!=0) {
            for (int i=0; i<ECMs.size(); i++) {
                num_of_atoms += ECMs[i].get_num_of_nodes();
                num_of_bonds += ECMs[i].get_num_of_node_pairs();
            }
        }
        if (GenConst::Num_of_Chromatins!=0) {
            for (int i=0; i<Chromatins.size(); i++) {
                num_of_atoms    += Chromatins[i].get_num_of_nodes();
                num_of_bonds    += Chromatins[i].get_num_of_bonds();
            }
        }
        
        
        if (GenConst::Num_of_Membranes!=0){
            if (GenConst::Num_of_Actins!=0){
                for (int i=0; i<GenConst::Num_of_Actins; i++) {
                    for (int j=0; j<GenConst::Num_of_Membranes; j++) {
                        Actin_Membrane_shared_Node_Identifier(Actins[i],Membranes[j],i,j);
                        num_of_bonds        += Actins[i].return_num_of_actin_membrane_shared_nodes(j);
                    }
                    
                } //for (int i=0; i<GenConst::Num_of_Actins; i++)
            }
        } // End of if (Include_Membrane)
    }
    
    
    float progressp=0;
    
    int progress=0;
    double MC_Acceptance_Rate=0;
    int MC_total_tries=0;
    int Accepted_Try_Counter=0;
    
    cout<<TCYAN<<"\nBeginnig the OpenMM section:\n"<<TRESET;
    std::string   platformName;
    int atom_count=0;
    int bond_count=0;
    int dihe_count=0;
    
    int mem_atom_count=0;
    //int act_atom_count=0;
    
    //The +1 is for the last member of the list that is set to -1 to indicate the end of list.
    MyAtomInfo* all_atoms     = new MyAtomInfo[num_of_atoms+1];
    Bonds*      all_bonds     = new Bonds[num_of_bonds+1];
    Dihedrals*  all_dihedrals = new Dihedrals[num_of_dihedrals+1];
    
    cout<<num_of_atoms<<" Atoms"<<endl;
    cout<<num_of_bonds<<" Bonds"<<endl;
    
    all_atoms[num_of_atoms].type         =EndOfList;
    all_bonds[num_of_bonds].type         =EndOfList;
    all_dihedrals[num_of_dihedrals].type =EndOfList;
    
    
    if (GenConst::Num_of_Membranes!=0) {
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
    
    
    if (GenConst::Num_of_Actins!=0) {
        OpenMM_Actin_info_relay(Actins,
                                actin_set,
                                all_atoms,
                                all_bonds,
                                all_dihedrals,
                                atom_count,
                                bond_count,
                                dihe_count);
    }
    
    
    if (GenConst::Num_of_Membranes!=0  && GenConst::Num_of_Actins!=0) {
        OpenMM_ActMem_info_relay(Actins,
                                 Membranes,
                                 all_bonds,
                                 mem_atom_count,
                                 bond_count);
        
    }
    
    if (GenConst::Num_of_ECMs!=0) {
        OpenMM_ECM_info_relay(ECMs,
                              ecm_set,
                              all_atoms,
                              all_bonds,
                              all_dihedrals,
                              atom_count,
                              bond_count,
                              dihe_count);
    }
    if (GenConst::Num_of_Chromatins!=0) {
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
        assign_project_directories(buffer);
        cout<< "\nFile name: "<<TFILE<<GenConst::trajectory_file_name<<TRESET<<endl<<endl;
        
        ofstream write_hardwareReport;
        string hardwarereportpath =GenConst::trajectory_file_name+"_hardware_runtime.txt";
        write_hardwareReport.open(hardwarereportpath.c_str(),std::ios_base::app);
        write_hardwareReport<<GenConst::hardwareReport;
        write_hardwareReport<<endl;
        
        //export generated chromatin coordinates; The export happens only if the chromatin is set to export coordinates in it's configuration.
        for (int ich=0;ich<Chromatins.size();ich++) {
            Chromatins[ich].export_coordinates();
        }
        // Run the simulation:
        //  (1) Write the first line of the PDB file and the initial configuration.
        //  (2) Run silently entirely within OpenMM between reporting intervals.
        //  (3) Write a PDB frame when the time comes.
        
        cout<<"REMARK  Using OpenMM platform "<<TBOLD;
        if (platformName=="OpenCL") {
            cout<<TOCL;
        } else if (platformName=="CUDA"){
            cout<<TCUD;
        } else {
            cout<<TCPU;
        }
        cout<<platformName.c_str()<<TRESET<<endl;
        
        
        std::string Reporname=GenConst::trajectory_file_name+"_report.txt";
        
        generate_report(config_lines);
        
        
        std::filebuf wfb;
        wfb.open (ckeckpoint_name.c_str(),std::ios::out);
        std::ostream wcheckpoint(&wfb);
        
        //Time the programme
        tStart = clock();
        chrono_clock_start = chrono::steady_clock::now();
        chrono_sys_clock_start = chrono::system_clock::now();
        
        
        
        //int SavingStep     = (int)(GenConst::Report_Interval_In_Fs / GenConst::Step_Size_In_Fs + 0.5);
        int MCCalcstep = (int)(GenConst::Mem_fluidity * GenConst::Step_Size_In_Fs + 0.5);
        int NumSilentSteps = (int)(GenConst::Report_Interval_In_Fs / GenConst::Step_Size_In_Fs + 0.5);
        int Savingstep = NumSilentSteps;
        if ( (MCCalcstep < NumSilentSteps) && MCCalcstep != 0 ) {
            Savingstep = (int)(NumSilentSteps / MCCalcstep )* MCCalcstep;
            NumSilentSteps = MCCalcstep;
        } else if ( (NumSilentSteps < MCCalcstep) && MCCalcstep != 0 ){
            int rate =  (int)(MCCalcstep / NumSilentSteps );
            Savingstep = int(MCCalcstep/rate);
            NumSilentSteps = Savingstep;
        }
        
        int savetime = 0;
        int MCCalcTime = MCCalcstep;
        
        cout<<"Savingstep "<<Savingstep<<"\nNumSilentSteps "<<NumSilentSteps<<endl;
        
        int total_step_num = 0;
        
        
        bool expanding = false;
        bool set_spring = false;
        bool changeTemp = false;
        double TempStep = 0.2;
        double initTemp = GenConst::temperature;
        
        myWritePDBFrame(0, 0, 0, 0, all_atoms, all_bonds);
        myWritePSF(num_of_atoms, num_of_bonds, all_atoms, all_bonds);
        
        for (int frame=1; ; ++frame) {
            
            double time, energyInKJ, potential_energyInKJ;
            
            myGetOpenMMState(omm, time, energyInKJ, potential_energyInKJ, all_atoms);
            
            if ( int(time*1000/GenConst::Step_Size_In_Fs) > savetime ) {
                Update_classes(Membranes, Actins, ECMs, Chromatins, time, all_atoms);
                
                collect_data(omm, all_atoms, interaction_map, Membranes, time);
                myWritePDBFrame(frame, time, energyInKJ, potential_energyInKJ, all_atoms, all_bonds);
                
                //Begin: Exporting congiguration of classes for simulation .
                
                
                savetime += Savingstep;
                
            }
            
            
            if (changeTemp) {
                initTemp -= TempStep;
                omm->Lintegrator->setTemperature(initTemp);
                //                    cout<<"temp: "<<omm->Lintegrator->getTemperature()<<endl;
                //                    map <string, double> params;
                //                    params = omm->context->getParameters();
                //
                //                    cout<<"\nframe: "<<frame<<endl;
                //                    cout<<params.size()<<endl;
                //                    for(auto elem : params)
                //                    {
                //                       cout << elem.first << " " << elem.second << "\n";
                //                    }
                //                    cout<<"\n";
                //
            }
            //                cout<<"CreateCheckpoint\n";
            
            if (GenConst::CreateCheckpoint) {
                omm->context->createCheckpoint(wcheckpoint);
                //End: Exporting congiguration of classes for simulation resume.
            }
            
            if (time >= GenConst::Simulation_Time_In_Ps){
                break;
            }
            //                cout<<"myStepWithOpenMM\n";
            
            myStepWithOpenMM(omm,time_dependant_data, all_atoms, NumSilentSteps, total_step_num);
            
            if (100*time/GenConst::Simulation_Time_In_Ps>progressp){
                printf("[ %2.1f ] time: %4.1f Ps [out of %4.1f Ps]    \r",100*time/GenConst::Simulation_Time_In_Ps, time, GenConst::Simulation_Time_In_Ps);
                cout<< std::flush;
                //                    cout<<"[ "<<int(progressp*10)/10.0<<"% ]   \t time: "<<time<<" Ps [out of "<<GenConst::Simulation_Time_In_Ps<<" Ps]    \r" << std::flush;
                progressp =  int(1000*time/GenConst::Simulation_Time_In_Ps)/10. + 0.1;
                //                    cout<<progressp<<endl;
                progress++;
            }
            
            
            if (check_for_membrane_update(Membranes, time)) {
                updateOpenMMforces(Membranes, Chromatins, omm, time, all_atoms, all_bonds, membrane_set, interaction_map);
            }
            
            
            //                if (time>100 && set_spring) {
            //                    int atom1, atom2 ;
            //                    double length, stiffness;
            //
            //                    for(int i=0; i<Membranes[0].get_num_of_node_pairs() ; i++)
            //                    {
            //                        omm->harmonic->getBondParameters(i, atom1, atom2, length, stiffness);
            //                        stiffness = 0.5 ;
            //
            //                        omm->harmonic->setBondParameters(i, atom1, atom2, length, stiffness);
            //                    }
            //                    omm->harmonic->updateParametersInContext(*omm->context);
            //                    set_spring=false;
            //                }
            //
            //                if (time>800 && expanding) {
            //                    expand(Chromatins, omm);
            //                    expanding=false;
            //                }
            
            
            
            //the monte_carlo part
            
            //if(progress==0 or progress==25 or progress==50 or progress==75){
            
            if ( int(time*1000/GenConst::Step_Size_In_Fs) >= MCCalcTime and GenConst::Mem_fluidity !=0){
                
                //Membranes[0].check_the_flip(omm, all_bonds , all_dihedrals);
                //                    cout<<"Frame  "<<frame<<endl;
                
                Monte_Carlo_Reinitialize(omm, all_bonds , all_dihedrals, Membranes[0], all_atoms, MC_total_tries,Accepted_Try_Counter, MC_Acceptance_Rate);
                
                MCCalcTime += MCCalcstep;
            }
            if(frame%50==2 and  GenConst::Mem_fluidity !=0){
                
                cout<<"\n total monte_carlo tries  "<<MC_total_tries<<endl;
                cout<<"total accepted tries"<<Accepted_Try_Counter<<endl;
                cout<<"acceptance_rate  "<<MC_Acceptance_Rate<<endl;
            }

        }
        
        cout<<"[ 100% ]\t time: "<<GenConst::Simulation_Time_In_Ps<<"Ps\n";
        
        
        print_wall_clock_time((double)((clock() - tStart)/CLOCKS_PER_SEC));
        print_real_time(chrono_clock_start, chrono::steady_clock::now());
        print_system_time(chrono_sys_clock_start, chrono::system_clock::now());
        
        
        
        
        
        // Clean up OpenMM data structures.
        myTerminateOpenMM(omm,time_dependant_data);
        cout<<"MC_Acceptance_Rate   "<<MC_Acceptance_Rate<<endl;
        cout<<"\nDone!!!!!"<<endl;
        return 0; // Normal return from main.
    }
    
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 0;
    }
}


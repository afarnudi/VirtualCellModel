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
#include <cstdio>
#include <boost/filesystem.hpp>

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

GeneralParameters generalParameters;
/**A struct that stores differnet potential indexies.*/
PotentialModelIndex potentialModelIndex;

namespace GenConst {

int    MC_step;
int    Mem_fluidity;

string force_file_name;
double Buffer_temperature; //***********OLDCODE
double Bussi_tau;
double Actin_Membrane_Bond_Coefficient;

double sigma_LJ_12_6;
double epsilon_LJ_12_6;


//bool   CreateCheckpoint;
bool   Load_from_checkpoint;
string Checkpoint_path;
string Checkpoint_file_name;
bool   ChromatinVirtualSites;


std::vector<double> data_colection_times;


}



const int EndOfList=-1;


int main(int argc, char **argv)
{
//    generalParameters.CBP=true;
    cout<<TRESET;
    // get the current time.
    time_t t = time(0);
    auto chrono_clock_start = chrono::steady_clock::now();
    auto chrono_sys_clock_start = chrono::system_clock::now();
    
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_time_%H_%M",now);
    
    ArgStruct_VCM userinputs = cxxparser_vcm(argc, argv);
    
    if (userinputs.configfilename == "None") {
        cout<<TBOLD<<"\nHi!\nPlease enter the path + name of the configuration file. If you do not have a configuration file, run \""<<argv[0]<<" -h\" for more options.\n"<<TRESET<<"Example:\t../../myconfigfile.txt\n\nPath to configuration file: ";
        string configfilename;
        cout<<TFILE;
        cin>>configfilename;
        cout<<TRESET;
        userinputs.configfilename=check_if_file_exists(configfilename);
    } else {
        userinputs.configfilename=check_if_file_exists(userinputs.configfilename);
    }
    
    clock_t tStart = clock();//Time the programme

    map<string, vector<string> > config_lines =read_configfile(userinputs.configfilename);
    get_class_numbers(config_lines);
    parse_genconfig_parameters(config_lines["-GeneralParameters"]);

    NonBondInteractionMap interaction_map(config_lines["-InteractionTable"]);
    
    
    
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
    int num_of_angles=0;
    
//    if (!generalParameters.Resume) {
        if (generalParameters.Num_of_Membranes!=0) {
            Membranes.resize(generalParameters.Num_of_Membranes);
            membrane_set.resize(generalParameters.Num_of_Membranes);
            for (int i=0; i<generalParameters.Num_of_Membranes; i++) {
                string label=generalParameters.Membrane_label+to_string(i);
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
        
        if (generalParameters.Num_of_Actins!=0) {
            Actins.resize(generalParameters.Num_of_Actins);
            actin_set.resize(generalParameters.Num_of_Actins);
            for (int i=0; i<generalParameters.Num_of_Actins; i++) {
                string label=generalParameters.Actin_label+to_string(i);
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
        
        if (generalParameters.Num_of_ECMs!=0){
            ECMs.resize(generalParameters.Num_of_ECMs);
            ecm_set.resize(generalParameters.Num_of_ECMs);
            for (int i=0; i<generalParameters.Num_of_ECMs; i++) {
                string label=generalParameters.ECM_label+to_string(i);
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
        
        
        if (generalParameters.Num_of_Chromatins!=0) {
            Chromatins.resize(generalParameters.Num_of_Chromatins);
            chromatin_set.resize(generalParameters.Num_of_Chromatins);
            for (int i=0; i<generalParameters.Num_of_Chromatins; i++) {
                string label=generalParameters.Chromatin_label+to_string(i);
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
        
        
        if (generalParameters.Num_of_Membranes!=0) {
            
            for (int i=0; i<Membranes.size(); i++) {
                num_of_atoms        += Membranes[i].get_num_of_nodes();
                num_of_bonds        += Membranes[i].get_num_of_node_pairs();
                num_of_dihedrals    += Membranes[i].get_num_of_triangle_pairs();
            }
        }
        
        if (generalParameters.Num_of_Actins!=0) {
            for (int i=0; i<Actins.size(); i++) {
                num_of_atoms        += Actins[i].get_num_of_nodes();
//                num_of_bonds        += 4*Actins[i].get_num_of_node_pairs() + 4*Actins[i].get_num_of_abp_pairs() + 4*Actins[i].get_num_of_MT_pairs();
                num_of_bonds        += Actins[i].get_num_of_node_pairs();
            }
        }
        if (generalParameters.Num_of_ECMs!=0) {
            for (int i=0; i<ECMs.size(); i++) {
                num_of_atoms += ECMs[i].get_num_of_nodes();
                num_of_bonds += ECMs[i].get_num_of_node_pairs();
            }
        }
        if (generalParameters.Num_of_Chromatins!=0) {
            for (int i=0; i<Chromatins.size(); i++) {
                num_of_atoms    += Chromatins[i].get_num_of_nodes();
                num_of_bonds    += Chromatins[i].get_num_of_bonds();
                num_of_angles   += Chromatins[i].get_num_of_angle_bonds();
            }
        }
        
        
        if (generalParameters.Num_of_Membranes!=0){
            if (generalParameters.Num_of_Actins!=0){
                for (int i=0; i<generalParameters.Num_of_Actins; i++) {
                    for (int j=0; j<generalParameters.Num_of_Membranes; j++) {
                        Actin_Membrane_shared_Node_Identifier(Actins[i],Membranes[j],i,j);
                        num_of_bonds        += Actins[i].return_num_of_actin_membrane_shared_nodes(j);
                    }
                    
                } //for (int i=0; i<GenConst::Num_of_Actins; i++)
            }
        } // End of if (Include_Membrane)
//    } End of Resume
    
//    cout<<"num_of_bonds "<<num_of_bonds<<endl;
    float progressp=0;
    
//    int progress=0;
    double MC_Acceptance_Rate=0;
    int MC_total_tries=0;
    int Accepted_Try_Counter=0;
    
    cout<<TCYAN<<"\nBeginning the OpenMM section:\n"<<TRESET;
    std::string   platformName;
    int atom_count=0;
    int bond_count=0;
    int dihe_count=0;
    int angle_count=0;
    
    int mem_atom_count=0;
    //int act_atom_count=0;
    
    //The +1 is for the last member of the list that is set to -1 to indicate the end of list.
    MyAtomInfo* all_atoms     = new MyAtomInfo[num_of_atoms+1];
    Bonds*      all_bonds     = new Bonds[num_of_bonds+1];
    Dihedrals*  all_dihedrals = new Dihedrals[num_of_dihedrals+1];
    Angles*     all_angles    = new Angles[num_of_angles+1];
    
    
//    if (!generalParameters.Resume) {
        all_atoms[num_of_atoms].type         =EndOfList;
        all_bonds[num_of_bonds].type         =EndOfList;
        all_dihedrals[num_of_dihedrals].type =EndOfList;
        all_angles[num_of_angles].type       =EndOfList;
        
        if (generalParameters.Num_of_Membranes!=0) {
            OpenMM_membrane_info_relay(Membranes,
                                       membrane_set,
                                       all_atoms,
                                       all_bonds,
                                       all_dihedrals,
                                       atom_count,
                                       bond_count,
                                       dihe_count,
                                       interaction_map);
        }
        
        mem_atom_count = atom_count;
        
        
        if (generalParameters.Num_of_Actins!=0) {
            OpenMM_Actin_info_relay(Actins,
                                    actin_set,
                                    all_atoms,
                                    all_bonds,
                                    all_dihedrals,
                                    atom_count,
                                    bond_count,
                                    dihe_count);
        }
        
        
        if (generalParameters.Num_of_Membranes!=0  && generalParameters.Num_of_Actins!=0) {
            OpenMM_ActMem_info_relay(Actins,
                                     Membranes,
                                     all_bonds,
                                     mem_atom_count,
                                     bond_count);
            
        }
        
        if (generalParameters.Num_of_ECMs!=0) {
            OpenMM_ECM_info_relay(ECMs,
                                  ecm_set,
                                  all_atoms,
                                  all_bonds,
                                  all_dihedrals,
                                  atom_count,
                                  bond_count,
                                  dihe_count);
        }
        
        if (generalParameters.Num_of_Chromatins!=0) {
            OpenMM_Chromatin_info_relay(Chromatins,
                                        chromatin_set,
                                        all_atoms,
                                        all_bonds,
                                        all_angles,
                                        atom_count,
                                        bond_count,
                                        angle_count,
                                        interaction_map);
        }
        
        print_statistics(num_of_atoms,
                         num_of_bonds,
                         num_of_dihedrals,
                         num_of_angles,
                         Membranes,
                         Chromatins);
    
//    } End of Resume
    //autocorrelation calculations:
    //        GenConst::velocity_save.resize(6);
    
    
    // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
    // usage and runtime errors are caught and reported.
    
    try {
        
        int savetime = 0;
        MyOpenMMData* omm = new MyOpenMMData();
        TimeDependantData* time_dependant_data = new TimeDependantData();
        
//        if (!generalParameters.Resume) {
            omm = myInitializeOpenMM(all_atoms, generalParameters.Step_Size_In_Fs, platformName, time_dependant_data, all_bonds, all_dihedrals, all_angles, membrane_set, actin_set, ecm_set, chromatin_set, userinputs, interaction_map);
//        } else {
            
            
//        }
        
        
        
        if (!generalParameters.Resume) {
            assign_project_directories(buffer);
            cout<< "\nFile name: "<<TFILE<<generalParameters.trajectory_file_name<<TRESET<<endl<<endl;
        
        }
        string hardwarereportpath = generalParameters.trajectory_file_name+"_hardware_runtime.txt";
        if (!generalParameters.Resume) {
            ofstream write_hardwareReport;
            write_hardwareReport.open(hardwarereportpath.c_str(),std::ios_base::app);
            write_hardwareReport<<generalParameters.hardwareReport;
            write_hardwareReport<<endl;
        }
        
        string hardwareReportHeader;
        ifstream read_hardwareReport;
        read_hardwareReport.open(hardwarereportpath.c_str());
        string line;
        while(getline (read_hardwareReport,line)){
            hardwareReportHeader+=line+"\n";
        }
        if (generalParameters.Resume) {
            hardwareReportHeader+="\n Resume ================================================\n";
        }
        
        
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
        
        if (!generalParameters.Resume) {
        
            std::string Reporname=generalParameters.trajectory_file_name+"_report.txt";
            generate_report(config_lines);
        }
        string ckeckpoint_name=generalParameters.trajectory_file_name+"_Checkpoint";
        
        
        //Time the programme
        tStart = clock();
        chrono_clock_start = chrono::steady_clock::now();
        chrono_sys_clock_start = chrono::system_clock::now();
        
        
        
        //int SavingStep     = (int)(GenConst::Report_Interval_In_Fs / GenConst::Step_Size_In_Fs + 0.5);
        int MCCalcstep = (int)(GenConst::Mem_fluidity * generalParameters.Step_Size_In_Fs + 0.5);
        int NumSilentSteps = (int)(generalParameters.Report_Interval_In_Fs / generalParameters.Step_Size_In_Fs + 0.5);
        int Savingstep = NumSilentSteps;
        if ( (MCCalcstep < NumSilentSteps) && MCCalcstep != 0 ) {
            Savingstep = (int)(NumSilentSteps / MCCalcstep )* MCCalcstep;
            NumSilentSteps = MCCalcstep;
        } else if ( (NumSilentSteps < MCCalcstep) && MCCalcstep != 0 ){
            int rate =  (int)(MCCalcstep / NumSilentSteps );
            Savingstep = int(MCCalcstep/rate);
            NumSilentSteps = Savingstep;
        }
        
        
        int MCCalcTime = MCCalcstep;
        
        cout<<"Savingstep "<<Savingstep<<"\nNumSilentSteps "<<NumSilentSteps<<endl;
        
        int total_step_num = 0;
        
        bool expanding = false;
        bool set_spring = false;
        bool changeTemp = false;
        double TempStep = 0.2;
        double initTemp = generalParameters.temperature;
        
        if (!generalParameters.Resume) {
            if (generalParameters.WantPSF) {
                myWritePSF(num_of_atoms, num_of_bonds, all_atoms, all_bonds);
            }
            if (generalParameters.WantPDB) {
                myWritePDBFrame(0, 0, 0, 0, all_atoms, all_bonds);
            }
            if (generalParameters.WantXYZ) {
                writeXYZFrame(atom_count,all_atoms,0, 0, 0, false);
            }
            std::filebuf wfb;
            wfb.open (ckeckpoint_name.c_str(),std::ios::out);
            std::ostream wcheckpoint(&wfb);
            omm->context->createCheckpoint(wcheckpoint);
            wfb.close();
        } else {
            int infoMask = 0;
            infoMask = OpenMM::State::Positions;
            infoMask += OpenMM::State::Velocities;  // for kinetic energy (cheapm)
            const OpenMM::State state = omm->context->getState(infoMask);
            
            double timeInPs = state.getTime(); // OpenMM time is in ps already
            savetime=    int(timeInPs*1000/generalParameters.Step_Size_In_Fs);
            progressp =  int(1000*timeInPs/generalParameters.Simulation_Time_In_Ps)/10. + 0.1;
            
        }
        
//        cout<<progressp<<"   "<<savetime<<endl;exit(0);
        
        
        if (generalParameters.Minimise && !generalParameters.Resume) {
            string traj_name= generalParameters.trajectory_file_name+"_init.xyz";
            ofstream writexyz(traj_name.c_str());
            for (int n=0; all_atoms[n].type != -1; n++) {
                writexyz<<all_atoms[n].posInNm[0]<<"\t"<<all_atoms[n].posInNm[1]<<"\t"<<all_atoms[n].posInNm[2]<<"\n";
            }
            writexyz.close();
            minimisation(omm,
                         all_atoms,
                         all_bonds);
        }
        
        double time, energyInKJ, potential_energyInKJ;
        
        for (int frame=1; ; ++frame) {
            
            time=0; energyInKJ=0; potential_energyInKJ=0;
            
            myGetOpenMMState(omm->context, time, energyInKJ, potential_energyInKJ, all_atoms);
            
            if ( int(time*1000/generalParameters.Step_Size_In_Fs) > savetime ) {
//                Update_classes(Membranes, Actins, ECMs, Chromatins, time, all_atoms);
                
                collect_data(omm, all_atoms, interaction_map, Membranes, time);
                if (generalParameters.WantPDB) {
                    myWritePDBFrame(frame, time, energyInKJ, potential_energyInKJ, all_atoms, all_bonds);
                }
                if (generalParameters.WantXYZ) {
                    if (generalParameters.usingBackupCheckpoint) {
                        //This will only be run once after the checkpoint back up is loaded. generalParameters.usingBackupCheckpoint will become false after the new checkpoint is made.
                        writeXYZFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ,false);
                        
                    } else {
                        writeXYZFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ,true);
                    }
                }
                //Begin: Exporting congiguration of classes for simulation .
                
                string ckeckpoint_name_backup = ckeckpoint_name + "Backup";
                if (!generalParameters.usingBackupCheckpoint) {
                    boost::filesystem::copy(ckeckpoint_name, ckeckpoint_name_backup);
                } else {
                    generalParameters.usingBackupCheckpoint = false;
                }
                
                    
                
                std::filebuf wfb;
                wfb.open (ckeckpoint_name.c_str(),std::ios::out);
                std::ostream wcheckpoint(&wfb);
                omm->context->createCheckpoint(wcheckpoint);
                wfb.close();
                
                savetime += Savingstep;
                print_time(generalParameters.trajectory_file_name+"_hardware_runtime.txt", hardwareReportHeader, false,
                           tStart,
                           chrono_clock_start,
                           chrono_sys_clock_start
                           );
            }
            
            
//            if (changeTemp) {
//                initTemp -= TempStep;
//                omm->Lintegrator->setTemperature(initTemp);
//                //                    cout<<"temp: "<<omm->Lintegrator->getTemperature()<<endl;
//                //                    map <string, double> params;
//                //                    params = omm->context->getParameters();
//                //
//                //                    cout<<"\nframe: "<<frame<<endl;
//                //                    cout<<params.size()<<endl;
//                //                    for(auto elem : params)
//                //                    {
//                //                       cout << elem.first << " " << elem.second << "\n";
//                //                    }
//                //                    cout<<"\n";
//                //
//            }
            //                cout<<"CreateCheckpoint\n";
            
//            if (GenConst::CreateCheckpoint) {
            
                //End: Exporting congiguration of classes for simulation resume.
//            }
            
            if (time >= generalParameters.Simulation_Time_In_Ps){
                break;
            }
            //                cout<<"myStepWithOpenMM\n";
            
            myStepWithOpenMM(omm,time_dependant_data, all_atoms, NumSilentSteps, total_step_num);
            
            if (100*time/generalParameters.Simulation_Time_In_Ps>progressp){
                printf("[ %2.1f ] time: %4.1f Ps [out of %4.1f Ps]    \r",100*time/generalParameters.Simulation_Time_In_Ps, time, generalParameters.Simulation_Time_In_Ps);
                cout<< std::flush;
                //                    cout<<"[ "<<int(progressp*10)/10.0<<"% ]   \t time: "<<time<<" Ps [out of "<<GenConst::Simulation_Time_In_Ps<<" Ps]    \r" << std::flush;
                progressp =  int(1000*time/generalParameters.Simulation_Time_In_Ps)/10. + 0.1;
                //                    cout<<progressp<<endl;
//                progress++;
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
            
            if ( int(time*1000/generalParameters.Step_Size_In_Fs) >= MCCalcTime and GenConst::Mem_fluidity !=0){
                
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
        if (generalParameters.WantXYZ) {
            writeXYZFrame(atom_count,all_atoms,time, energyInKJ, potential_energyInKJ,true);
        }
        cout<<"[ 100% ]\t time: "<<generalParameters.Simulation_Time_In_Ps<<"Ps\n";
        
        print_time(generalParameters.trajectory_file_name+"_hardware_runtime.txt", hardwareReportHeader, true,
                   tStart,
                   chrono_clock_start,
                   chrono_sys_clock_start
                   );
//        print_wall_clock_time((double)((clock() - tStart)/CLOCKS_PER_SEC), true);
//        print_real_time(chrono_clock_start, chrono::steady_clock::now(), true);
//        print_system_time(chrono_sys_clock_start, chrono::system_clock::now(), true);
        

        
        
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


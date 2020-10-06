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
#include "Arg_pars.hpp"

//#include "Tests.hpp"

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
#pragma warning(disable:4996)   // sprintf is unsafe
#endif

#include "OpenMM.h"


namespace GenConst {
double Simulation_Time_In_Ps;
int    MD_traj_save_step;
double Report_Interval_In_Fs;
double Step_Size_In_Fs;
double MD_T;
double K;
int    MD_thrmo_step;
int    MC_step;
int    Mem_fluidity;
bool   Periodic_box;
double Lbox;
bool   Periodic_condtion_status;
int    Num_of_Membranes;
int    Num_of_Chromatins;
int    Num_of_Actins;
int    Num_of_ECMs;
string trajectory_file_name;
string ProjectName;
string force_file_name;
double Bussi_tau;
double Actin_Membrane_Bond_Coefficient;
bool   Interaction_map;
string Interaction_map_file_name;
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

bool   write_bonds_to_PDB;
bool   WantEnergy;
bool   WantForce;
bool   WriteVelocitiesandForces;
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
}


int main(int argc, char **argv)
{
    ArgStruct_Analysis args;
    args = cxxparser_analysis(argc, argv);
    
    // get the current time.
    time_t t = time(0);
    
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_time_%H_%M",now);
    
    string general_file_name="General_param_map.txt";
    
    vector<string> membrane_config_list;
    vector<string> chromatin_config_list;
    vector<string> actin_config_list;
    vector<string> ecm_config_list;
    vector<string> pointparticle_config_list;
    
    read_general_parameters(general_file_name, membrane_config_list, chromatin_config_list, actin_config_list, ecm_config_list, pointparticle_config_list);
    
    ofstream Trajectory;
    string traj_file_name;//="Results/"+GenConst::trajectory_file_name+buffer+".xyz";
    string ckeckpoint_name;//=GenConst::Checkpoint_path+GenConst::trajectory_file_name+buffer;
    
    
    vector<Membrane> Membranes;
    vector<std::set<int> > membrane_set;
    
    vector<Actin> Actins;
    vector<std::set<int> > actin_set;
    
    vector<ECM> ECMs;
    vector<std::set<int> > ecm_set;
    
    vector<Chromatin> Chromatins;
    //    vector<std::set<int> > chromatin_set;
    vector<vector<std::set<int> > > chromatin_set;
    
    
    bool Include_Membrane  = false;
    bool Include_Actin     = false;
    
    int num_of_bonds=0;
    
    Include_Membrane = true;
    GenConst::Num_of_Membranes=1;
    Membranes.resize(GenConst::Num_of_Membranes);
    membrane_set.resize(GenConst::Num_of_Membranes);
    
    
    if (Include_Membrane){
        if (Include_Actin){
            for (int i=0; i<GenConst::Num_of_Actins; i++) {
                for (int j=0; j<GenConst::Num_of_Membranes; j++) {
                    num_of_bonds        += Actins[i].return_num_of_actin_membrane_shared_nodes(j);
                }
                
            } //for (int i=0; i<GenConst::Num_of_Actins; i++)
        }
    } // End of if (Include_Membrane)
    
    
    
    cout<<TBOLD<<"Entering analysis mode:\n"<<TRESET;
    if (args.membane_labels.size()==0) {
        args.membane_labels.push_back(get_pdb_first_label(args.analysis_filename));
    }
    args.num_atoms_per_frame = get_pdb_num_of_atoms(args.analysis_filename);
    if (args.framelimits_end==0) {
        args.framelimits_end = get_pdb_num_of_frames(args.analysis_filename, args.num_atoms_per_frame);
    } else {
        //c++ arrays start from 0
        args.framelimits_end--;
    }
    if (args.framelimits_beg!=0) {
        //c++ arrays start from 0
        args.framelimits_beg--;
    }
    cout<<args.num_atoms_per_frame<<" ATOMs in each frame"<<endl;
    for (int i=0; i<args.output_filename.size(); i++) {
        
        if (args.analysis_dim == 3 ) {
            cout<<TPINK;
            cout<<"3D Analysis"<<endl;
            Membranes[i].import_pdb_frames(args, i);
            
            if (!GenConst::Testmode) {
                cout<<"num of frames = "<<args.framelimits_end-args.framelimits_beg<<endl;
            }
            for (int i=0; i<args.framelimits_end-args.framelimits_beg; i++) {
                
                Membranes[0].load_pdb_frame(i, args);
                
                for (int runs=0; runs<args.num_ang_avg; runs++) {
                    
                    Membranes[0].calculate_real_ulm(args);
                    //                Membranes[0].myWritePDBFrame(runs,temp_pdb_name+".pdb");
                    //                Membranes[0].calculate_ulm_sub_particles(ell_max, analysis_averaging_option);
                }
                
                cout<<"frame "<<i+args.framelimits_beg+1<<", End=[ "<<args.framelimits_end<<" ]\r"<< std::flush;
            }
            cout<<TRESET<<endl;
            Membranes[0].write_ulm(args, i);
            return 3;
        } else if(args.analysis_dim==2){
            cout<<TCYAN;
            cout<<"2D Analysis"<<endl;
            Membranes[i].import_pdb_frames(args, i);
            
            if (!GenConst::Testmode) {
                cout<<"num of frames = "<<args.framelimits_end-args.framelimits_beg<<endl;
            }
            
            for (int i=0; i<args.framelimits_end-args.framelimits_beg; i++) {
                Membranes[0].load_pdb_frame(i, args);
                
                
                
                cout<<"frame "<<i+args.framelimits_beg+1<<", End=[ "<<args.framelimits_end<<" ]\r"<< std::flush;
            }
            cout<<endl;
            return 2;
        }
    }
    
    
    
}


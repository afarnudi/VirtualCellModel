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

//#include "Tests.hpp"

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


bool   CreateCheckpoint;
bool   Load_from_checkpoint;
string Checkpoint_path;
string Checkpoint_file_name;
bool   ChromatinVirtualSites;


std::vector<double> data_colection_times;


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
    
    vector<Membrane> Membranes;
    
    generalParameters.Num_of_Membranes=int(args.Mesh_files.size());
    Membranes.resize(generalParameters.Num_of_Membranes);
    
    
    
    cout<<TBOLD<<"Entering analysis mode:\n"<<TRESET;
    if (args.membane_labels.size()==0) {
        args.membane_labels.push_back(get_pdb_first_label(args.analysis_filename));
        generalParameters.Num_of_Membranes=int(args.membane_labels.size());
        Membranes.resize(generalParameters.Num_of_Membranes);
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
    for (int memindex=0; memindex<args.output_filename.size(); memindex++) {
        
        if (args.analysis_dim == 3 ) {
            cout<<TPINK;
            cout<<"3D Analysis"<<endl;
            Membranes[memindex].import_pdb_frames(args, memindex);
            
            if (!generalParameters.Testmode) {
                cout<<"num of frames = "<<args.framelimits_end-args.framelimits_beg<<endl;
            }
            if (args.MeshMinimisation) {
                Membranes[memindex].get_ground_state_from_mesh(args);
            } else if (args.FrameMinimisation>-0.5) {
                Membranes[memindex].get_ground_state_from_frame(args);
                args.MeshMinimisation= true;
            }
            for (int frame=0; frame<args.framelimits_end-args.framelimits_beg; frame++) {
                
                Membranes[memindex].load_pdb_frame(frame, args);
                
                for (int runs=0; runs<args.num_ang_avg; runs++) {
                    
                    Membranes[memindex].calculate_real_ulm(args);
                    //                Membranes[memindex].myWritePDBFrame(runs,temp_pdb_name+".pdb");
                    //                Membranes[memindex].calculate_ulm_sub_particles(ell_max, analysis_averaging_option);
                }
                
                cout<<"frame "<<frame+args.framelimits_beg+1<<", End=[ "<<args.framelimits_end<<" ]\r"<< std::flush;
            }
            cout<<TRESET<<endl;
            Membranes[memindex].write_ulm(args, memindex);
            return 3;
        } else if(args.analysis_dim==2){
            cout<<TFILE;
            cout<<"2D Analysis"<<endl;
            Membranes[memindex].import_pdb_frames(args, memindex);
            if (!generalParameters.Testmode) {
                cout<<"num of frames = "<<args.framelimits_end-args.framelimits_beg<<endl;
            }
//
            for (int frame=0; frame<args.framelimits_end-args.framelimits_beg; frame++) {
                Membranes[memindex].load_pdb_frame(frame, args);
                
                Membranes[memindex].get_ring(args);

                cout<<"frame "<<frame+args.framelimits_beg+1<<", End=[ "<<args.framelimits_end<<" ]\r"<< std::flush;
            }
            cout<<endl;
            return 2;
        }
    }
    
    
    
}


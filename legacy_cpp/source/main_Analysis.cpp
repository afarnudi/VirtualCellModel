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
    
    if (args.undulateMesh) {
        if (args.extension=="") {
            args.extension="_undulated.ply";
        }
//        cout<<args.meshpathinput<<endl;
        
        
        string label=generalParameters.Membrane_label+to_string(0);
        Membranes[0].set_label(label);
        Membranes[0].set_file_time(buffer);
        Membranes[0].set_index(0);
        
        vector<string> configlines;
        configlines.push_back("-Membrane");
        configlines.push_back("MeshFile  " + args.meshpathinput);
        Membranes[0].import_config(configlines);
        
        int seed = args.SeedNumLmaxUmax[0];
        int numberOfRandomModes = args.SeedNumLmaxUmax[1];
        int ellMaxOfRandomModes = args.SeedNumLmaxUmax[2];
        double ulmOfRandomModes = args.SeedNumLmaxUmax[3];
        
        string filename = args.meshpathinput;
        filename.erase(filename.end()-4,filename.end() );
        args.output_filename.push_back(filename + "_S" + to_string(seed) + "N" + to_string(numberOfRandomModes) + "Ellmax" + to_string(ellMaxOfRandomModes) + "U" + to_string(ulmOfRandomModes) +args.extension);
//        cout<<args.output_filename[0]<<endl;
        srand(seed);
        
        Membranes[0].randomundulationgenerator(numberOfRandomModes,
                                               ellMaxOfRandomModes,
                                               ulmOfRandomModes);
        
//        Membranes[0].calculate_volume_and_surface_area();
//        double volume = Membranes[0].get_volume();
//        Membranes[0].rescale_membrane(cbrt(1./volume));
        
        ifstream read;
        read.open(args.meshpathinput.c_str());
        ofstream write;
        write.open(args.output_filename[0]);
        
        int temp_int; // This is just a temp intiger charachter that we use to read unnecessary ply file intigers. We never use these intigers in the actual programme.
        std::string temp_string;
        getline(read, temp_string);
        write<<temp_string<<endl;
        getline(read, temp_string);
        write<<temp_string<<endl;
        getline(read, temp_string);
        write<<temp_string<<endl;
        getline(read, temp_string);
        write<<temp_string<<endl;
        
        getline(read, temp_string);
        write<<temp_string<<endl;
        getline(read, temp_string);
        write<<temp_string<<endl;
        getline(read, temp_string);
        write<<temp_string<<endl;
        getline(read, temp_string);
        write<<temp_string<<endl;
        
        getline(read, temp_string);
        write<<temp_string<<endl;
        getline(read, temp_string);
        write<<temp_string<<endl;
        // In this section the Node coordinates are read from the ply file.
        
        for(int i=0;i<Membranes[0].get_num_of_nodes();i++)
        {
            getline(read, temp_string);
            write<< Membranes[0].get_node_position(i, 0) << " ";
            write<< Membranes[0].get_node_position(i, 1) << " ";
            write<< Membranes[0].get_node_position(i, 2) << " "<<endl;
        }
        
        while (getline(read, temp_string)) {
            write<<temp_string<<endl;
        }
        
        return 0;
    }
    
    cout<<TBOLD<<"Entering analysis mode:\n"<<TRESET;
    if (args.membane_labels.size()==0) {
        if (args.fileExtension=="pdb") {
            args.membane_labels.push_back(get_pdb_first_label(args.analysis_filename));
        } else if (args.fileExtension=="xyz") {
            args.membane_labels.push_back(get_xyz_first_label(args.analysis_filename));
        }
        
        if (args.analysis == "3") {
            if (args.extension=="") {
                args.extension="_ulmt.txt";
            }
        } else if (args.analysis == "2") {
            if (args.extension=="") {
                args.extension="_f2d.txt";
            }
        } else if (args.analysis == "E") {
            if (args.extension=="") {
                args.extension="_MechanicalEnergy.txt";
            }
        }
        
        string filename = args.analysis_filename;
        filename.erase(filename.end()-4,filename.end() );
        args.output_filename[0]=filename +"_"+args.membane_labels[0]+args.extension;
        
        generalParameters.Num_of_Membranes=int(args.membane_labels.size());
        Membranes.resize(generalParameters.Num_of_Membranes);
    }
    if (args.fileExtension=="pdb") {
        args.num_atoms_per_frame = get_pdb_num_of_atoms(args.analysis_filename);
        if (args.framelimits_end==0) {
            args.framelimits_end = get_pdb_num_of_frames(args.analysis_filename, args.num_atoms_per_frame);
        }
    } else if (args.fileExtension=="xyz") {
        args.num_atoms_per_frame = get_xyz_num_of_atoms(args.analysis_filename);
        if (args.framelimits_end==0) {
            args.framelimits_end = get_xyz_num_of_frames(args.analysis_filename, args.num_atoms_per_frame);
            args.framelimits_end++;
        }
//        cout<<args.framelimits_end<<endl;exit(0);
    }
    
    
    if (args.framelimits_end!=0) {
        //c++ arrays start from 0
        args.framelimits_end--;
    }
    if (args.framelimits_beg!=0) {
        //c++ arrays start from 0
        args.framelimits_beg--;
    }
    
    cout<<args.num_atoms_per_frame<<" ATOMs in each frame"<<endl;
    for (int memindex=0; memindex<args.output_filename.size(); memindex++) {
        
        if (args.analysis == "3" ) {
            cout<<TPINK;
            cout<<"3D Analysis"<<endl;
            if (args.fileExtension=="pdb") {
                Membranes[memindex].import_pdb_frames(args, memindex);
            } else if (args.fileExtension=="xyz") {
                Membranes[memindex].import_xyz_frames(args, memindex);
            }
            
            if (!generalParameters.Testmode) {
                cout<<args.framelimits_end<<" "<<args.framelimits_beg<<endl;
                cout<<"num of frames = "<<args.framelimits_end-args.framelimits_beg<<endl;
            }
            if (args.MeshMinimisation) {
                Membranes[memindex].get_ground_state_from_mesh(args);
            } else if (args.FrameMinimisation>-0.5) {
                Membranes[memindex].get_ground_state_from_frame(args);
                args.MeshMinimisation= true;
            }
            for (int frame=0; frame<args.framelimits_end-args.framelimits_beg; frame++) {
                Membranes[memindex].load_analysis_coord_frame(frame, args);
                for (int runs=0; runs<args.num_ang_avg; runs++) {
                    
                    Membranes[memindex].calculate_real_ulm(args);
                    //                Membranes[memindex].myWritePDBFrame(runs,temp_pdb_name+".pdb");
                    //                Membranes[memindex].calculate_ulm_sub_particles(ell_max, analysis_averaging_option);
                }
                
                cout<<"frame "<<frame+args.framelimits_beg+1<<", End=[ "<<args.framelimits_end<<" ]\r"<< std::flush;
            }
            
            Membranes[memindex].write_ulm(args, memindex);
            
            cout<<TBOLD<<"\n\nAnalysis finished successfully."<<TRESET<<endl;
            
        } else if(args.analysis== "2"){
            cout<<TFILE;
            cout<<"2D Analysis"<<endl;
            if (args.fileExtension=="pdb") {
                Membranes[memindex].import_pdb_frames(args, memindex);
            } else if (args.fileExtension=="xyz") {
                Membranes[memindex].import_xyz_frames(args, memindex);
            }
            if (!generalParameters.Testmode) {
                cout<<"num of frames = "<<args.framelimits_end-args.framelimits_beg<<endl;
            }
            
//            Membranes[memindex].updatepos(100);
            Membranes[memindex].get_ring(args);
            
            for (int frame=0; frame<args.framelimits_end-args.framelimits_beg; frame++) {
                Membranes[memindex].load_analysis_coord_frame(frame, args);
                
                Membranes[memindex].calculate_freqs(args);
                
//                Membranes[memindex].calculate_freqs_usingSH(args);
//                Membranes[memindex].calculate_freqs_alexandra(args);

//                cout<<"frame "<<frame+args.framelimits_beg+1<<", End=[ "<<args.framelimits_end<<" ]\r"<< std::flush;
            }
            Membranes[memindex].write_un_uq(args, memindex);
            
            
            for (int frame=0; frame<args.framelimits_end-args.framelimits_beg; frame++) {
                Membranes[memindex].load_analysis_coord_frame(frame, args);
                
               
                Membranes[memindex].calculate_freqs_projection(args);
//                Membranes[memindex].calculate_freqs_usingSH(args);
//                Membranes[memindex].calculate_freqs_alexandra(args);

//                cout<<"frame "<<frame+args.framelimits_beg+1<<", End=[ "<<args.framelimits_end<<" ]\r"<< std::flush;
            }
            
            Membranes[memindex].write_un_uq_projection(args, memindex);
//            Membranes[memindex].write_uq_SH(args, memindex);
//            Membranes[memindex].write_uq_Alexandra(args, memindex);
            cout<<TBOLD<<"\n\nAnalysis finished successfully."<<TRESET<<endl;
            
        } else if (args.analysis == "E"){
            cout<<TPURPLE;
            cout<<"Mechanical Energy Analysis"<<endl;

            map<string, vector<string> > config_lines = read_configfile(args.Mesh_files[memindex]);
            get_class_numbers(config_lines);
            vector<vector<string> > membrane_configs =sort_class_configs(config_lines["-Membrane"]);
            Membranes[memindex].import_config(membrane_configs[ args.membane_label_indecies[memindex] ]);
            cout<<TPURPLE;
            Membranes[memindex].import_pdb_frames(args, memindex);
            Membranes[memindex].load_analysis_coord_frame(0, args);
            Membranes[memindex].set_bond_nominal_length();
            Membranes[memindex].set_bending_nominal_angle();
            FILE* pFile;
            pFile = fopen (args.output_filename[memindex].c_str(),"w");
            fprintf(pFile,"time ps;   BendingEnergy      BondEnergy     TotalEnergy KJ/mole \n");
            fclose (pFile);
            
            
            if (!generalParameters.Testmode) {
                cout<<"num of frames = "<<args.framelimits_end-args.framelimits_beg<<endl;
            }
            for (int frame=0; frame<args.framelimits_end-args.framelimits_beg; frame++) {
                
                Membranes[memindex].load_analysis_coord_frame(frame, args);
                
                Membranes[memindex].Mechanical_Energy_calculator();
                Membranes[memindex].Export_Mechanical_Energy(args.output_filename[memindex], frame+args.framelimits_beg);
                
                cout<<"frame "<<frame+args.framelimits_beg+1<<", End=[ "<<args.framelimits_end<<" ]\r"<< std::flush;
            }
            
            cout<<TBOLD<<"\n\nAnalysis finished successfully."<<TRESET<<endl;
            
            
        }
    }

    return 0;
    
}


//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
#include <iostream>
#include <fstream>
#include "General_constants.h"
//#include <stdlib.h>

#include "Arg_pars.hpp"
#include "cxxopts.hpp"

using namespace std;

void consistency_check(ArgStruct &args);
void build_output_names(ArgStruct &args);
void check_membrane_labels(ArgStruct &args);


ArgStruct cxxparser(int argc, char **argv){
    ArgStruct args;
    string programme_name =argv[0];
    try
    {
        
        cxxopts::Options options(programme_name,
                                 "\n Welcome to VCM.\n\n"
                                 " If the configuration files are ready you can run the simulation by running the"
                                 " programme without any option arguments by following this example:"
                                 " \n"
                                 " "+programme_name+"\n"
                                 " \n"
                                 " How to utilise VCM's analysis options:\n"
                                 " 3D analysis of the Membrane surface undulations:\n"
                                 " \n "+programme_name+
                                 " --analysis 3"
                                 " --pdbpath <path to pdb>"
                                 " --lmax <maximum mode number>"
                                 " --angavg <number of random rotations>"
                                 " --align_axes <node on the Z axis>,<node on the Z-Y plane>"
                                 " --ext <output file extension>"
                                 " --framelimits <first frame>,<last frame>"
                                 " --memlabels <firstlabel>,<secondlabel>,...");
        options
        .positional_help("")
        .show_positional_help();
        
        options
        //                .allow_unrecognised_options()
        .add_options()
        ("h,help", "Print help")
        ("analysis", "3 for 3D analysis and 2 for 2D", cxxopts::value<int>(),"int")
        ("pdbfile", "Path to Membrane pdb file that contain ONLY a single Membrane. Example: path/to/my/pdbfile.pdb", cxxopts::value<std::vector<std::string>>(),"Path+file")
        
        ("lmax", "The maximum mode number (l) the RSH amplitudes are measured for. This is a very expensive analysis since the number of angular modes (m) grow very rapidly (2*l+1) with l. Default 20.", cxxopts::value<int>(), "int")
        
        ("angavg", "[optional] Determin the number of random orientations the amplitude of a U_lm will be measured for. For each pdb frame, the Membrane will be randomly rotated (with respect to the Membrane's centre of mass) and the average amplitude of the measured U_lms will be written to the output file. Default 0. Note: This option cannot be selected simultaniusly with \"--align_axes\".", cxxopts::value<int>(), "int")
        
        ("alignaxes", "[optional] Takes two integers, \"first index\" and \"second index\" to specify the orientation of the material frame with respect to which the RSH analysis will be done. During the analysis, the origin of the reference frame will be at the Membrane centre of mass and the frame of reference will be rotated so that the Z axis goes through the first node (first index) and the second node (second index) lies on the Z-Y plane. Note: This option cannot be selected simultaniusly with \"--angavg\".", cxxopts::value<std::vector<int>>(),"int,int")
        
        ("ext", "[optional] The Ulm amplitudes will be saved in a file with the same name as the pdb file provided in the \"pdbpath\" + what you put in after this flag. Example (Note: irelevant flags are omited):   ./VCM --pdbpath mydirectory/mypdb.pdb --ext _myUlms.myextension --memlabels mem0                 Output files: mydirectory/mypdb_mem0_myUlms.myextension", cxxopts::value<std::string>(),"example.txt")
        
        ("framelimits", "[optional] Takes the start and end frame of the pdb. The analysis will be performed only for these frames and all pdb frames in between. A 0 value for the first or last frame will be interpreted as the biginning frame and the las frame respectfully. Default is 0 for both limits. Example (Note: irelevant flags are omited):"
         "                                             "
         "\n  1)  ./VCM --framelimits 0,0""                 "
         "\n  2)  ./VCM --framelimits 3,0""                 "
         "\n  3)  ./VCM --framelimits 0,293""               "
         "\n  4)  ./VCM --framelimits 11,73""               "
         "\n  interpretation:""                             "
         "\n  1)  All frames.""                             "
         "\n  2)  From frame 3 to the end.""                "
         "\n  3)  From the beginning upto frame 293         "
         "(Not including 293)                               "
         "\n  4)  From frame 11 upto frame 73 (Not including"
         "73)                                               "
         , cxxopts::value<std::vector<int>>(),"int,int")
        ("memlabels", "The label(s) used to represent the memebrane(s) in the pdb file. The label(s) will also be used to distinguish between output files. Example (Note: irelevant flags are omited):               ./VCM --pdbpath mydirectory/mypdb.pdb --ext _myUlms.myextension --memlabels mem0,mem1            Output files: mydirectory/mypdb_mem0_myUlms.myextension and mydirectory/mypdb_mem1_myUlms.myextension", cxxopts::value<std::vector<std::string>>(),"mem0")
        ;
        //      I have put this here so that parser throuws an exception when an extra argument is put in the command line that is not associated with any flags
        options.parse_positional({""});
        
        auto result = options.parse(argc, argv);
        
        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            exit(0);
        }
        
        if (result.count("analysis"))
        {
            args.analysis_dim=result["analysis"].as<int>();
        }
        
        if (result.count("pdbfile"))
        {
            auto& ff = result["pdbfile"].as<std::vector<std::string>>();
            string filename;
            for (const auto& f : ff)
            {
                filename+= f ;
            }
            args.analysis_filename = filename;
        }
        if (result.count("lmax"))
        {
            args.ell_max = result["lmax"].as<int>();
        }
        
        if (result.count("angavg"))
        {
            args.analysis_averaging_option = 1;
            args.num_ang_avg = result["angavg"].as<int>();
        }
        if (result.count("alignaxes"))
        {
            args.analysis_averaging_option=2;
            const auto values = result["alignaxes"].as<std::vector<int>>();
            if (values.size()!=2) {
                cout<<TWWARN<<"Error: "<<TRESET<<"Expected 2 int arguments for \"alignaxes\", got "<<values.size()<<"\n Input format example: --alignaxes 2,4\nRun help (-h, --help) for more information."<<endl;
                exit(1);
            }
            args.z_node = values[0];
            args.zy_node= values[1];
        }
        if (result.count("ext"))
        {
            args.extension = result["ext"].as<std::string>();
        }
        if (result.count("framelimits"))
        {
            const auto values = result["framelimits"].as<std::vector<int>>();
            if (values.size()!=2) {
                cout<<TWWARN<<"Error: "<<TRESET<<"Expected 2 int arguments for \"framelimits\", got "<<values.size()<<"\n Run help (-h, --help) for more information."<<endl;
                exit(1);
            } else if (values[0] > values[1] && values[1]!=0){
                cout<<TWWARN<<"Error: "<<TRESET"Got "<<values[0]<<" as the beginning frame and "<<values[1]<<" as the end frame. Beginning frame cannot be larger than end frame.\n Run help (-h, --help) for more information."<<endl;
                exit(1);
            } else if (values[0]!=0 && values[1]!=0 ){
                if (values[0]==values[1]) {
                    cout<<TWWARN<<"Error: "<<TRESET<<"Beginning frame ("<<values[0]<<") and end frame ("<<values[0]<<") cannot point to the same frame number. If you wish to analyse frame number N, set the beginning frame to N and end frame to N+1.\n Run help (-h, --help) for more information."<<endl;
                    exit(1);
                }
            }
            args.framelimits_beg = values[0];
            args.framelimits_end = values[1];
        }
        if (result.count("memlabels"))
        {
            auto& ff = result["memlabels"].as<std::vector<std::string>>();
            for (const auto& f : ff)
            {
                args.membane_labels.push_back(f) ;
            }
        }
        
        consistency_check(args);
        return args;
        
    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(6);
    }
    
    
}

void consistency_check(ArgStruct &args){
    cout<<endl;
    
    if (args.analysis_dim != 2 && args.analysis_dim !=3) {
        cout<<"Analysis dimension option "<<TWWARN<<"not supported"<<TRESET<<".\nUse -h for more information"<<endl;
        exit(2);
    } else {
        if (args.analysis_dim == 3) {
            cout<<TWARN<<"3D"<<TRESET" analysis mode: "<<TON<<"ON\n"<<TRESET;
            if (args.analysis_filename == "" ) {
                cout<<"Analysis file (path+file) "<<TWWARN<<"not provided"<<TRESET<<". Use -h for more information.\n";
                exit(3);
            } else {
                cout<<"Will analyse \""<<TFILE<<args.analysis_filename<<TRESET<<"\""<<endl;
            }
            cout<<"Real Spherical harmonics, U_lm, will go up to l="<<args.ell_max<<endl;
            
            if (args.num_ang_avg !=1 && args.z_node !=-1) {
                cout<<"\"angavg\" and \"alignaxes\" "<<TWWARN<<"cannot be used simultaneously"<<TRESET<<". Run help (-h, --help) for more inforamtion."<<endl;
                exit(4);
            }else {
                if (args.analysis_averaging_option == 1) {
                    cout<<"Average over random rotations:"<<TON<<" ON"<<TRESET<<endl;
                    cout<<"Number of random rotations: "<<args.num_ang_avg<<endl;
                }else if (args.analysis_averaging_option == 2) {
                    cout<<"Material frame reference: "<<TON<<"ON"<<TRESET<<endl;
                    cout<<"Node index on the Z axis: "<<args.z_node<<endl;
                    cout<<"Node inedx on the Z-Y plane: "<<args.zy_node<<endl;
                    
                }
            }
            
        }
        if (args.analysis_dim == 2) {
            cout<<TWARN<<"2D"<<TRESET" analysis mode: "<<TON<<"ON\n"<<TRESET;
            //                cout<<"This option is not available at the moment\nUse -h for more information.\n";
            //                exit(5);
            
            cout<<TWWARN<<"This option is under development. Prceed with the trial version.\n"<<TRESET;
            
            if (args.analysis_filename == "" ) {
                cout<<"Analysis file (path+file) "<<TWWARN<<"not provided"<<TRESET<<". Use -h for more information.\n";
                exit(3);
            } else {
                cout<<"Will analyse \""<<TFILE<<args.analysis_filename<<TRESET<<"\""<<endl;
            }
            
            
            
            
        }
        
        if (args.membane_labels.size()==0) {
            cout<<TWARN<<"\n!!!Warning"<<TRESET<<", no membrane labels assigned. If the pdb strictly containes one Membrane class (and no other classes), ignore this warning. If it contains multiple membranes or other classes as well, assign the membrane label to be analysed.\n"<<endl;
        } else {
            check_membrane_labels(args);
        }
        build_output_names(args);
        
        cout<<"Output file(s):\n";
        for (int i=0; i<args.output_filename.size(); i++) {
            cout<<"\t"<<TFILE<<args.output_filename[i]<<TRESET<<endl;
        }
        
        if (args.analysis_dim==3) {
            string meshpath;
            cout<<"Mesh information is "<<TWARN<<"required"<<TRESET<<" for 3D analysis."<<endl;
            if (args.membane_labels.size()<2) {
                cout<<"Please enter the path+file for the Membrane mesh.\nExample:\n\tpath/to/my/mesh.ply"<<endl;
                cout<<TFILE;
                cin>>meshpath;
                cout<<TRESET;
                std::ifstream read_mesh;
                read_mesh.open(meshpath.c_str());
                while (!read_mesh) {
                    std::cout << TWWARN<<"Unable to read"<<TRESET<<" \""<<TFILE<<meshpath<<TRESET<<"\" or it does not exist.\nPlease try again."<<std::endl;
                    cin>>meshpath;
                    read_mesh.open(meshpath.c_str());
                }
                args.Mesh_files.push_back(meshpath);
            } else {
                cout<<"Please enter the path+file for the following Membranes.\nExample:\n\tpath/to/my/mesh.ply"<<endl;
                for (int i=0; i<args.membane_labels.size(); i++) {
                    cout<<"For "<<TWARN<<args.membane_labels[i]<<TRESET<<":"<<endl;
                    cout<<TFILE;
                    cin>>meshpath;
                    cout<<TRESET;
                    std::ifstream read_mesh;
                    read_mesh.open(meshpath.c_str());
                    if (!read_mesh) {
                        std::cout << TWWARN<<"Unable to read"<<TRESET<<" \""<<TFILE<<meshpath<<TRESET<<"\" or it does not exist.\nPlease try again."<<std::endl;
                        i--;
                    } else {
                        args.Mesh_files.push_back(meshpath);
                    }
                }
            }
            
        }
        
        if (args.framelimits_beg==0 && args.framelimits_end==0) {
            cout<<"All frames will be analysied."<<endl;
        } else if (args.framelimits_beg==0) {
            cout<<"Frames from the beginning to frame "<<args.framelimits_end<<" will be analysed."<<endl;
            
        } else if (args.framelimits_end==0) {
            cout<<"Frames from "<<args.framelimits_beg<<" to the end will be analysed."<<endl;
        } else if (args.framelimits_end == args.framelimits_beg+1){
            cout<<"Frame "<<args.framelimits_beg<<" will be analysed."<<endl;
        } else {
            cout<<"Frames from "<<args.framelimits_beg<<" to "<<args.framelimits_end<<" will be analysed."<<endl;
        }
    }
    cout<<endl;
    
    
}

void build_output_names(ArgStruct &args){
    if (args.analysis_dim ==3) {
        if (args.extension=="") {
            args.extension="_ulmt.txt";
        }
    } else if (args.analysis_dim ==2) {
        if (args.extension=="") {
            args.extension="_f2d.txt";
        }
    }
    if (args.membane_labels.size()==0) {
        string filename = args.analysis_filename;
        filename.erase(filename.end()-4,filename.end() );
        args.output_filename.push_back(filename + args.extension);
    } else {
        string filename = args.analysis_filename;
        filename.erase(filename.end()-4,filename.end() );
        for (int i=0; i<args.membane_labels.size(); i++) {
            args.output_filename.push_back(filename +"_"+args.membane_labels[i]+args.extension);
        }
    }
    
}

void check_membrane_labels(ArgStruct &args){
    if (args.membane_labels.size()!=0) {
        for (int i=0; i<args.membane_labels.size()-1; i++) {
            string label = args.membane_labels[i];
            for (int j=i+1; j<args.membane_labels.size(); j++) {
                if (label==args.membane_labels[j]) {
                    args.membane_labels.erase(args.membane_labels.begin()+j);
                    j--;
                }
            }
            
        }
        for (int i=0; i<args.membane_labels.size(); i++) {
            if (args.membane_labels[i].size()!=4) {
                cout<<TWWARN<<"\n!!!Warning"<<TRESET<<", pdb labels are 4 charachters long. The Membrane label \""<<TWARN<<args.membane_labels[i]<<TRESET<<"\" is "<<args.membane_labels[i].size()<<" long. If there is a typo, end the programme and start again.\n\n";
            }
        }
    }
}

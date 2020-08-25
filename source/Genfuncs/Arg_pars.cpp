//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
#include <iostream>
//#include <stdlib.h>

#include "Arg_pars.hpp"
#include "cxxopts.hpp"

using namespace std;

void print_help(void);
void print_report_on_screen(ArgStruct args);


ArgStruct cxxparser(int argc, char **argv){
    ArgStruct args;
    string programme_name =argv[0];
    try
    {
        
        cxxopts::Options options(argv[0],
                                 "\n Welcome to VCM.\n\n"
                                 " If the configuration files are ready you can run the simulation by running the"
                                 " programme without any option arguments by following this example:"
                                 " \n"
                                 " "+programme_name+"\n"
                                 " \n"
                                 " How to utilise VCM's analysis options:\n"
                                 " 3D analysis of the Membrane surface undulations:\n"
                                 " \n "+programme_name+" -analysis 3 -pdbpath <path to pdb> -lmax <maximum mode number> -angavg "
                                 " <number of random rotations> -align_axes <node on the Z axis> <node on the "
                                 " Z-Y plane> -ext <output file extension> -framelimits <first frame> "
                                 " <last frame>");
        options
        .positional_help("")
        .show_positional_help();
        
        options
        //        .allow_unrecognised_options()
        .add_options()
        ("a,analysis", "Takes 2 or 3 to perform Membrane mode analysis in 3D or 2D. 3D analysis decomposes Membrnae surface undulations into Real Spherical Harmonic (RSH) modes, U_lm. l and m are integers where l is defined between 0 and lmax, and m between -l and l.", cxxopts::value<int>(), "int")
        
        ("pdbfile", "Path to Membrane pdb file that contain ONLY a single Membrane. Example: path/to/my/pdbfile.pdb", cxxopts::value<std::vector<std::string>>(),"Path+file")
        
        ("lmax", "The maximum mode number (l) the RSH amplitudes are measured for. This is a very expensive analysis since the number of angular modes (m) grow very rapidly (2*l+1) with l. Default 20.", cxxopts::value<int>(), "int")
        
        ("angavg", "[optional] Determin the number of random orientations the amplitude of a U_lm will be measured for. For each pdb frame, the Membrane will be randomly rotated (with respect to the Membrane's centre of mass) and the average amplitude of the measured U_lms will be written to the output file. Default 0. Note: This option cannot be selected simultaniusly with \"-align_axes\".", cxxopts::value<int>(), "int")
        
        ("alignaxes", "[optional] Takes two integers, \"first index\" and \"second index\" to specify the orientation of the material frame with respect to which the RSH analysis will be done. During the analysis, the origin of the reference frame will be at the Membrane centre of mass and the frame of reference will be rotated so that the Z axis goes through the first node (first index) and the second node (second index) lies on the Z-Y plane. Note: This option cannot be selected simultaniusly with \"-angavg\".", cxxopts::value<std::vector<int>>(),"int,int")
        
        ("ext", "The Ulm amplitudes will be saved in a file with the same name as the pdb file provided in the \"pdbpath\" + what you put in after this flag. Example (Note: irelevant flags are omited): ./VCM -pdbpath mydirectory/mypdb.pdb -ext _myUlms.myextension Output file: mydirectory/mypdb_myUlms.myextension", cxxopts::value<std::vector<std::string>>(),"example.txt")
        
        ("framelimits", "[optional] Takes the start and end frame of the pdb. The analysis will be performed only for these frames and all pdb frames in between. A 0 value for the first or last frame will be interpreted as the biginning frame and the las frame respectfully. Default is 0 for both limits. Example (Note: irelevant flags are omited):"
         "                                             "
         "\n  1)  ./VCM -framelimits 0,0""                  "
         "\n  2)  ./VCM -framelimits 3,0""                  "
         "\n  3)  ./VCM -framelimits 0,293""                "
         "\n  4)  ./VCM -framelimits 11,73""                "
         "\n  interpretation:""                             "
         "\n  1)  All frames.""                             "
         "\n  2)  From frame 3 to the end.""                "
         "\n  3)  From the beginning upto and including fram"
         "e 293                                             "
         "\n  4)  From and including frame 11 upto and inclu"
         "ding frame 73", cxxopts::value<std::vector<int>>(),"int,int")
        
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
            args.analysis_dim = result["analysis"].as<int>();
            args.analysis_mode = true;
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
                cout<<"Error: Expected 2 int arguments for \"alignaxes\", got "<<values.size()<<"\n Input format example: --alignaxes 2,4\nRun help (-h, --help) for more information."<<endl;
                exit(1);
            }
            args.z_node = values[0];
            args.zy_node= values[1];
        }
        if (result.count("ext"))
        {
            auto& ff = result["ext"].as<std::vector<std::string>>();
            string filename;
            for (const auto& f : ff)
            {
                filename+= f ;
            }
            args.output_filename = filename;
        }
        if (result.count("framelimits"))
        {
            const auto values = result["framelimits"].as<std::vector<int>>();
            if (values.size()!=2) {
                cout<<"Error: Expected 2 int arguments for \"framelimits\", got "<<values.size()<<"\n Run help (-h, --help) for more information."<<endl;
                exit(1);
            } else if (values[0] > values[1]){
                cout<<"Error: Got "<<values[0]<<" as the beginning frame and "<<values[1]<<" as the end frame. Beginning frame cannot be larger than end frame.\n Run help (-h, --help) for more information."<<endl;
                exit(1);
            }
            args.framelimits_beg = values[0];
            args.framelimits_end = values[1];
        }
        
        
        print_report_on_screen(args);
        return args;
        
    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(6);
    }
    
    
}

void print_report_on_screen(ArgStruct args){
    cout<<endl;
    if (args.analysis_mode) {
        if (args.analysis_dim != 2 && args.analysis_dim !=3) {
            cout<<"Analysis dimension option not understood.\nUse -h for more information"<<endl;
            exit(2);
        } else {
            if (args.analysis_dim == 3) {
                cout<<"3D analysis mode: ON\n";
                if (args.analysis_filename == "" ) {
                    cout<<"Analysis file (path+file) not provided. Use -h for more information.\n";
                    exit(3);
                } else {
                    cout<<"Will analyse \""<<args.analysis_filename<<"\""<<endl;
                }
                cout<<"Real Spherical harmonics, U_lm, will go up to l="<<args.ell_max<<endl;
                
                if (args.num_ang_avg !=1 && args.z_node !=-1) {
                    cout<<"\"angavg\" and \"alignaxes\" cannot be used simultaneously. Run help (-h, --help) for more inforamtion."<<endl;
                    exit(4);
                }else {
                    if (args.analysis_averaging_option == 1) {
                        cout<<"Average over random rotations: ON"<<endl;
                        cout<<"Number of random rotations: "<<args.num_ang_avg<<endl;
                    }else if (args.analysis_averaging_option == 2) {
                        cout<<"Material frame reference: ON"<<endl;
                        cout<<"Node index on the Z axis: "<<args.z_node<<endl;
                        cout<<"Node inedx on the Z-Y plane: "<<args.zy_node<<endl;
                        
                    }
                }
                string filename = args.analysis_filename;
                string extension = args.output_filename;
                filename.erase(filename.end()-4,filename.end() );
                args.output_filename = filename + extension;
                cout<<"Output file extension: "<<args.output_filename<<endl;
                if (args.framelimits_beg==0 && args.framelimits_end==0) {
                    cout<<"All frames will be analysied."<<endl;
                } else if (args.framelimits_beg==0) {
                    cout<<"Frames from the beginning to frame "<<args.framelimits_end<<" will be analysed."<<endl;
                    
                } else if (args.framelimits_end==0) {
                    cout<<"Frames from "<<args.framelimits_beg<<" to the end will be analysed."<<endl;
                } else {
                    cout<<"Frames from "<<args.framelimits_beg<<" to "<<args.framelimits_end<<" will be analysed."<<endl;
                }
            }
            if (args.analysis_dim == 2) {
                cout<<"2D analysis mode: ON\n";
                cout<<"This option is not available at the moment\nUse -h for more information.\n";
                exit(5);
            }
        }
        cout<<endl;
        
    }
}


ArgStruct arg_parse_old(int argc, char **argv){
    ArgStruct args;
    string read;
    for (int i=1; i<argc; i++) {
        string argvcase = argv[i];
        if ( (argvcase == "-h") || (argvcase == "-help") || (argvcase == "--h") || (argvcase == "--help")) {
            print_help();
            exit(0);
        } else if (argvcase == "-analysis"){
            if (i+1<argc) {
                read = argv[i+1];
                args.analysis_dim = stoi(read);
                args.analysis_mode = true;
                i++;
            }
        } else if (argvcase == "-pdbpath"){
            args.analysis_filename = argv[i+1];
            
            
            string extension = "_ulmt_cpp.txt";
            string filename = args.analysis_filename;
            filename.erase(filename.end()-4,filename.end() );
            args.output_filename = filename + extension;
            i++;
        } else if (argvcase == "-lmax"){
            read = argv[i+1];
            args.ell_max = stoi(read);
            cout<<"lmax "<<args.ell_max<<endl;
            i++;
        } else if (argvcase == "-angavg"){
            args.analysis_averaging_option = 1;
            cout<<"Average over random rotations: ON"<<endl;
            read = argv[i+1];
            args.num_ang_avg = stoi(read);
            cout<<"Number of random rotations: "<<args.num_ang_avg<<endl;
            i++;
        } else if (argvcase == "-align_axes"){
            args.analysis_averaging_option=2;
            cout<<"Material frame reference: ON"<<endl;
            read = argv[i+1];
            args.z_node = stoi(read);
            cout<<"Node index on the Z axis: "<<args.z_node<<endl;
            read = argv[i+2];
            args.zy_node = stoi(read);
            cout<<"Node inedx on the Z-Y plane: "<<args.zy_node<<endl;
            i+=2;
        } else if (argvcase == "-ext"){
            read = argv[i+1];
            string filename = args.analysis_filename;
            filename.erase(filename.end()-4,filename.end() );
            args.output_filename = filename + read;
            cout<<"Output file extension: "<<args.output_filename<<endl;
            i++;
        } else if (argvcase == "-framelimits"){
            read = argv[i+1];
            args.framelimits_beg = stoi(read);
            read = argv[i+2];
            args.framelimits_end = stoi(read);
            cout<<"Will analyse frames "<<args.framelimits_beg<<" through "<<args.framelimits_end<<endl;
            i+=2;
        }
    }
    if (args.analysis_mode) {
        if (args.analysis_dim != 2 && args.analysis_dim !=3) {
            cout<<"Analysis dimension option not understood.\nUse -h for more information"<<endl;
            exit(3);
        } else {
            if (args.analysis_dim == 3) {
                cout<<"3D analysis mode: ON\n";
                if (args.analysis_filename == "" ) {
                    cout<<"Analysis file (path+file) not provided. Use -h for more information.\n";
                    exit(2);
                } else {
                    cout<<"Will analyse \""<<args.analysis_filename<<"\""<<endl;
                }
                cout<<"Real Spherical harmonics, U_lm, will go up to l="<<args.ell_max<<endl;
            }
            if (args.analysis_dim == 2) {
                cout<<"2D analysis mode: ON\n";
                cout<<"This option is not available at the moment\nUse -h for more information.\n";
                exit(1);
            }
        }
        cout<<endl;
        
    }
    
    return args;
}


void print_help(void){
    cout<<
    "\n Welcome to VCM.\n\n"
    " If the configuration files are ready you can run the simulation by running the\n"
    " programme without any arguments by following this example:\n"
    " \n"
    " ./VCM\n"
    " \n"
    " How to utilise VCM's analysis options:\n"
    " \n"
    " 3D analysis of the Membrane surface undulations:\n"
    " ./VCM -analysis 3 -pdbpath <path to pdb> -lmax <maximum mode number> -angavg \n"
    " <number of random rotations> -align_axes <node on the Z axis> <node on the \n"
    " Z-Y plane> -ext <output file extension> -framelimits <first frame> \n"
    " <last frame>\n"
    "\n"
    " -analysis  Takes integers 2 or 3 to perform Membrane mode analysis in 3D or 2D.\n"
    "            3D analysis decomposes Membrnae surface undulations into Real Spher-\n"
    "            ical Harmonic (RSH) modes, U_lm. l and m are integers where l is \n"
    "            defined between 0 and lmax, and m between -l and l.\n"
    " \n"
    " -pdbpath   Path to Membrnae pdb file that contain ONLY a single Membrane.\n"
    " \n"
    " -lmax      An intiger determining the maximum mode (l) the RSH amplitudes\n"
    "            are measured for. This is a very expensive analysis since the \n"
    "            number of angular modes (m) grow very rapidly with l. Default 20.\n"
    " \n"
    " -angavg    [optional] An integer determining the number of random orientations \n"
    "            the amplitude of a U_lm will be measured for. For each pdb frame, \n"
    "            the Membrane is a randomply rotated in space (with respect to the \n"
    "            membrane centre of mass) and the average amplitude of the measured \n"
    "            U_lms will be written to the output file. Default 0\n "
    "            Note: This option cannot be selected simultaniusly with \n"
    "           \"-align_axes\".\n"
    "\n"
    " -align_axes\n"
    "            [optional] Takes two integers, \"first index\" and \"second index\" \n"
    "            to specify the orientation of the material frame with respect to \n"
    "            which the RSH analysis will be done. During the analysis, the origin\n"
    "            of the reference frame will be at the Membrane centre of mass and \n"
    "            the frame of reference will be rotated so that the Z axis goes \n"
    "            through the first node (first index) and the second node (second \n"
    "            index) lies on the Z-Y plane.\n"
    "            Note: This option cannot be selected simultaniusly with \"-angavg\".\n"
    " \n"
    " -ext       The Ulm amplitudes will be saved in a file with the same name as the\n"
    "            pdb file provided in the \"pdbpath\" + what you put in after this \n"
    "            flag.\n"
    "            Example (Note: irelevant flags are omited): \n"
    "                ./VCM -pdbpath mydirectory/mypdb.pdb -ext _myUlms.myextension \n"
    "            Output file: \n"
    "                mydirectory/mypdb_myUlms.myextension\n"
    " \n"
    " -framelimits\n"
    "            [optional] Takes the start and end frame of the pdb. The analysis \n"
    "            will be performed only for these frames and all pdb frames in \n"
    "            between. A 0 value for the first or last frame will be interpreted \n"
    "            as the biginning frame and the las frame respectfully. Default is 0\n"
    "            for both limits.\n"
    "            Example (Note: irelevant flags are omited): \n"
    "            1)  ./VCM -framelimits 0 0\n"
    "            2)  ./VCM -framelimits 3 0\n"
    "            3)  ./VCM -framelimits 0 293\n"
    "            4)  ./VCM -framelimits 11 73\n"
    "            interpretation:\n"
    "            1)  All frames.\n"
    "            2)  From frame 3 to the end.\n"
    "            3)  From the beginning to upto and including frame 293.\n"
    "            4)  From and including frame 11 upto and including frame 73.\n"
    " \n"
    <<endl;
}

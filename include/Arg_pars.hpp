//
//  OpenMM_structs.h
//  Membrae
//
//  Created by Ali Farnudi on 10/06/2019.
//  Copyright © 2019 Ali Farnudi. All rights reserved.
//

#ifndef Arg_pars_hpp
#define Arg_pars_hpp

#include <string>
#include <vector>
#include <map>
#include "OpenMM_structs.h"

struct ArgStruct_Analysis
{
    string analysis = "0";
    int analysis_averaging_option = 0;
    int num_ang_avg= 1;
    int z_node=-1;
    int zy_node=-1;
    std::string analysis_filename;
    std::string fileExtension;
    std::vector<std::string> membane_labels;
    std::vector<int>         membane_label_indecies;
    int ell_max =20;
    int q_max   =70;
    std::string extension;
    std::vector<std::string> output_filename;
    int framelimits_beg=0;
    int framelimits_end=0;
    int num_atoms_per_frame=0;
    std::vector<std::string> Mesh_files;
//    std::vector<bool> MeshMinimisation;
    bool MeshMinimisation  = false;
    int  FrameMinimisation = -1;
    
    std::string meshpathinput = "None";
    
    bool undulateMesh = false;
    std::vector<double> SeedNumLmaxUmax;
    
    std::string reference_radius;
    bool dontSquareAmps=true;
    bool use_bin_xyz = false;
};

struct ArgStruct_VCM
{
    std::string  configfilename;
    std::string  resumePath;
    PlatformInfo platforminfo;
    bool platforminput = false;
    bool platformPluginInput = false;
    bool use_voronoi = false;
    bool write_at_end = false;
};

ArgStruct_Analysis  cxxparser_analysis(int argc, char **argv);
ArgStruct_VCM       cxxparser_vcm(int argc, char **argv);


#endif /* OpenMM_structs_h */


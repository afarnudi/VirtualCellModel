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

struct ArgStruct_Analysis
{
    int analysis_dim = 0;
    int analysis_averaging_option = 0;
    int num_ang_avg= 1;
    int z_node=-1;
    int zy_node=-1;
    std::string analysis_filename;
    std::vector<std::string> membane_labels;
    int ell_max =20;
    std::string extension;
    std::vector<std::string> output_filename;
    int framelimits_beg=0;
    int framelimits_end=0;
    int num_atoms_per_frame=0;
    std::vector<std::string> Mesh_files;
};


ArgStruct_Analysis  cxxparser_analysis(int argc, char **argv);
std::string         cxxparser_vcm(int argc, char **argv);


#endif /* OpenMM_structs_h */


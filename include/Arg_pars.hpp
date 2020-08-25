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

using namespace std;

struct ArgStruct
{
    bool analysis_mode=false;
    int analysis_dim = 0;
    int analysis_averaging_option = 0;
    int num_ang_avg= 1;
    int z_node=-1;
    int zy_node=-1;
    string analysis_filename;
    int ell_max =20;
    string output_filename="_ulmt_cpp.txt";
    int framelimits_beg=0;
    int framelimits_end=0;
};

ArgStruct arg_parse_old(int argc, char **argv);
ArgStruct cxxparser(int argc, char **argv);


#endif /* OpenMM_structs_h */


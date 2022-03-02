//
//  General_constants.h
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 26/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef General_constants_h
#define General_constants_h
#include <string>
#include <vector>
#include <map>

#include "ConfigfileStructs.hpp"

struct PotentialModelIndex{
    std::map<std::string, int> Model;
    
    PotentialModelIndex(){
        Model["None"] = 0;
        Model["KremerGrest"] = 1;
        Model["Harmonic"] = 2;
        Model["Kelvin-Voigt"] = 3;
        Model["Maxwell"] = 4;
        Model["HarmonicX4"] = 5;
        Model["Contractile"] = 6;
        Model["Harmonic_minmax"] = 7;
        Model["hill"] = 8;
        Model["KFs"] = 9;
        Model["Gompper"] = 10;
        Model["angleCOS"] = 11;
        Model["Abraham1989"] = 12;
        Model["LocalConstraint"] = 13;
        Model["GlobalConstraint"] = 14;
        
        Model["SmoothTheta4"]=18;
        Model["SmoothEXP"]=19;
        Model["ItzyksonBarycentric"] = 20;
        Model["JulicherVoronoi"] = 21;
        Model["Espiru1987"] = 22;
        Model["Julicher1996"] = 23;
        Model["Itzykson1986"] = 24;
        Model["ExpDihedral"] = 25;
        Model["RealHarmonic"] = 26;
        Model["Dihedral"] = 27;
        Model["Virtual"] = 28;
        
        
    }
};

extern GeneralParameters generalParameters;
extern PotentialModelIndex potentialModelIndex;

namespace GenConst {

    extern int MC_step;
    extern int Mem_fluidity;

    extern std::string force_file_name;
    
//    /**Set the integrator type. This flag is for the OpenMM integrators and will not function if the OpenMM engine is not selected.
//     * Type 0 : Verlet
//     * Type 1 : Brownian, temperature and frictionCoeff need to be set as well.
//     * Type 2 : Langevin, temperature and frictionCoeff need to be set as well.
//     *Default 0*/
//    extern int Integrator_type;

    /**Create  a checkpoint during each saving step. default true*/
    extern bool CreateCheckpoint;
    /**Load the context from OpenMM checkpoint (binary file). Default flase*/
    extern bool Load_from_checkpoint;
    /**Load the context from OpenMM checkpoint (binary file. Default ./Results/Resumes/OpenMM/ */
    extern std::string Checkpoint_path;
    /**name of the OpenMM checkpoint file (binary file)*/
    extern std::string Checkpoint_file_name;
    


    /**Specify if Virtual Sites are used  in the chromatin. Default false*/
    extern bool ChromatinVirtualSites;
    
    
    extern std::vector<double> data_colection_times;

    
//    /**A string containing the selected platform information during runtime.*/
//    extern std::string hardwareReport;

}

#define TRESET   "\e[0m"             /* Reset */
#define TBLINK   "\e[5m"             /* Blink */
#define TBOLD    "\e[1m"             /* Bold */
#define TGRAY    "\e[38;5;249m"      /* Gray */
#define TRED     "\e[38;5;196m"      /* Red */
#define TGREEN   "\e[38;5;118m"      /* Green */
#define TYELLOW  "\e[38;5;226m"      /* Yellow */
#define TBLUE    "\e[38;5;33m"       /* Blue */
#define TPINK    "\e[38;5;200m"      /* Pink */
#define TCYAN    "\e[38;5;123m"      /* Cyan */
#define TWHITE   "\e[38;5;255m"      /* White */
#define TORANGE  "\e[38;5;214m"      /* Orange */
#define TGB      "\e[38;5;50m"       /* Greenish Blue */
#define TDARKGB  "\e[38;5;6m"        /* Dark Greenish Blue */
#define TPURPLE  "\e[38;5;93m"       /* Purple */
#define TLGREEN  "\e[38;5;155m"      /* Light Green */
#define TBORANGE "\e[38;5;166m"      /* Blood Orange */

#define TSUCCESS TGREEN
#define TON      TGREEN
#define TFAILED  TRED
#define TOFF     TRED
#define TFILE    TDARKGB
#define TWARN    TORANGE
#define TWWARN   TBLINK<<TRED

#define TMEM     TBLUE
#define TACT     TCYAN
#define TECM     TDARKGB
#define TCHR     TPINK
#define TOMM     TGB

#define TOCL     TPURPLE
#define TCUD     TLGREEN
#define TCPU     TBORANGE


#endif /* General_constants_h */

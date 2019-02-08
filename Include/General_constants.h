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


namespace GenConst {
    extern int MD_num_of_steps;
    extern int MD_traj_save_step;
    extern double MD_Time_Step;
    extern double MD_T;
    extern double K;
    extern int MD_thrmo_step;
    extern int MC_step;
    extern double Mem_fluidity;
    extern double Lbox;
    extern bool Periodic_condtion_status;
    extern int Num_of_Membranes;
    extern int Num_of_Chromatins;
    extern std::string trajectory_file_name;
    extern bool File_header;
    extern bool Relaxation;
    extern double Buffer_temperature;
    extern double Bussi_tau;
}

//#define MD_num_of_steps  300000//35000// number of MD stps
//#define MD_traj_save_step    2000//The step on which the trajector of the membrane is saved.
//#define MD_Time_Step     0.001 // time length of steps in MD
//#define MD_KT     1.0  // KT the quanta of energy
//#define MD_thrmo_step   100 //
//#define MC_step 1

//#define Mem_fluidity       0.002 //Used in the MC step

//#define Lbox  1000.0///    (size of square periodic box-1)
//#define Periodic_condtion_status 0.0 //status 0.0 = off (The Periodic update will not be executed in the 'Main MD' loop). status = 1.0 = on

#endif /* General_constants_h */

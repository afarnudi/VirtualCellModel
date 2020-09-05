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

namespace GenConst {
    /**Simulation time (in picoseconds). If this parameter is not set in the general config file by the user, or the value is set to zero, it will be calculate during runtime by multiplying the 'step size' by the 'total number of steps'.*/
    extern double Simulation_Time_In_Ps;
    /**Trajectory saving rate in femtoseconds. If this parameter is not set by the user in the general config file, or the value is set to zero, the interval will be calculated by multiplying the MD_traj_save_step by the Step_Size_In_Fs.*/
    extern double Report_Interval_In_Fs;
    /** Integration step size (fs).  Default 0.001*/
    extern double Step_Size_In_Fs;
    /**During every 'MC_step' the MC step will be applied to the membrane mesh. Default 100*/
    extern int MC_step;
    extern int Mem_fluidity;
    /**The simulation uses periodic boundary condition. Default False*/
    extern bool Periodic_box;
    /**he size of the simulation box (cube). If Periodic_box  == False the default value will be set to zero (here -1 will trigger this proccess).*/
    extern double Lbox;
    extern double Simulation_box_length;
    /**Number of membranes in the system followed by the directory of their respective configuration files. All membranes will be configured using the first configuration file if only one is provided. Default 0.*/
    extern int Num_of_Membranes;
    /**Number of chromatins in the system followed by the directory of their respective configuration files. All chromatins will be configured using the first configuration file if only one is provided. Default 0.*/
    extern int Num_of_Chromatins;
    /**Number of actins in the system followed by the directory of their respective configuration files. All actins will be configured using the first configuration file if only one is provided. Default 0.*/
    extern int Num_of_Actins;
    /**Number of ECMs in the system followed by the directory of their respective configuration files. All ECMs will be configured using the first configuration file if only one is provided. Default 0.*/
    extern int Num_of_ECMs;
    /**Name of the output file. Please note that the date and time the file is generated will be attached to this name.*/
    extern std::string trajectory_file_name;
    extern std::string force_file_name;
    
    /**Set the integrator type. This flag is for the OpenMM integrators and will not function if the OpenMM engine is not selected.
     * Type 0 : Verlet
     * Type 1 : Brownian, temperature and frictionCoeff need to be set as well.
     * Type 2 : Langevin, temperature and frictionCoeff need to be set as well.
     *Default 0*/
    extern int Integrator_type;
    /**Set the friction coefficient which couples the system to the heat bath (in inverse picoseconds). Default 5*/
    extern double frictionInPs;
    /**Set the temperature of the heat bath (in Kelvin). Default 300*/
    extern double temperature;
    /**Specifies if an interaction map is provided by the user or not. deafult false.*/
    extern bool Interaction_map;
    /**Path to the interaction map (including the file name). if set to zero no class instances will interact with one another. Set to 1 and provide a path to the "interaction_map.txt" file.  Default = "interaction_map.txt".*/
    extern std::string Interaction_map_file_name;
    /**Set the excluded volume interaction for nodes of class instances. 0 for no repulsion and 1 for excluded volume interaction. Default 0.*/
    extern bool Excluded_volume_interaction;
    /**Set Membrane class label. An index will be assigned during runtime. Default mem*/
    extern std::string Membrane_label;
    /**Set Actin class label. An index will be assigned during runtime. Default act*/
    extern std::string Actin_label;
    /**Set Chromatin class label. An index will be assigned during runtime. Default chr*/
    extern std::string Chromatin_label;
    /**Set ECM class label. An index will be assigned during runtime. Default ecm*/
    extern std::string ECM_label;
    /**Create  a checkpoint during each saving step. default true*/
    extern bool CreateCheckpoint;
    /**Load the context from OpenMM checkpoint (binary file). Default flase*/
    extern bool Load_from_checkpoint;
    /**Load the context from OpenMM checkpoint (binary file. Default ./Results/Resumes/OpenMM/*/
    extern std::string Checkpoint_path;
    /**name of the OpenMM checkpoint file (binary file)*/
    extern std::string Checkpoint_file_name;
    
    extern bool write_bonds_to_PDB;
    /**Collect energy parameters for the potentials (expensive) during each Report_Interval_In_Fs time point. Default true*/
    extern bool   WantEnergy;
    /**Collect forces (cheap) during each Report_Interval_In_Fs time point. Default true*/
    extern bool   WantForce;
    /**Writes velocities and forces (cheap) of particles during each Report_Interval_In_Fs time point to the disk. Default false*/
    extern bool   WriteVelocitiesandForces;
    /**Make the velocity of the centre of mass (COM) zero by subtracting the COM velocity from all the particles' velocity after every CMMotionRemoverStep step. Default false*/
    extern bool CMMotionRemover;

    /**Specify if Virtual Sites are used  in the chromatin. Default false*/
    extern bool ChromatinVirtualSites;
    /**The number of steps where the centre of mass velocity is set to zero using the CMMotionRemover. Default 10*/
    extern int CMMotionRemoverStep;
    
    //non config file parameters
//    extern std::vector<std::vector<std::vector<double> > > data;
    extern std::vector<double> data_colection_times;
    
    //1/0 true/false if you want/don't want to calculate and write the voronoi area associated with Membrane nodes in the properties output file. Default 0
    extern bool Wantvoronoi;

    //When in test mode, most of the console prints (std::cout) of the programme will be turned off for a better viewing of the test report. Deafault 0 (off)
    extern bool Testmode;
    
    /**The pressure acting on the system (in bar) through OpenMM's MonteCarloBarostat. Default 0*/
    extern double MCBarostatPressure;
    /**The temperature at which OpenMM's MonteCarloBarostat will think the system is being maintained (in Kelvin). Default Thermostat temperature*/
    extern double MCBarostatTemperature;
    /**The frequency at which the MonteCarloBarostat pressure changes should be attempted (in time steps). If zero MCBarostat will be disabled. Default 0.*/
    extern int MCBarostatFrequency;
    
    /**The temperature at which OpenMM's MonteCarloBarostat will think the system is being maintained (in Kelvin). Default Thermostat temperature*/
    extern std::vector<std::vector<double> > Lboxdims;
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
#define TWARN    TYELLOW
#define TWWARN   TBLINK<<TRED

#define TMEM     TBLUE
#define TACT     TCYAN
#define TECM     TORANGE
#define TCHR     TPINK
#define TOMM     TGB

#define TOCL     TPURPLE
#define TCUD     TLGREEN
#define TCPU     TBORANGE


#endif /* General_constants_h */

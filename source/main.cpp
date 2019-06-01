/** @file doxygen_example.cpp
 @author Lastname:Firstname:A00123456:cscxxxxx
 @version Revision 1.1
 @brief Illustrates doxygen-style comments for documenting a C++
 program file and the functions in that file.
 @details If you want to add any further detailed description of
 what is in the file, then place it here (after the first statement)
 and it will appear in the detailed description section of the HTML
 output description for the file.
 @date Monday, September 19, 2011
 */

/// \file

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

#include "Membrane.h"
#include "Chromatin.h"
#include "Actin.h"
#include "ECM.h"

#include "General_functions.hpp"
#include "write_functions.hpp"
#include "interaction.hpp"
#include "maps.hpp"
#include "Global_functions.hpp"




namespace GenConst {
    int MD_num_of_steps;
    int MD_traj_save_step;
    double MD_Time_Step;
    double MD_T;
    double K;
    int MD_thrmo_step;
    int MC_step;
    double Mem_fluidity;
    double Lbox;
    bool Periodic_condtion_status;
    int Num_of_Membranes;
    int Num_of_Chromatins;
    int Num_of_Actins;
    int Num_of_ECMs;
    string trajectory_file_name;
    bool File_header;
    double Buffer_temperature;
    double Bussi_tau;
    double Actin_Membrane_Bond_Coefficient;
}


// -----------------------------------------------------------------------------
//                                 MOCK MD CODE
// -----------------------------------------------------------------------------
// The code starting here and through main() below is meant to represent in
// simplified form some pre-existing molecular dynamics code, which defines its
// own data structures for force fields, the atoms in this simulation, and the
// simulation parameters, and takes care of recording the trajectory. All this
// has nothing to do with OpenMM; the OpenMM-dependent code comes later and is
// clearly marked below.
// -----------------------------------------------------------------------------

//                   MODELING AND SIMULATION PARAMETERS
const bool   UseConstraints      = false;   // Should we constrain C-H bonds?
const double StepSizeInFs        = 2;       // integration step size (fs)
const double ReportIntervalInFs  = 10;      // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 100;     // total simulation time (ps)
static const bool   WantEnergy   = true;

//                            FORCE FIELD DATA
// For this example we're using a tiny subset of the Amber99 force field.
// We want to keep the data in the original unit system to avoid conversion
// bugs; this requires conversion on the way in and out of OpenMM.

// Amber reduces nonbonded forces between 1-4 bonded atoms.
//const double Coulomb14Scale      = 0.5;
//const double LennardJones14Scale = 0.5;

struct AtomType {
    double mass, charge, vdwRadiusInAngstroms, vdwEnergyInKcal;
} atomType[] = {
    // mass,     charge, vdwRadius, vdwEnergy
    /*1 C*/12.011, 0.0,   1.9080,    0.1094,
    /*1 C*/12.011, 0.0,   1.9080,    0.1094};

const int H = 0, C = 1;

struct BondType {
    double nominalLengthInAngstroms, stiffnessInKcalPerAngstrom2;
    bool   canConstrain;
} bondType[] = {
    // nominalLength  stiffness   canConstrain;
    /*0 CC*/1.526,       310.,       false,
    /*1 CH*/1.09 ,       340.,       true};

const int CC = 0, CH = 1;

//                                MOLECULE DATA
/** Now describe an ethane molecule by listing its atoms, bonds, angles, and
 * torsions. We'll provide a default configuration which centers the
 * at (0,0,0) with the C-C bond along the x axis.
 *
 * Use this as an "end of list" marker so that we do not have to count; let the
 * computer do that!
 */
const int EndOfList=-1;

struct MyAtomInfo
{
    int type;
    const char* pdb;
    double initPosInAng[3];
    double posInAng[3];
}atoms[] ={
    // type        pdb           initPosInAng            posInAng
    /*0*/C,       " C1 ",       { 0,   0,   0   },     {0,0,0},
    /*1*/C,       " C2 ",       { 0,   0,   1   },     {0,0,0},
    /*2*/C,       " C3 ",       { 0.7, 0,   0.5 },     {0,0,0},
    /*3*/C,       " C4 ",       {-0.5, 0.5, 0.5},     {0,0,0},
    EndOfList};

static struct {int type; int atoms[2];} bonds[] =
{
    // type  atoms[0] atoms[1]
    CC,     0,      1,
    CC,     0,      2,
    CC,     0,      3,
    CC,     1,      2,
    CC,     1,      3,
    EndOfList};

//                               PDB FILE WRITER
// Given state data, output a single frame (pdb "model") of the trajectory.
static void
myWritePDBFrame(int frameNum, double timeInPs, double energyInKcal,
                const MyAtomInfo atoms[])
{
    // Write out in PDB format -- printf is so much more compact than formatted cout.
    printf("MODEL     %d\n", frameNum);
    printf("REMARK 250 time=%.3f ps; energy=%.3f kcal/mole\n",
           timeInPs, energyInKcal);
    for (int n=0; atoms[n].type != EndOfList; ++n)
        printf("ATOM  %5d %4s ETH     1    %8.3f%8.3f%8.3f  1.00  0.00\n",
               n+1, atoms[n].pdb,
               atoms[n].posInAng[0], atoms[n].posInAng[1], atoms[n].posInAng[2]);
    printf("ENDMDL\n");
}



// -----------------------------------------------------------------------------
//                           INTERFACE TO OpenMM
// -----------------------------------------------------------------------------
// These four functions and an opaque structure are used to interface our main
// program with OpenMM without the main program having any direct interaction
// with the OpenMM API. This is a clean approach for interfacing with any MD
// code, although the details of the interface routines will differ. This is
// still just "locally written" code and is not required by OpenMM.
struct MyOpenMMData;
/** This function and an opaque structure are used to interface our main
 * programme with OpenMM without the main programme having any direct interaction
 * with the OpenMM API.
 * This function initilises the openmm system + contex + forces.
 */
static MyOpenMMData* myInitializeOpenMM(const MyAtomInfo atoms[],
                                        double stepSizeInFs,
                                        std::string& platformName);
static void          myStepWithOpenMM(MyOpenMMData*, int numSteps);
static void          myGetOpenMMState(MyOpenMMData*, bool wantEnergy,
                                      double& time, double& energy,
                                      MyAtomInfo atoms[]);
static void          myTerminateOpenMM(MyOpenMMData*);




int main(int argc, char **argv)
{
    //time
    
    clock_t tStart = clock();//Time the programme
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_time_%H_%M",now);
    
    string general_file_name="general-config.txt";
    cout<<"Please enter the path (relative to the binary file) + name of the config file:\n";
    cin>>general_file_name;
    vector<string> membrane_config_list;
    vector<string> chromatin_config_list;
    vector<string> actin_config_list;
    vector<string> ecm_config_list;
    
    read_general_parameters(general_file_name, membrane_config_list, chromatin_config_list, actin_config_list, ecm_config_list);
    
    ofstream Trajectory;
    string traj_file_name="Results/"+GenConst::trajectory_file_name+buffer+".xyz";
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    
    vector<Membrane> Membranes;
    vector<Chromatin> Chromatins;
    vector<Actin> Actins;
    vector<ECM> ECMs;
    
    bool Include_Membrane  = false;
    bool Include_Chromatin = false;
    bool Include_Actin     = false;
    bool Include_ECM       = false;
    
    if (GenConst::Num_of_Membranes!=0) {
        Include_Membrane = true;
        
        Membranes.resize(GenConst::Num_of_Membranes);
        for (int i=0; i<GenConst::Num_of_Membranes; i++) {
            string label="Membrane_"+to_string(i);
            Membranes[i].set_label(label);
            Membranes[i].set_file_time(buffer);
            Membranes[i].set_index(i);
            Membranes[i].import_config(membrane_config_list[i]);
            Membranes[i].generate_report();
        }
    }
    
    if (GenConst::Num_of_Chromatins!=0) {
        Include_Chromatin = true;
        Chromatins.resize(GenConst::Num_of_Chromatins);
        for (int i=0; i<GenConst::Num_of_Chromatins; i++) {
            Chromatins[i].set_file_time(buffer);
            Chromatins[i].set_index(i);
            if (GenConst::Num_of_Membranes == GenConst::Num_of_Chromatins) {
                ///put a flag for chromatin inside membrane
                Chromatins[i].import_config(chromatin_config_list[i], Membranes[i].return_min_radius_after_relaxation());
            } else {
                Chromatins[i].import_config(chromatin_config_list[i]);
            }
            
            Chromatins[i].generate_report();
        }
    }
    
    if (GenConst::Num_of_Actins!=0) {
        Include_Actin = true;
        Actins.resize(GenConst::Num_of_Actins);
        for (int i=0; i<GenConst::Num_of_Actins; i++) {
            Actins[i].set_file_time(buffer);
            Actins[i].set_index(i);
            Actins[i].import_config(actin_config_list[i]);
            //            Actins[i].generate_report();
        }
    }
    
    if (GenConst::Num_of_ECMs!=0){
        Include_ECM=true;
        ECMs.resize(GenConst::Num_of_ECMs);
        for (int i=0; i<GenConst::Num_of_ECMs; i++) {
            ECMs[i].set_file_time(buffer);
            ECMs[i].set_index(i);
            ECMs[i].import_config(ecm_config_list[i]);
        }
        
    }
    
    
    if (Include_Membrane && Include_ECM) {
        for (int i=0; i<GenConst::Num_of_Membranes; i++) {
            for (int j=0; j<GenConst::Num_of_ECMs; j++) {
                Membrane_ECM_neighbour_finder(ECMs[j], Membranes[i]);
            }
        }
    }
    
    
    
    int num_of_elements=0;
    if (Include_Membrane) {
        for (int i=0; i<Membranes.size(); i++) {
            num_of_elements+=Membranes[i].return_num_of_nodes();
        }
    }
    if (Include_Chromatin) {
        for (int i=0; i<Chromatins.size(); i++) {
            num_of_elements+=Chromatins[i].return_num_of_nodes();
        }
    }
    if (Include_Actin) {
        for (int i=0; i<Actins.size(); i++) {
            num_of_elements+=Actins[i].return_num_of_nodes();
        }
    }
    if (Include_ECM) {
        for (int i=0; i<ECMs.size(); i++) {
            num_of_elements+=ECMs[i].return_num_of_nodes();
        }
    }
    
    if (Include_Membrane){
        if (Include_Actin){
            for (int i=0; i<GenConst::Num_of_Actins; i++) {
                for (int j=0; j<Membranes.size(); j++) {
                    Actin_Membrane_shared_Node_Identifier(Actins[i], Membranes[j] , j);
                    if (Membranes[j].return_relax_with_actin_flag()) {
                        Membranes[j].Relax_1();
                    }
                }
                
            } //for (int i=0; i<GenConst::Num_of_Actins; i++)
        } //if (Include_Actin)
        else{
            for (int i=0; i<Membranes.size(); i++){
                Membranes[i].Relax_1();
            }// End of for (int i=0; i<Membranes.size(); i++)
        }//end else
    } // End of if (Include_Membrane)
    
    bool openmm=false;
    
    if (openmm) {
        cout<<"\nBeginnig the OpenMM section:\n";
        // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
        // usage and runtime errors are caught and reported.
        try {
            std::string   platformName;
            
            // Set up OpenMM data structures; returns OpenMM Platform name.
            MyOpenMMData* omm = myInitializeOpenMM(atoms, StepSizeInFs, platformName);
            
            // Run the simulation:
            //  (1) Write the first line of the PDB file and the initial configuration.
            //  (2) Run silently entirely within OpenMM between reporting intervals.
            //  (3) Write a PDB frame when the time comes.
            printf("REMARK  Using OpenMM platform %s\n", platformName.c_str());
            
            const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);
            for (int frame=1; ; ++frame) {
                double time, energy;
                myGetOpenMMState(omm, WantEnergy, time, energy, atoms);
                myWritePDBFrame(frame, time, energy, atoms);
                
                if (time >= SimulationTimeInPs)
                    break;
                
                myStepWithOpenMM(omm, NumSilentSteps);
            }
            
            // Clean up OpenMM data structures.
            myTerminateOpenMM(omm);
            return 0; // Normal return from main.
        }
        
        // Catch and report usage and runtime errors detected by OpenMM and fail.
        catch(const std::exception& e) {
            printf("EXCEPTION: %s\n", e.what());
            return 1;
        }
    }
    
    int progress=0;
    cout<<"\nBeginnig the MD\nProgress:\n";
    for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_steps ; MD_Step++){
//        cout<<Membranes[0].return_node_position(0, 0);
        
        //Thermostat step first step
        if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {
            if (Include_Membrane) {
                for (int i=0; i<Membranes.size(); i++) {
                    Membranes[i].Thermostat_Bussi(GenConst::Buffer_temperature);
                }
            }
            if (Include_Actin) {
                for (int i=0; i<Actins.size(); i++) {
                    //                    Actins[i].Thermostat_Bussi(GenConst::MD_T);
                }
            }
            if (Include_Chromatin) {
                for (int i=0; i<Chromatins.size(); i++) {
                    Chromatins[i].Thermostat_Bussi(GenConst::MD_T*0.01);
                }
            }
            if (Include_ECM) {
                for (int i=0; i<ECMs.size(); i++) {
                    //                    ECMs[i].Thermostat_Bussi(GenConst::MD_T);
                }
            }
        }
        
        
        //Velocity Verlet first step
        if (Include_Membrane)
        {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].MD_Evolution_beginning(GenConst::MD_Time_Step);
            }
        }
        if (Include_Chromatin)
        {
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatins[i].MD_Evolution_beginning(GenConst::MD_Time_Step);
            }
        }
        if (Include_Actin)
        {
            for (int i=0; i<Actins.size(); i++) {
                Actins[i].MD_Evolution_beginning(GenConst::MD_Time_Step);
            }
        }
        if (Include_ECM)
        {
            for (int i=0; i<ECMs.size(); i++) {
                ECMs[i].MD_Evolution_beginning(GenConst::MD_Time_Step);
            }
        }
        
        
        
        //force implamentation
        if (Include_Membrane)
        {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].Elastic_Force_Calculator(0);
            }
        }
        if (Include_Chromatin)
        {
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatins[i].Force_Calculator_2();
            }
        }
        if (Include_Actin)
        {
            for (int i=0; i<Actins.size(); i++) {
                Actins[i].Elastic_Force_Calculator();
            }
        }
        if (Include_ECM)
        {
            for (int i=0; i<ECMs.size(); i++) {
                //                ECMs[i].Elastic_Force_Calculator();
            }
        }
        
        //Shared Forces
        if (Include_Chromatin && Include_Membrane) {
            if (MD_Step%2000==0) {
                for (int i=0; i<Chromatins.size(); i++) {
                    Chromatin_Membrane_neighbour_finder(Chromatins[i], Membranes[i]);
                    Chromatin_Membrane_hard_sphere(Chromatins[i], Membranes[i]);
                }
            }
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatin_Membrane_hard_sphere(Chromatins[i], Membranes[i]);
            }
            
            
        }
        
        if (Include_Membrane && Include_Actin) {
            for (int i=0; i<Actins.size(); i++) {
                for (int j=0; j<Membranes.size(); j++) {
                    Actin_Membrane_shared_Node_Force_calculator(Actins[i], Membranes[j], j);
                }
            }
        }
        
        if (Include_Membrane && Include_ECM) {
            for (int i=0; i<Membranes.size(); i++) {
                for (int j=0; j<ECMs.size(); j++) {
                    Membrane_ECM_shared_node_force (ECMs[j], Membranes[i]);
                    if (MD_Step%2000==0) {
                        update_ecm_mem_neighbour_list (ECMs[j], Membranes[i]);
                    }
                }
            }
        }
        
        
        
        
        //Velocity Verlet second step
        if (Include_Membrane) {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].MD_Evolution_end(GenConst::MD_Time_Step);
            }
        }
        if (Include_Chromatin) {
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatins[i].MD_Evolution_end(GenConst::MD_Time_Step);
            }
        }
        if (Include_Actin) {
            for (int i=0; i<Actins.size(); i++) {
                Actins[i].MD_Evolution_end(GenConst::MD_Time_Step);
            }
        }
        if (Include_ECM) {
            for (int i=0; i<ECMs.size(); i++) {
                ECMs[i].MD_Evolution_end(GenConst::MD_Time_Step);
            }
        }
        
        
        
        //Thermostat second step
        if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {
            if (Include_Membrane) {
                for (int i=0; i<Membranes.size(); i++) {
                    Membranes[i].Thermostat_Bussi(GenConst::Buffer_temperature);
                }
            }
            if (Include_Actin) {
                for (int i=0; i<Actins.size(); i++) {
                    //                    Actins[i].Thermostat_Bussi(GenConst::MD_T);
                }
            }
            if (Include_Chromatin) {
                for (int i=0; i<Chromatins.size(); i++) {
                    Chromatins[i].Thermostat_Bussi(GenConst::MD_T*0.01);
                }
            }
            if (Include_ECM) {
                for (int i=0; i<ECMs.size(); i++) {
                    //                    ECMs[i].Thermostat_Bussi(GenConst::MD_T);
                }
            }
        }
        
        //saving Results
        if (MD_Step%GenConst::MD_traj_save_step == 0)
        {
            Trajectory << num_of_elements<<endl;
            Trajectory << " nodes  "<<endl;
            
            
            if (Include_Membrane) {
                for (int i=0; i<Membranes.size(); i++) {
//                    string label="Membrane_"+to_string(i);
                    Membranes[i].write_traj(traj_file_name);
                    Membranes[i].export_for_resume(MD_Step);
                }
            }
            
            if (Include_Chromatin) {
                for (int i=0; i<Chromatins.size(); i++) {
                    string label="Chromatin_"+to_string(i);
                    Chromatins[i].write_traj(traj_file_name, label);
                    Chromatins[i].export_for_resume(MD_Step);
                }
            }
            if (Include_Actin) {
                for (int i=0; i<Actins.size(); i++) {
                    string label="Actin_"+to_string(i);
                    Actins[i].write_traj(traj_file_name, label);
                    //                    Actins[i].export_for_resume(MD_Step);
                }
            }
            if (Include_ECM) {
                for (int i=0; i<ECMs.size(); i++) {
                    string label="ECM_"+to_string(i);
                    ECMs[i].write_traj(traj_file_name, label);
                    //                    Actins[i].export_for_resume(MD_Step);
                }
            }
        }// End of if (MD_Step%100==0)
        
        
        if (int(100*MD_Step/GenConst::MD_num_of_steps)>progress){
            cout<<"[ "<<progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;
            progress+=5;
        }
        
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    cout<<"[ 100% ]\t step: "<<GenConst::MD_num_of_steps<<"\n";
    cout<<"\nDone!"<<endl;
    printf("Time taken: %.2f Minutes\n", (double)((clock() - tStart)/CLOCKS_PER_SEC)/60.0);
    return 0;
}

// -----------------------------------------------------------------------------
//                           OpenMM-USING CODE
// -----------------------------------------------------------------------------
// The OpenMM API is visible only at this point and below. Normally this would
// be in a separate compilation module; we're including it here for simplicity.
// -----------------------------------------------------------------------------

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
#pragma warning(disable:4996)   // sprintf is unsafe
#endif

#include "OpenMM.h"
using OpenMM::Vec3; // so we can just say "Vec3" below

// This is our opaque "handle" class containing all the OpenMM objects that
// must persist from call to call during a simulation. The main program gets
// a pointer to one of these but sees it as essentially a void* since it
// doesn't know the definition of this class.
struct MyOpenMMData {
    MyOpenMMData() : system(0), context(0), integrator(0) {}
    ~MyOpenMMData() {delete context; delete integrator; delete system;}
    OpenMM::System*         system;
    OpenMM::Integrator*     integrator;
    OpenMM::Context*  context;
};

// -----------------------------------------------------------------------------
//                      INITIALIZE OpenMM DATA STRUCTURES
// -----------------------------------------------------------------------------
// We take these actions here:
// (1) Load any available OpenMM plugins, e.g. Cuda and Brook.
// (2) Allocate a MyOpenMMData structure to hang on to OpenMM data structures
//     in a manner which is opaque to the caller.
// (3) Fill the OpenMM::System with the force field parameters we want to
//     use and the particular set of atoms to be simulated.
// (4) Create an Integrator and a Context associating the Integrator with
//     the System.
// (5) Select the OpenMM platform to be used.
// (6) Return the MyOpenMMData struct and the name of the Platform in use.
//
// Note that this function must understand the calling MD code's molecule and
// force field data structures so will need to be customized for each MD code.
static MyOpenMMData*
myInitializeOpenMM( const MyAtomInfo    atoms[],
                   double              stepSizeInFs,
                   std::string&        platformName)
{
    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory
    (OpenMM::Platform::getDefaultPluginsDirectory());
    
    // Allocate space to hold OpenMM objects while we're using them.
    MyOpenMMData* omm = new MyOpenMMData();
    
    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the System owns
    // the force objects and will take care of deleting them; don't do it yourself!
    OpenMM::System&                 system      = *(omm->system = new OpenMM::System());
//    OpenMM::NonbondedForce&         nonbond     = *new OpenMM::NonbondedForce();
    OpenMM::HarmonicBondForce&      bondStretch = *new OpenMM::HarmonicBondForce();
//    OpenMM::HarmonicAngleForce&     bondBend    = *new OpenMM::HarmonicAngleForce();
//    OpenMM::PeriodicTorsionForce&   bondTorsion = *new OpenMM::PeriodicTorsionForce();
//    system.addForce(&nonbond);
    system.addForce(&bondStretch);
//    system.addForce(&bondBend);
//    system.addForce(&bondTorsion);
    
    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    std::vector<Vec3> initialPosInNm;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atype.mass);
//        nonbond.addParticle(atype.charge,
//                            atype.vdwRadiusInAngstroms * OpenMM::NmPerAngstrom
//                            * OpenMM::SigmaPerVdwRadius,
//                            atype.vdwEnergyInKcal      * OpenMM::KJPerKcal);
        // Convert the initial position to nm and append to the array.
        const Vec3 posInNm(atoms[n].initPosInAng[0] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[1] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[2] * OpenMM::NmPerAngstrom);
        initialPosInNm.push_back(posInNm);
    }
    
    // Process the bonds:
    //  (1) If we're using constraints, tell System about constrainable bonds;
    //      otherwise, tell HarmonicBondForce the bond stretch parameters
    //      (tricky units!).
    //  (2) Create a list of bonds for generating nonbond exclusions.
    std::vector< std::pair<int,int> > bondPairs;
    for (int i=0; bonds[i].type != EndOfList; ++i) {
        const int*      atom = bonds[i].atoms;
        const BondType& bond = bondType[bonds[i].type];
        
        if (UseConstraints && bond.canConstrain) {
            system.addConstraint(atom[0], atom[1],
                                 bond.nominalLengthInAngstroms * OpenMM::NmPerAngstrom);
        } else {
            // Note factor of 2 for stiffness below because Amber specifies the constant
            // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
            // it as used in the force term kx, with energy kx^2/2.
            bondStretch.addBond(atom[0], atom[1],
                                bond.nominalLengthInAngstroms
                                * OpenMM::NmPerAngstrom,
                                bond.stiffnessInKcalPerAngstrom2
                                * 2 * OpenMM::KJPerKcal
                                * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
        }
        
        bondPairs.push_back(std::make_pair(atom[0], atom[1]));
    }
    // Exclude 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms.
//    nonbond.createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);
    
    // Create the 1-2-3 bond angle harmonic terms.
//    for (int i=0; angles[i].type != EndOfList; ++i) {
//        const int*       atom  = angles[i].atoms;
//        const AngleType& angle = angleType[angles[i].type];
//
//        // See note under bond stretch above regarding the factor of 2 here.
//        bondBend.addAngle(atom[0],atom[1],atom[2],
//                          angle.nominalAngleInDegrees     * OpenMM::RadiansPerDegree,
//                          angle.stiffnessInKcalPerRadian2 * 2 * OpenMM::KJPerKcal);
//    }
    
    // Create the 1-2-3-4 bond torsion (dihedral) terms.
//    for (int i=0; torsions[i].type != EndOfList; ++i) {
//        const int*         atom = torsions[i].atoms;
//        const TorsionType& torsion = torsionType[torsions[i].type];
//        bondTorsion.addTorsion(atom[0],atom[1],atom[2],atom[3],
//                               torsion.periodicity,
//                               torsion.phaseInDegrees  * OpenMM::RadiansPerDegree,
//                               torsion.amplitudeInKcal * OpenMM::KJPerKcal);
//    }
    
    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero.
    omm->integrator = new OpenMM::VerletIntegrator(StepSizeInFs * OpenMM::PsPerFs);
    omm->context    = new OpenMM::Context(*omm->system, *omm->integrator);
    omm->context->setPositions(initialPosInNm);
    
    platformName = omm->context->getPlatform().getName();
    return omm;
}

// -----------------------------------------------------------------------------
//                     COPY STATE BACK TO CPU FROM OPENMM
// -----------------------------------------------------------------------------
static void
myGetOpenMMState(MyOpenMMData* omm, bool wantEnergy,
                 double& timeInPs, double& energyInKcal,
                 MyAtomInfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    if (wantEnergy) {
        infoMask += OpenMM::State::Velocities; // for kinetic energy (cheap)
        infoMask += OpenMM::State::Energy;     // for pot. energy (expensive)
    }
    // Forces are also available (and cheap).
    
    const OpenMM::State state = omm->context->getState(infoMask);
    timeInPs = state.getTime(); // OpenMM time is in ps already
    
    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    for (int i=0; i < (int)positionsInNm.size(); ++i)
        for (int j=0; j < 3; ++j)
            atoms[i].posInAng[j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;
    
    // If energy has been requested, obtain it and convert from kJ to kcal.
    energyInKcal = 0;
    if (wantEnergy)
        energyInKcal = (state.getPotentialEnergy() + state.getKineticEnergy())
        * OpenMM::KcalPerKJ;
}

// -----------------------------------------------------------------------------
//                     TAKE MULTIPLE STEPS USING OpenMM
// -----------------------------------------------------------------------------
static void
myStepWithOpenMM(MyOpenMMData* omm, int numSteps) {
    omm->integrator->step(numSteps);
}

// -----------------------------------------------------------------------------
//                     DEALLOCATE OpenMM OBJECTS
// -----------------------------------------------------------------------------
static void
myTerminateOpenMM(MyOpenMMData* omm) {
    delete omm;
}


#include <iostream>
#include <gtest/gtest.h>
#include <string>

#include "maps.hpp"
#include "Membrane.h"
#include "General_functions.hpp"
#include "Arg_pars.hpp"


namespace GenConst {
int    MD_num_of_steps;
double Simulation_Time_In_Ps;
int    MD_traj_save_step;
double Report_Interval_In_Fs;
double Step_Size_In_Fs;
double MD_T;
double K;
int    MD_thrmo_step;
int    MC_step;
int    Mem_fluidity;
bool   Periodic_box;
double Lbox;
bool   Periodic_condtion_status;
int    Num_of_Membranes;
int    Num_of_Chromatins;
int    Num_of_Actins;
int    Num_of_ECMs;
int    Num_of_pointparticles;
string trajectory_file_name;;
string force_file_name;;
double Buffer_temperature;
double Bussi_tau;
double Actin_Membrane_Bond_Coefficient;
bool   Interaction_map;
string Interaction_map_file_name;
bool   Excluded_volume_interaction;
bool   OpenMM;
double sigma_LJ_12_6;
double epsilon_LJ_12_6;
string Membrane_label;
string Actin_label;
string Chromatin_label;
string ECM_label;
int    Integrator_type;
double frictionInPs;
double temperature;
bool   CreateCheckpoint;
bool   Load_from_checkpoint;
string Checkpoint_path;
string Checkpoint_file_name;
bool   ChromatinVirtualSites;


bool   write_bonds_to_PDB;
bool   WantEnergy;
bool   WantForce;
bool   WriteVelocitiesandForces;
bool   CMMotionRemover;
int    CMMotionRemoverStep;
bool   Wantvoronoi;
bool   Testmode;


double MCBarostatPressure;
double MCBarostatTemperature;
int    MCBarostatFrequency;


//    std::vector<std::vector<std::vector<double> > > data;
std::vector<double> data_colection_times;
std::vector<std::vector<double> > Lboxdims;
}

using namespace std;
struct MemClassTest : testing::Test{
    //setup
    Membrane* Mem = new Membrane;
    
    MemClassTest() {
        
        string general_file_name="test_conf.txt";
        vector<string> membrane_config_list;
        //The following configfiles are just needed for the structure of the "read_general_parameters" function
        vector<string> chromatin_config_list;
        vector<string> actin_config_list;
        vector<string> ecm_config_list;
        vector<string> pointparticle_config_list;
        
        read_general_parameters(general_file_name, membrane_config_list, chromatin_config_list, actin_config_list, ecm_config_list, pointparticle_config_list);
        
        
        
        string label=GenConst::Membrane_label+to_string(0);
        Mem->set_label(label);
        Mem->set_index(0);
        Mem->import_config(membrane_config_list[0]);
        
    }
    ~MemClassTest() {
        delete Mem;
    }
    
};


TEST_F( MemClassTest, NumOfNodes){
    ASSERT_EQ(Mem->get_num_of_nodes() , 1002);
}

TEST_F( MemClassTest, NumOfTriangles){
    ASSERT_EQ(Mem->get_num_of_triangle() , 2000);
}

TEST_F( MemClassTest, NumOfBonds){
    ASSERT_EQ(Mem->get_num_of_node_pairs() , 3000);
}

TEST_F( MemClassTest, NumOfTrianglePairs){
    ASSERT_EQ(Mem->get_num_of_triangle_pairs() , 3000);
}

int main (int argc, char **argv){
    testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}

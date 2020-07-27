#include <iostream>
#include <gtest/gtest.h>
#include <string>

#include "maps.hpp"
#include "Membrane.h"
#include "General_functions.hpp"
//#include "Tests.hpp"


namespace GenConst {
    int MD_num_of_steps;
    double Simulation_Time_In_Ps;
    int MD_traj_save_step;
    double Report_Interval_In_Fs;
    double Step_Size_In_Fs;
    double MD_T;
    double K;
    int MD_thrmo_step;
    int MC_step;
    int Mem_fluidity;
    double Lbox;
    bool Periodic_condtion_status;
    int Num_of_Membranes;
    int Num_of_Chromatins;
    int Num_of_Actins;
    int Num_of_ECMs;
    int Num_of_pointparticles;
    std::string trajectory_file_name;;
    double Buffer_temperature;
    double Bussi_tau;
    double Actin_Membrane_Bond_Coefficient;
    bool Interaction_map;
    std::string Interaction_map_file_name;
    bool Excluded_volume_interaction;
    bool OpenMM;
    double sigma_LJ_12_6;
    double epsilon_LJ_12_6;
    std::string Membrane_label;
    std::string Actin_label;
    std::string Chromatin_label;
    std::string ECM_label;
    int Integrator_type;
    double frictionInPs;
    double temperature;
    bool CreateCheckpoint;
    bool Load_from_checkpoint;
    std::string Checkpoint_path;
    std::string Checkpoint_file_name;
    bool ChromatinVirtualSites;


    bool   write_bonds_to_PDB;
    bool   WantEnergy;
    bool   WantForce;
    bool   WriteVelocitiesandForces;
    bool   CMMotionRemover;
    int    CMMotionRemoverStep;
    bool   Wantvoronoi;
    bool   Testmode;



    //    std::vector<std::vector<std::vector<double> > > data;
    std::vector<double> data_colection_times;
}

using namespace std;
struct MemClassTest : protected testing::Test{
    //setup
    
    
    vector<Membrane*> Membranes;
    
    MemClassTest() {
        
        string general_file_name="test_conf.txt";
        vector<string> membrane_config_list;
        //The following configfiles are just needed for the structure of the "read_general_parameters" function
        vector<string> chromatin_config_list;
        vector<string> actin_config_list;
        vector<string> ecm_config_list;
        vector<string> pointparticle_config_list;
        
        read_general_parameters(general_file_name, membrane_config_list, chromatin_config_list, actin_config_list, ecm_config_list, pointparticle_config_list);
        
        vector<vector<int> > interaction_map;
        read_interaction_map(interaction_map);
        
        bool Include_Membrane  = false;
        
        int num_of_atoms=0;
        int num_of_bonds=0;
        int num_of_dihedrals=0;
        
        if (!GenConst::Load_from_checkpoint) {
            if (GenConst::Num_of_Membranes!=0) {
                Include_Membrane = true;
                
                for (int i=0; i<GenConst::Num_of_Membranes; i++) {
                    Membranes.push_back(new Membrane);
                    string label=GenConst::Membrane_label+to_string(i);
                    Membranes[i]->set_label(label);
                    Membranes[i]->set_index(i);
                    Membranes[i]->import_config(membrane_config_list[i]);
                }
            }
            
            
        }
    }
    ~MemClassTest() {
        for (auto pointer : Membranes){
         delete pointer;
       }
       Membranes.clear();
    }
    
};


TEST_F( MemClassTest, Num_of_nodes){
    ASSERT_EQ(Membranes[0]->get_num_of_nodes() , 1002);
}

int main (int argc, char **argv){
    testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}

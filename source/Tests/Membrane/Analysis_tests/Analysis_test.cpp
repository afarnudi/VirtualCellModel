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
struct MemAnalysisTest : protected testing::Test{
    //setup
    
    
    vector<Membrane*> Membranes;
    
    MemAnalysisTest() {
        bool Include_Membrane = false;
        bool analysis_mode=true;
        int analysis_averaging_option = 0;
        int num_ang_avg= 1;
        int z_node=-1;
        int y_node=-1;
        
        int ell_max =20;
        std::string analysis_extension = "_ulmt_cpp.txt";
        //    cout<<"argc "<<argc<<endl;
        
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
        
//        int max_frame = Membranes[0].import_pdb_frames(analysis_filename);
//        int  L=2,M=2;
//        double U=0.00001;
//
//        string temp_pdb_name = analysis_filename;
//        temp_pdb_name.pop_back();
//        temp_pdb_name.pop_back();
//        temp_pdb_name.pop_back();
//        temp_pdb_name.pop_back();
//
//        string lmtrajname ="Results/real_real/ulm_"+to_string(U)+"_"+to_string(L)+"_"+to_string(M)+".pdb";
//        vector<int> Ells;
//        vector<int> Ms;
//        Ells.resize(superposes,0);
//        Ms.resize(superposes,0);
//
//        for (int sup =0 ; sup<superposes; sup++) {
//
//        }
//        for (int i=2; i<max_frame; i++) {
//
//
//            Membranes[0].load_pdb_frame(i, analysis_averaging_option, z_node, y_node);
//            Membranes[0].generate_ulm_mode_real(L, M, U);
//            //                    Membranes[0].generate_ulm_mode(L, M, U);
//            for (int runs=0; runs<num_ang_avg; runs++) {
//                //                Membranes[0].calculate_ulm(ell_max, analysis_averaging_option);
//                //                        Membranes[0].surface_integral_test();
//                Membranes[0].calculate_ulm_radiustest_real(ell_max, analysis_averaging_option);
//                //                        Membranes[0].calculate_ulm_radiustest(ell_max, analysis_averaging_option);
//                //                        Membranes[0].myWritePDBFrame(runs,lmtrajname);
//                //                Membranes[0].calculate_ulm_sub_particles(ell_max, analysis_averaging_option);
//            }
//
//        }
        
        
    }
    ~MemAnalysisTest() {
        for (auto pointer : Membranes){
            delete pointer;
        }
        Membranes.clear();
    }
    
};


TEST_F( MemAnalysisTest, NumOfframes){
    std::string analysis_filename = "10_1002n3frames.pdb";
    ASSERT_EQ(count_pdb_frames(analysis_filename,Membranes[0]->get_num_of_nodes()) , 3);
}


#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>
using namespace std::complex_literals;

TEST_F( MemAnalysisTest, Ylm00Theta0Phi0){
    
    std::complex<double> ylm;
    int ell=0,m=0;
    double theta =0, phi = 0;
    ylm = boost::math::spherical_harmonic(ell,m,theta,phi);
    ASSERT_EQ(ylm, 0.5*sqrt(1/M_PI));
}

TEST_F( MemAnalysisTest, Ylm00ThetaPiPhi2Pi){
    
    std::complex<double> ylm;
    int ell=0,m=0;
    double theta =M_PI, phi = 2*M_PI;
    
    ylm = boost::math::spherical_harmonic(ell,m,theta,phi);
    ASSERT_EQ(ylm, 0.5*sqrt(1/M_PI));
}
TEST_F( MemAnalysisTest, Ylm00ThetaPiD3Phi2PiD3){
    
    std::complex<double> ylm;
    int ell=0,m=0;
    double theta =M_PI/3, phi = 2*M_PI/7;
    ylm = boost::math::spherical_harmonic(ell,m,theta,phi);
    ASSERT_EQ(ylm, 0.5*sqrt(1/M_PI));
}
TEST_F( MemAnalysisTest, Ylm00ThetaRandPhiRand){
    
    std::complex<double> ylm;
    int ell=0,m=0;
    double theta =M_PI*((double) rand() / (RAND_MAX)), phi = 2*M_PI*((double) rand() / (RAND_MAX));
    ylm = boost::math::spherical_harmonic(ell,m,theta,phi);
    ASSERT_EQ(ylm, 0.5*sqrt(1/M_PI));
}

int main (int argc, char **argv){
    testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}

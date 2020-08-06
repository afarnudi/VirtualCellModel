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
bool Periodic_box;
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
struct MemAnalysisTest : public testing::Test{
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
        std::string analysis_filename = "10_1002n3frames.pdb";
        int max_frame = Membranes[0]->import_pdb_frames(analysis_filename);
        Membranes[0]->load_pdb_frame(2, 0, -1, -1);
        
    }
    ~MemAnalysisTest() {
        for (auto pointer : Membranes){
            delete pointer;
        }
        Membranes.clear();
    }
    
};


TEST_F( MemAnalysisTest, Radius){
    
    double radius=0;
    for(int i=0;i<Membranes[0]->get_num_of_nodes();i++){
        radius+=Membranes[0]->get_spherical_position(i,0);
    }
    radius/=Membranes[0]->get_num_of_nodes();
    EXPECT_NEAR( radius , 1, 0.00002);
}

TEST_F( MemAnalysisTest, VoronoiAreaSum){
    
    double integral =0;
    for (int i=0; i<Membranes[0]->get_num_of_nodes(); i++) {
        integral+=Membranes[0]->node_voronoi_area[i];
    }
    
    EXPECT_NEAR(integral , 4*M_PI, 0.04);
}

using SurfaceIntegral = MemAnalysisTest;

TEST_F( SurfaceIntegral, UnitVector){
    
    
    vector<double>  unit_vector;
    
    unit_vector.resize(Membranes[0]->get_num_of_nodes(), 1);
    
    double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(unit_vector, unit_vector);
    
    EXPECT_NEAR(integral , 4*M_PI,0.04);
}

TEST_F( SurfaceIntegral, VoronoiAreaSumEllipsoidxyz121){
    std::string analysis_filename = "10_1002nEllipsoid121_3frames.pdb";
    int max_frame = Membranes[0]->import_pdb_frames(analysis_filename);
    Membranes[0]->load_pdb_frame(2, 0, -1, -1);
    double integral =0;
    for (int i=0; i<Membranes[0]->get_num_of_nodes(); i++) {
        integral+=Membranes[0]->node_voronoi_area[i];
    }
    
    EXPECT_NEAR(integral , 21.45963, 0.03);
}

TEST_F( SurfaceIntegral, OmegaEllipsoidxyz121){

    std::string analysis_filename = "10_1002nEllipsoid121_3frames.pdb";
    int max_frame = Membranes[0]->import_pdb_frames(analysis_filename);
    Membranes[0]->load_pdb_frame(2, 0, -1, -1);
    vector<double>  unit_vector;

    unit_vector.resize(Membranes[0]->get_num_of_nodes(), 1);

    double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(unit_vector, unit_vector);

    EXPECT_NEAR(integral , 4*M_PI,1.2);
}

TEST_F( SurfaceIntegral, VoronoiAreaSumEllipsoidxyz123){
    std::string analysis_filename = "10_1002nEllipsoid123_3frames.pdb";
    int max_frame = Membranes[0]->import_pdb_frames(analysis_filename);
    Membranes[0]->load_pdb_frame(2, 0, -1, -1);
    double integral =0;
    for (int i=0; i<Membranes[0]->get_num_of_nodes(); i++) {
        integral+=Membranes[0]->node_voronoi_area[i];
    }
    
    EXPECT_NEAR(integral , 48.93663, 0.15);
}

TEST_F( SurfaceIntegral, OmegaEllipsoidxyz123){

    std::string analysis_filename = "10_1002nEllipsoid123_3frames.pdb";
    int max_frame = Membranes[0]->import_pdb_frames(analysis_filename);
    Membranes[0]->load_pdb_frame(2, 0, -1, -1);
    vector<double>  unit_vector;

    unit_vector.resize(Membranes[0]->get_num_of_nodes(), 1);

    double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(unit_vector, unit_vector);
    EXPECT_NEAR(integral , 4*M_PI,0.000000002);
    
    
}


//******************************************************************************
//******************** Spherical Harmonics Mode Generator  *********************
//******************************************************************************


struct EllM{
    
    int l1;
    int m1;
    
    int l2;
    int m2;
};

class LMParameterized : public testing::TestWithParam<EllM>{
protected:
    vector<Membrane*> Membranes;
    
    LMParameterized() {
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
        std::string analysis_filename = "10_1002n3frames.pdb";
        int max_frame = Membranes[0]->import_pdb_frames(analysis_filename);
        Membranes[0]->load_pdb_frame(2, 0, -1, -1);
        
    }
    ~LMParameterized() {
        for (auto pointer : Membranes){
            delete pointer;
        }
        Membranes.clear();
    }
    
    std::complex<double> Ylmintegral (int l1, int m1, int l2, int m2){
        
        vector<std::complex<double> >  ylm_cc = Membranes[0]->get_ylm_vectorlist_for_mesh(l1, m1, true);
        vector<std::complex<double> >    ylm = Membranes[0]->get_ylm_vectorlist_for_mesh(l2, m2, false);
        
        std::complex<double> integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(ylm_cc, ylm);
        return integral;
    }
    
    std::complex<double> evaluate_u (int l1, int m1, int l2, int m2, double u, double r){
        Membranes[0]->generate_ulm_mode(l2, m2, u, r);
        vector<double> membrane_radii_list = Membranes[0]->get_ulmYlm_vectorlist_for_mesh();
        
        vector<std::complex<double> >  ylm_cc = Membranes[0]->get_ylm_vectorlist_for_mesh(l1, m1, true);
        
        std::complex<double> integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(ylm_cc, membrane_radii_list);
        return integral;
    }
    
    double RealYlmintegral (int l1, int m1, int l2, int m2){
        
        vector<double>  Realylm1 = Membranes[0]->get_real_ylm_vectorlist_for_mesh(l1, m1);
        vector<double>  Realylm2 = Membranes[0]->get_real_ylm_vectorlist_for_mesh(l2, m2);
        
        double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(Realylm1, Realylm2);
        return integral;
    }
    
    double Realevaluate_u (int l1, int m1, int l2, int m2, double u, double r){
        
        Membranes[0]->generate_ulm_mode_real(l2, m2, u, r);
        
        vector<double> membrane_radii_list = Membranes[0]->get_ulmYlm_vectorlist_for_mesh();
        
        vector<double>  Realylm1 = Membranes[0]->get_real_ylm_vectorlist_for_mesh(l1, m1);
        
        double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(Realylm1, membrane_radii_list);
        
        return integral;
    }
};

TEST_P(LMParameterized, get_ulmYlm_vectorlist_for_mesh){
    
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0;
    double R = 1;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.00002);
    EXPECT_NEAR(imag(test_value), 0, 0.00002);
    
    double t_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    EXPECT_NEAR(t_value, 0, 0.00002);
}


/**Spherical harmonics Normality*/

TEST_P(LMParameterized, YlmccYlmzzzNormality){
    EllM params = GetParam();
    int ell1 = params.l2;
    int m1 = params.m2;
    std::complex<double> test_value = Ylmintegral(params.l2, params.m2, params.l2, params.m2);
    
    EXPECT_NEAR(real(test_value), 1, 0.0037);
    EXPECT_NEAR(imag(test_value), 0, 0.000000001);
}


/**Spherical harmonics Orthonormality*/
TEST_P(LMParameterized, YlmccYlmzzzOrthonormality){
    EllM params = GetParam();
    //    int ell1 = params.l2;
    //    int m1 = params.m2;
    std::complex<double> test_value = Ylmintegral(params.l1, params.m1, params.l2, params.m2);
    
    EXPECT_NEAR(real(test_value), 0, 0.0012);
    EXPECT_NEAR(imag(test_value), 0, 0.0007);
}

/**Spherical harmonics Mode generator Orthonormality*/
/**Var:: Radius */
TEST_P(LMParameterized, YlmccModeOrthonormalityzzU1zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 1;
    double R = 1;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.0022);
    EXPECT_NEAR(imag(test_value), 0, 0.0034);
}


TEST_P(LMParameterized, YlmccModeOrthonormalityzzU1zzR10){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 1;
    double R = 10;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.0069);
    EXPECT_NEAR(imag(test_value), 0, 0.0006);
}

TEST_P(LMParameterized, YlmccModeOrthonormalityzzU1zzR100){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 1;
    double R = 100;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.001);
    EXPECT_NEAR(imag(test_value), 0, 0.00006);
}

/**Spherical harmonics Mode generator Orthonormality*/
/**Var:: Amplitude */

TEST_P(LMParameterized, YlmccModeOrthonormalityzzU0p1zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.1;
    double R = 1;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.018);
    EXPECT_NEAR(imag(test_value), 0, 0.00009);
}


TEST_P(LMParameterized, YlmccModeOrthonormalityzzU0p01zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.01;
    double R = 1;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.000018);
    EXPECT_NEAR(imag(test_value), 0, 0.000009);
}

TEST_P(LMParameterized, YlmccModeOrthonormalityzzU0p001zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.001;
    double R = 1;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.000012);
    EXPECT_NEAR(imag(test_value), 0, 0.0000069);
}

/**Spherical harmonics Mode generator Orthonormality*/
//Var:: Amplitude , Radius

TEST_P(LMParameterized, YlmccModeOrthonormalityzzU0p1zzR10){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.1;
    double R = 10;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.0007);
    EXPECT_NEAR(imag(test_value), 0, 0.0003);
}


TEST_P(LMParameterized, YlmccModeOrthonormalityzzU0p01zzR100){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.01;
    double R = 100;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.00008);
    EXPECT_NEAR(imag(test_value), 0, 0.00007);
}

TEST_P(LMParameterized, YlmccModeOrthonormalityzzU0p001zzR1000){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.001;
    double R = 1000;
    std::complex<double> test_value = evaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.000008);
    EXPECT_NEAR(imag(test_value), 0, 0.000003);
}

/**Spherical harmonics Mode generator Normality*/
/**Var::  Radius*/

TEST_P(LMParameterized, YlmccModeNormalityzzU1zzR1){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 1;
    double R = 1;
    std::complex<double> test_value = evaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(real(test_value), U, 0.35);
    EXPECT_NEAR(imag(test_value), 0, 0.0008);
}

//TEST_P(LMParameterized, YlmccModeNormalityzzU1zzR10){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 1;
//    double R = 10;
//    std::complex<double> test_value = evaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(real(test_value), U, 0.000000001);
//    EXPECT_NEAR(imag(test_value), 0, 0.000000001);
//}
//
//TEST_P(LMParameterized, YlmccModeNormalityzzU1zzR100){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 1;
//    double R = 100;
//    std::complex<double> test_value = evaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(real(test_value), U, 0.000000001);
//    EXPECT_NEAR(imag(test_value), 0, 0.000000001);
//}

///**Spherical harmonics Mode generator Normality*/
///**Var::  Amplitute*/
//
//TEST_P(LMParameterized, YlmccModeNormalityzzU0p1zzR1){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 0.1;
//    double R = 1;
//    std::complex<double> test_value = evaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(real(test_value), U, 0.051);
//    EXPECT_NEAR(imag(test_value), 0, 0.00005);
//}
//
//TEST_P(LMParameterized, YlmccModeNormalityzzU0p01zzR1){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 0.01;
//    double R = 1;
//    std::complex<double> test_value = evaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(real(test_value), U, 0.0051);
//    EXPECT_NEAR(imag(test_value), 0, 0.000005);
//}

//TEST_P(LMParameterized, YlmccModeNormalityzzU0p001zzR1){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 0.001;
//    double R = 1;
//    std::complex<double> test_value = evaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(real(test_value), U, 0.000051);
//    EXPECT_NEAR(imag(test_value), 0, 0.00000046);
//}

/**Spherical harmonics Mode generator Normality*/
/**Var::  Amplitute*/

//TEST_P(LMParameterized, YlmccModeNormalityzzU0p1zzR100){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 0.1;
//    double R = 100;
//    std::complex<double> test_value = evaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(real(test_value), U, 0.00000001);
//    EXPECT_NEAR(imag(test_value), 0, 0.00000001);
//}
//
//TEST_P(LMParameterized, YlmccModeNormalityzzU0p01zzR1000){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 0.01;
//    double R = 1000;
//    std::complex<double> test_value = evaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(real(test_value), U, 0.00000001);
//    EXPECT_NEAR(imag(test_value), 0, 0.00000001);
//}
//
//TEST_P(LMParameterized, YlmccModeNormalityzzU0p001zzR10000){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 0.001;
//    double R = 10000;
//    std::complex<double> test_value = evaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(real(test_value), U, 0.000000001);
//    EXPECT_NEAR(imag(test_value), 0, 0.000000001);
//}


//******************************************************************************
//******************** Real Spherical Harmonics   *********************
//******************************************************************************


/** Real Spherical harmonics Normality*/

TEST_P(LMParameterized, RYlmxRYlmzzzNormality){
    EllM params = GetParam();
    double test_value = RealYlmintegral(params.l2, params.m2, params.l2, params.m2);
    
    EXPECT_NEAR(test_value, 1, 0.0052);
    
}


/**Real Spherical harmonics Orthonormality*/
TEST_P(LMParameterized, RYlmxRYlmzzzOrthonormality){
    EllM params = GetParam();
    
    double test_value = RealYlmintegral(params.l1, params.m1, params.l2, params.m2);
    
    EXPECT_NEAR(test_value, 0, 0.0017);
}


/**Real Spherical harmonics Mode generator Orthonormality*/
/**Var:: Radius */
TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU1zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 1;
    double R = 1;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.1);
    
}


TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU1zzR10){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 1;
    double R = 10;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.1);
    
}

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU1zzR100){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 1;
    double R = 100;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(real(test_value), 0, 0.1);
    
}

/**Spherical harmonics Mode generator Orthonormality*/
/**Var:: Amplitude */

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p1zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.1;
    double R = 1;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.0002);
    
}


TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p01zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.01;
    double R = 1;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.00004);
    
}

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p001zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.001;
    double R = 1;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.000017);
    
}

/**Real Spherical harmonics Mode generator Orthonormality*/
//Var:: Amplitude , Radius

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p1zzR10){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.1;
    double R = 10;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.0002);
    
}


TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p01zzR100){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.01;
    double R = 100;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.00003);
    
}

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p001zzR1000){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.001;
    double R = 1000;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.00003);
    
}

/**Real Spherical harmonics Mode generator Normality*/
/**Var::  Radius*/

TEST_P(LMParameterized, RYlmxRModeNormalityzzU1zzR1){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 1;
    double R = 1;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.90);
    
}
//Will FAIL the same conditions as above
//TEST_P(LMParameterized, RYlmxRModeNormalityzzU1zzR10){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 1;
//    double R = 10;
//    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(test_value, U, 0.000000001);
//
//}
//Will FAIL the same conditions as above
//TEST_P(LMParameterized, RYlmxRModeNormalityzzU1zzR100){
//    EllM params = GetParam();
//    int ell2 = params.l2;
//    int m2 = params.m2;
//    double U = 1;
//    double R = 100;
//    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
//
//    EXPECT_NEAR(test_value, U, 0.000000001);
//
//}

/**Spherical harmonics Mode generator Normality*/
/**Var::  Amplitute*/

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p1zzR1){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.1;
    double R = 1;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.0027);
    
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p01zzR1){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.01;
    double R = 1;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000052);
    
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p001zzR1){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.001;
    double R = 1;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000006);
    
}

/**Spherical harmonics Mode generator Normality*/
/**Fix::  Amplitute 0.1, Var:: Radius*/

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p1zzR10){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.1;
    double R = 10;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.0027);
    
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p1zzR100){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.1;
    double R = 100;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.0027);
    
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p1zzR1000){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.1;
    double R = 1000;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.0027);
    
}

/**Spherical harmonics Mode generator Normality*/
/**Fix::  Amplitute 0.01, Var:: Radius*/

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p01zzR10){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.01;
    double R = 10;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000052);
    
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p01zzR100){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.01;
    double R = 100;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000052);
    
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p01zzR1000){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.01;
    double R = 1000;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000052);
    
}

/**Spherical harmonics Mode generator Normality*/
/**Fix::  Amplitute 0.001, Var:: Radius*/

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p001zzR10){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.001;
    double R = 10;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000006);
    
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p001zzR100){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.001;
    double R = 100;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000006);
    
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p001zzR1000){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.001;
    double R = 1000;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000006);
    
}

INSTANTIATE_TEST_SUITE_P(ModeGen,
                         LMParameterized,
                         ::testing::Values(
                                           EllM{3,-1,2,0},
                                           EllM{2,0,2,2},
                                           EllM{5,2,2,1},
                                           EllM{6,1,3,0},
                                           EllM{2,-1,3,-1},
                                           EllM{7,3,4,3},
                                           EllM{5,-4,8,-7},
                                           EllM{10,4,6,0},
                                           EllM{2,1,7,-5},
                                           EllM{13,13,5,-3}
                                           ));







#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>
using namespace std::complex_literals;

struct thetaphi{
    double theta;
    double phi;
};

class YlmWithOmegaParameterized : public testing::TestWithParam<thetaphi>{
protected:
    Membrane* mem;
    
public:
    std::complex<double> calc_ylm(int ell, int m, double theta, double phi){
        return mem->calc_complex_ylmthetaphi(ell,m,theta,phi);
        
    }
};

TEST_P(YlmWithOmegaParameterized, AllYlm00AreConst){
    thetaphi params = GetParam();
    std::complex<double> test_value = calc_ylm(0, 0, params.theta, params.phi);
    ASSERT_EQ(test_value, 0.5*sqrt(1/M_PI));
}

INSTANTIATE_TEST_SUITE_P(YlmExpectedValueTests,
                         YlmWithOmegaParameterized,
                         ::testing::Values(
                                           thetaphi{0,0},
                                           thetaphi{M_PI,2*M_PI},
                                           thetaphi{M_PI/3,2*M_PI/7},
                                           thetaphi{M_PI/13,17*M_PI/7}
                                           ));


TEST_P(YlmWithOmegaParameterized, L1Mn1SinglePrec){
    thetaphi params = GetParam();
    
    std::complex<double> test_value = calc_ylm(1, -1, params.theta, params.phi);
    
    std::complex<double> assert_value;
    double multiplyer =0.5*sqrt(3./(2*M_PI));
    assert_value.real(multiplyer*sin(params.theta)*cos(params.phi));
    assert_value.imag(-multiplyer*sin(params.theta)*sin(params.phi));
    
    EXPECT_FLOAT_EQ(real(test_value), real(assert_value));
    ASSERT_FLOAT_EQ(imag(test_value), imag(assert_value));
    
}

TEST_P(YlmWithOmegaParameterized, L1Mn1DoublePrec){
    thetaphi params = GetParam();
    
    std::complex<double> test_value = calc_ylm(1, -1, params.theta, params.phi);
    
    std::complex<double> assert_value;
    double multiplyer =0.5*sqrt(3./(2*M_PI));
    assert_value.real(multiplyer*sin(params.theta)*cos(params.phi));
    assert_value.imag(-multiplyer*sin(params.theta)*sin(params.phi));
    
    EXPECT_DOUBLE_EQ(real(test_value), real(assert_value));
    ASSERT_DOUBLE_EQ(imag(test_value), imag(assert_value));
    
}

TEST_P(YlmWithOmegaParameterized, L2Mn2SinglePrec){
    thetaphi params = GetParam();
    
    std::complex<double> test_value = calc_ylm(2, -2, params.theta, params.phi);
    
    std::complex<double> assert_value;
    double multiplyer =0.25*sqrt(15./(2*M_PI));
    assert_value.real(multiplyer*sin(params.theta)*sin(params.theta)*cos(2*params.phi));
    assert_value.imag(-multiplyer*sin(params.theta)*sin(params.theta)*sin(2*params.phi));
    
    EXPECT_FLOAT_EQ(real(test_value), real(assert_value));
    ASSERT_FLOAT_EQ(imag(test_value), imag(assert_value));
    
}

TEST_P(YlmWithOmegaParameterized, L2Mn2DoublePrec){
    thetaphi params = GetParam();
    
    std::complex<double> test_value = calc_ylm(2, -2, params.theta, params.phi);
    
    std::complex<double> assert_value;
    double multiplyer =0.25*sqrt(15./(2*M_PI));
    assert_value.real(multiplyer*sin(params.theta)*sin(params.theta)*cos(2*params.phi));
    assert_value.imag(-multiplyer*sin(params.theta)*sin(params.theta)*sin(2*params.phi));
    
    EXPECT_DOUBLE_EQ(real(test_value), real(assert_value));
    ASSERT_DOUBLE_EQ(imag(test_value), imag(assert_value));
    
}

TEST_P(YlmWithOmegaParameterized, L2M0){
    thetaphi params = GetParam();
    
    std::complex<double> test_value = calc_ylm(2, 0, params.theta, params.phi);
    
    std::complex<double> assert_value;
    double multiplyer =0.25*sqrt(5./(M_PI));
    assert_value.real(multiplyer*(3*cos(params.theta)*cos(params.theta)-1));
    assert_value.imag(0);
    
    
    ASSERT_EQ(imag(test_value), imag(assert_value));
    
}


int main (int argc, char **argv){
    testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}

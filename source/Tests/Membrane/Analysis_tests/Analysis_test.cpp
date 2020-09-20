#include <iostream>
#include <gtest/gtest.h>
#include <string>

#include "maps.hpp"
#include "Membrane.h"
#include "General_functions.hpp"
#include "Arg_pars.hpp"


namespace GenConst {
double Simulation_Time_In_Ps;
double Report_Interval_In_Fs;
double Step_Size_In_Fs;
int    MC_step;
int    Mem_fluidity;
bool   Periodic_box;
double Lbox;
double Simulation_box_length;
bool   Periodic_condtion_status;
int    Num_of_Membranes;
int    Num_of_Chromatins;
int    Num_of_Actins;
int    Num_of_ECMs;
string trajectory_file_name;;
string force_file_name;;
double Buffer_temperature; //***********OLDCODE
double Bussi_tau;
double Actin_Membrane_Bond_Coefficient;
bool   Interaction_map;
string Interaction_map_file_name;
bool   Excluded_volume_interaction;
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

std::vector<double> PeriodicBoxVector0;
std::vector<double> PeriodicBoxVector1;
std::vector<double> PeriodicBoxVector2;


//    std::vector<std::vector<std::vector<double> > > data;
std::vector<double> data_colection_times;
std::vector<std::vector<double> > Lboxdims;
}


ArgStruct_Analysis args;

using namespace std;
struct MemAnalysisTest : public testing::Test{
    
    vector<Membrane*> Membranes;
    
    MemAnalysisTest() {
        
        args.analysis_dim = 3;
        args.analysis_averaging_option = 0;
        args.num_ang_avg= 1;
        args.z_node=-1;
        args.zy_node=-1;
        args.analysis_filename = "10_1002n3frames.pdb";
        args.membane_labels.push_back("mem0");
        args.ell_max =20;
        args.output_filename.push_back("10_1002n3frames_ulmt_cpp.txt");
        args.framelimits_beg=2;
        args.framelimits_end=3;
        args.num_atoms_per_frame=1002;
        args.Mesh_files.push_back("10_1002n.ply");
        
        string general_file_name="General_param_map.txt";
        vector<string> membrane_config_list;
        vector<string> chromatin_config_list;
        vector<string> actin_config_list;
        vector<string> ecm_config_list;
        vector<string> pointparticle_config_list;
        
        
        GenConst::Testmode=true;
        read_general_parameters(general_file_name, membrane_config_list, chromatin_config_list, actin_config_list, ecm_config_list, pointparticle_config_list);
        
        
        GenConst::Num_of_Membranes=1;
        Membranes.push_back(new Membrane);
        
        args.framelimits_beg--;
        args.framelimits_end--;
        
        Membranes[0]->import_pdb_frames(args, 0);
        Membranes[0]->load_pdb_frame(0, args);
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

TEST_F( MemAnalysisTest, dOmegaSum){
    
    double integral =0;
    for (int i=0; i<Membranes[0]->get_num_of_nodes(); i++) {
        integral+=Membranes[0]->node_dOmega[i];
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
    
    args.analysis_filename = "10_1002nEllipsoid121_3frames.pdb";
    
    Membranes[0]->import_pdb_frames(args, 0);
    Membranes[0]->load_pdb_frame(0, args);
    
    double integral =0;
    for (int i=0; i<Membranes[0]->get_num_of_nodes(); i++) {
        integral+=Membranes[0]->node_voronoi_area[i];
    }
    EXPECT_NEAR(integral , 21.45963, 0.02);
}

TEST_F( SurfaceIntegral, OmegaEllipsoidxyz121){
    args.analysis_filename = "10_1002nEllipsoid121_3frames.pdb";
    
    Membranes[0]->import_pdb_frames(args, 0);
    Membranes[0]->load_pdb_frame(0, args);
    
    vector<double>  unit_vector;
    unit_vector.resize(Membranes[0]->get_num_of_nodes(), 1);
    double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(unit_vector, unit_vector);

    EXPECT_NEAR(integral , 4*M_PI,0.06);
}

TEST_F( SurfaceIntegral, VoronoiAreaSumEllipsoidxyz123){
    
    args.analysis_filename = "10_1002nEllipsoid123_3frames.pdb";
    
    Membranes[0]->import_pdb_frames(args, 0);
    Membranes[0]->load_pdb_frame(0, args);
    
    double integral =0;
    for (int i=0; i<Membranes[0]->get_num_of_nodes(); i++) {
        integral+=Membranes[0]->node_voronoi_area[i];
    }
    
    EXPECT_NEAR(integral , 48.93663, 0.15);
}

TEST_F( SurfaceIntegral, OmegaEllipsoidxyz123){
    args.analysis_filename = "10_1002nEllipsoid123_3frames.pdb";
    
    Membranes[0]->import_pdb_frames(args, 0);
    Membranes[0]->load_pdb_frame(0, args);
    
    vector<double>  unit_vector;
    unit_vector.resize(Membranes[0]->get_num_of_nodes(), 1);
    double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(unit_vector, unit_vector);
    EXPECT_NEAR(integral , 4*M_PI,0.09);
}

TEST_F( SurfaceIntegral, Omegapotatoe){
    args.analysis_filename = "potatoe_3frames.pdb";
    args.Mesh_files[0]="potatoe.ply";
    
    Membranes[0]->import_pdb_frames(args, 0);
    Membranes[0]->load_pdb_frame(0, args);
    
    vector<double>  unit_vector;
    unit_vector.resize(Membranes[0]->get_num_of_nodes(), 1);
    double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(unit_vector, unit_vector);
    EXPECT_NEAR(integral , 4*M_PI,0.042);
}
using reconstruction = MemAnalysisTest;
TEST_F( reconstruction, Potatoereconstruction){
    args.analysis_filename = "potatoe2_3frames.pdb";
    args.Mesh_files[0]="potatoe2.ply";
    Membranes[0]->import_pdb_frames(args, 0);
    Membranes[0]->load_pdb_frame(0, args);
    
    double radius = sqrt(Membranes[0]->surface_area_voronoi/(4*M_PI));
    
    Membranes[0]->calculate_real_ulm(args);
    
    Membranes[0]->generate_ulm_mode_real(0, 0, 0, radius);
    
    for (int ell=0; ell<args.ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            Membranes[0]->add_ulm_mode_real(ell, m, sqrt(Membranes[0]->ulm_avg[ell][m+ell]), radius);
        }
    }
    
//    for (int i=0; i<Membranes[0]->get_num_of_nodes(); i++) {
//        cout<<Membranes[0]->get_node_position(i,0)<<" "<<
//              Membranes[0]->get_node_position(i,1)<<" "<<
//              Membranes[0]->get_node_position(i,2)<<endl;
//    }
//    exit(0);
}

TEST_F( reconstruction, modezero){
    args.analysis_filename = "10_1002n3frames.pdb";
    args.Mesh_files[0]="10_1002n.ply";
    Membranes[0]->import_pdb_frames(args, 0);
    Membranes[0]->load_pdb_frame(0, args);
    double radius = sqrt(Membranes[0]->surface_area_voronoi/(4*M_PI));
    
    for (auto &coord: Membranes[0]->spherical_positions){
        coord[0] = radius;
        
    }
    for (int i=0; i<Membranes[0]->get_num_of_nodes(); i++) {
        Membranes[0]->Node_Position[i]=convert_spherical_to_cartesian(Membranes[0]->spherical_positions[i]);
    }
    Membranes[0]->calculate_real_ulm(args);
    

    for (int ell=0; ell<args.ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            EXPECT_NEAR(Membranes[0]->ulm_avg[ell][m+ell], 0, 0.00000000000001);
        }
    }
}

TEST_F( reconstruction, CompositeMode){
    args.analysis_filename = "10_1002n3frames.pdb";
    args.Mesh_files[0]="10_1002n.ply";
    Membranes[0]->import_pdb_frames(args, 0);
    Membranes[0]->load_pdb_frame(0, args);
    
    double radius = sqrt(Membranes[0]->surface_area_voronoi/(4*M_PI));
    args.ell_max=20;
    vector<double> us;
    vector<double> ells;
    vector<double> ms;
    
    us.push_back(0.01); ells.push_back(1); ms.push_back(0);
    us.push_back(0.02); ells.push_back(1); ms.push_back(1);
    us.push_back(0.03); ells.push_back(1); ms.push_back(-1);
    us.push_back(0.01); ells.push_back(2); ms.push_back(-1);
    us.push_back(0.02); ells.push_back(3); ms.push_back(2);
    us.push_back(0.03); ells.push_back(3); ms.push_back(-3);
    us.push_back(0.01); ells.push_back(4); ms.push_back(3);
    us.push_back(0.03); ells.push_back(5); ms.push_back(-2);
    us.push_back(0.01); ells.push_back(6); ms.push_back(0);
    us.push_back(0.02); ells.push_back(7); ms.push_back(3);
    us.push_back(0.01); ells.push_back(8); ms.push_back(8);
    us.push_back(0.01); ells.push_back(9); ms.push_back(-5);
    us.push_back(0.01); ells.push_back(10); ms.push_back(2);
    
    Membranes[0]->generate_ulm_mode_real(0, 0, 0, radius);
    for (int i=0; i<ells.size(); i++) {
        Membranes[0]->add_ulm_mode_real(ells[i], ms[i], us[i], radius);
    }
    
    Membranes[0] -> calculate_real_ulm(args);
    
    for (int i=0; i<ells.size(); i++) {
        EXPECT_NEAR(Membranes[0]->ulm_avg[ells[i]][ells[i]+ms[i]], us[i]*us[i], us[i]*us[i]*0.024);
    }
    for (int ell=11; ell<args.ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            EXPECT_NEAR(Membranes[0]->ulm_avg[ell][ell+m], 0, 0.0000001);
            
        }
    }
    
    Membranes[0]->generate_ulm_mode_real(0, 0, 0, radius);
    for (int ell=0; ell<args.ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            Membranes[0]->add_ulm_mode_real(ell, m, sqrt(Membranes[0]->ulm_avg[ell][m+ell]), radius);
        }
    }
    
    for (int i=0; i<Membranes[0]->get_num_of_nodes(); i++) {
        cout<<Membranes[0]->get_node_position(i,0)<<" "<<
        Membranes[0]->get_node_position(i,1)<<" "<<
        Membranes[0]->get_node_position(i,2)<<endl;
    }
    exit(0);
}
//******************************************************************************
//******************** Spherical Harmonics Mode Generator  *********************
//******************************************************************************


struct EllM{
    int l1; int m1; int l2; int m2;
};

class LMParameterized : public testing::TestWithParam<EllM>{
protected:
    vector<Membrane*> Membranes;
    
    LMParameterized() {
        args.analysis_dim = 3;
        args.analysis_averaging_option = 0;
        args.num_ang_avg= 1;
        args.z_node=-1;
        args.zy_node=-1;
        args.analysis_filename = "10_1002n3frames.pdb";
        args.membane_labels.push_back("mem0");
        args.ell_max =20;
        args.output_filename.push_back("10_1002n3frames_ulmt_cpp.txt");
        args.framelimits_beg=2;
        args.framelimits_end=3;
        args.num_atoms_per_frame=1002;
        args.Mesh_files.push_back("10_1002n.ply");
        
        string general_file_name="General_param_map.txt";
        vector<string> membrane_config_list;
        vector<string> chromatin_config_list;
        vector<string> actin_config_list;
        vector<string> ecm_config_list;
        vector<string> pointparticle_config_list;
        
        
        GenConst::Testmode=true;
        read_general_parameters(general_file_name, membrane_config_list, chromatin_config_list, actin_config_list, ecm_config_list, pointparticle_config_list);
        
        
        GenConst::Num_of_Membranes=1;
        Membranes.push_back(new Membrane);
        
        args.framelimits_beg--;
        args.framelimits_end--;
        
        Membranes[0]->import_pdb_frames(args, 0);
        Membranes[0]->load_pdb_frame(0, args);
        
    }
    ~LMParameterized() {
        for (auto pointer : Membranes){
            delete pointer;
        }
        Membranes.clear();
    }
    
    double RealYlmintegral (int l1, int m1, int l2, int m2){
        
        vector<double>  Realylm1 = Membranes[0]->get_real_ylm_vectorlist_for_mesh(l1, m1);
        vector<double>  Realylm2 = Membranes[0]->get_real_ylm_vectorlist_for_mesh(l2, m2);
        
        double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(Realylm1, Realylm2);
        return integral;
    }
    
    double Realevaluate_u (int l1, int m1, int l2, int m2, double u, double r){
        
        Membranes[0]->generate_ulm_mode_real(l2, m2, 0, r);
        Membranes[0]->add_ulm_mode_real(l2, m2, u, r);
        
        vector<double> membrane_radii_list = Membranes[0]->get_ulmYlm_vectorlist_for_mesh();
        
        vector<double>  Realylm1 = Membranes[0]->get_real_ylm_vectorlist_for_mesh(l1, m1);
        
        double integral = Membranes[0]->calc_vectorlist_vectorlist_surface_integral(Realylm1, membrane_radii_list);
        
        return integral;
    }
};


/** Real Spherical harmonics Normality*/
TEST_P(LMParameterized, RYlmxRYlmzzzNormality){
    exit(0);
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
    
    EXPECT_NEAR(test_value, 0, 0.0017);
}

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU1zzR10){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 1;
    double R = 10;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.0017);
}

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU1zzR100){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 1;
    double R = 100;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.0017);
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
        
    EXPECT_NEAR(test_value, 0, 0.00015);
}

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p01zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.01;
    double R = 1;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.000021);
}

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p001zzR1){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.001;
    double R = 1;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.000016);
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
    
    EXPECT_NEAR(test_value, 0, 0.00015);
}

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p01zzR100){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.01;
    double R = 100;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.000021);
}

TEST_P(LMParameterized, RYlmxRModeOrthonormalityzzU0p001zzR1000){
    EllM params = GetParam();
    int ell1 = params.l1;
    int m1 = params.m1;
    double U = 0.001;
    double R = 1000;
    double test_value = Realevaluate_u(ell1, m1, params.l2, params.m2, U, R);
    
    EXPECT_NEAR(test_value, 0, 0.000016);
}

/**Real Spherical harmonics Mode generator Normality*/
/**Var::  Radius*/
TEST_P(LMParameterized, RYlmxRModeNormalityzzU1zzR1FAILURE){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 1;
    double R = 1;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_GT(U-test_value, 0.11);
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p1zzR1){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.1;
    double R = 1;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);

    EXPECT_NEAR(test_value, U, 0.0015);    
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p01zzR1){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.01;
    double R = 1;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000037);
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p001zzR1){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.001;
    double R = 1;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000005);
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
    
    EXPECT_NEAR(test_value, U, 0.0015);
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p1zzR100){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.1;
    double R = 100;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.0015);
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p1zzR1000){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.1;
    double R = 1000;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);

    EXPECT_NEAR(test_value, U, 0.0015);
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

    EXPECT_NEAR(test_value, U, 0.000037);
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p01zzR100){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.01;
    double R = 100;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000037);
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p01zzR1000){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.01;
    double R = 1000;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000037);
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
    
    EXPECT_NEAR(test_value, U, 0.000005);
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p001zzR100){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.001;
    double R = 100;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000005);
}

TEST_P(LMParameterized, RYlmxRModeNormalityzzU0p001zzR1000){
    EllM params = GetParam();
    int ell2 = params.l2;
    int m2 = params.m2;
    double U = 0.001;
    double R = 1000;
    double test_value = Realevaluate_u(ell2, m2, ell2, m2, U, R);
    
    EXPECT_NEAR(test_value, U, 0.000005);
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

int main (int argc, char **argv){
    testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}

//
//  OpenMM_structs.h
//  Membrae
//
//  Created by Ali Farnudi on 10/06/2019.
//  Copyright © 2019 Ali Farnudi. All rights reserved.
//

#ifndef OpenMM_structs_h
#define OpenMM_structs_h

#include "OpenMM.h"

/** OpenMM's atom information format*/
struct MyAtomInfo
{
    int type;
    int Vsite_particleindecies[2];
    double Vsite_weights[2];
    char* pdb;
    char symbol;
    double initPosInNm[3];
    double posInNm[3];
    double velocityInNmperPs[3];
    double mass;
    double radius;
    double sigma_LJ_12_6;
    double epsilon_LJ_12_6;
    double force[3];
    double energyInKJ;
    double stretching_energy;
    int ext_force_model;
    double ext_force_constants[3];
    std::string class_label;
    int vsite_atoms[2];
    double vsite_weights[2];
    double epsilonWCA   = 0;
    double sigmaWCA   = 0;
};

struct Bonds{
    int type;
    bool globalStat;
    int atoms[2];
    std::string class_label;
   // double nominalLengthInAngstroms, stiffnessInKcalPerAngstrom2, stiffnessInKcalPerAngstrom4;
   // double dampInKcalPsPerAngstrom2;
   // double FENE_lmax, FENE_lmin, FENE_le0, FENE_le1;
    double nominalLengthInNm, stiffnessInKJPerNm2, stiffnessInKJPerNm4;
    double dampInKJPsPerNm2;
    double FENER0inNm, epsilon_WCA_inKJpermol, k_FENE_inKJpermol;
    double F0;
    double r_min,r_max;
    double hill_co;
    double k_F0;
    bool   canConstrain;
    double gompperlmin;
    double gompperlc1;
    double gompperlc0;
    double gompperlmax;
    
    double ellipsoidLockXscale;
    double ellipsoidLockYscale;
    double ellipsoidLockZscale;
    double ellipsoidLockRscale;
    
    double LockOn_rigidity;
    double LockOnULM_amplitude;
};

struct Dihedrals{
    int type;
    std::string class_label;
    std::vector<int> atoms;
    double bendingStiffnessinKJ;
    double spontaneousBendingAngleInRad;
    double total_mem_area;
};

struct Angles{
    int type;
    std::string class_label;
    int atoms[3];
    double bendingStiffnessinKJpermol;
    double spontaneousBendingAngleInRad;
};

struct Triangles{
    int surface_type;
    int volume_type;
    std::string class_label;
    std::vector <int> atoms;
    double SurfaceConstraintStiffnessinKJpermolperNm2;
    double SurfaceConstraintValue;
    
    double VolumeConstraintStiffnessinKJpermolperNm3;
    double VolumeConstraintValue;
};

struct MeanCurvature{
    int node_order;
    int curvature_type;
    std::string class_label;
    std::vector<int> atoms;
    double curvatureStiffnessinKJpermol;
    double spontaneousCurvature;
};

struct PlatformInfo{
    int platform_id=0;
    int platform_device_id=0;
    
    std::vector<std::string> device_properties_report;
    std::vector<std::map<std::string, std::string> > device_properties;
};

/** -----------------------------------------------------------------------------
 *                           INTERFACE TO OpenMM
 * -----------------------------------------------------------------------------
 * These four functions and an opaque structure are used to interface our main
 * program with OpenMM without the main program having any direct interaction
 /Users/sajjad/virtual cell/Membrane_OBJ/source/neighbour_pool_constructor.cpp * with the OpenMM API. This is a clean approach for interfacing with any MD
 * code, although the details of the interface routines will differ. This is
 * still just "locally written" code and is not required by OpenMM.
 *
 * This is our opaque "handle" class containing all the OpenMM objects that
 * must persist from call to call during a simulation. The main program gets
 * a pointer to one of these but sees it as essentially a void* since it
 * doesn't know the definition of this class.
 */
struct MyOpenMMData {
    MyOpenMMData() : system(0), context(0), integrator(0) {}
    ~MyOpenMMData() {delete context; delete integrator; delete system;}
    OpenMM::System*         system;
    OpenMM::Integrator*           integrator;
    OpenMM::VerletIntegrator*     VerletIntegrator;
    OpenMM::BrownianIntegrator*   BrownianIntegrator;
    OpenMM::LangevinIntegrator*   LangevinIntegrator;
    
    OpenMM::CustomIntegrator*     CustomIntegrator;
//    OpenMM::CustomIntegrator*     LangevinMinimisation;
    
    OpenMM::Context*  context;
    OpenMM::Context*  minimisationContext;
    
    PlatformInfo platforminfo;
    
    OpenMM::HarmonicBondForce* harmonic;
    OpenMM::HarmonicBondForce* calcforce;
    std::vector<OpenMM::CustomBondForce*> x4harmonic;
    std::vector<OpenMM::CustomCompoundBondForce*> Dihedral;
    std::vector<OpenMM::CustomCompoundBondForce*> GlobalSurfaceConstraintForces;
    std::vector<OpenMM::CustomCompoundBondForce*> LocalSurfaceConstraintForces;
    std::vector<OpenMM::CustomCompoundBondForce*> GlobalVolumeConstraintForces;
    std::vector<OpenMM::CustomCompoundBondForce*> MeanCurvatureForces;
    std::vector<OpenMM::CustomAngleForce*> Angle;
    std::vector<OpenMM::CustomNonbondedForce*> LJ;
    std::vector<OpenMM::CustomNonbondedForce*> EV;
    std::vector<OpenMM::CustomNonbondedForce*> WCA;
    std::vector<OpenMM::CustomNonbondedForce*> WCAFC;
    
    
};


struct TimeDependantData {
    OpenMM::HarmonicBondForce*  Kelvin_VoigtBond;
    std::vector<OpenMM::CustomBondForce*> Hill_force;
    std::vector<OpenMM::CustomBondForce*> k_force;
    int hill_stepnum = 40;
    bool HillForce = false;
    int k_stepnum = 120;
    bool kForce = false;
    double COM[3];
    bool Kelvin_Voigt = false;
    int Kelvin_stepnum = 40;
    std::vector<double> Kelvin_Voigt_damp;
    std::vector<double> hill_const_force;
    std::vector<double> hill_coefficient;
    std::vector<std::vector<double> > Kelvin_Voigt_distInNm;
    std::vector<double> Kelvin_Voigt_initNominal_length_InNm;
    std::vector<std::vector<std::vector<double> > > hill_distInNm;
    //std::vector<std::vector<double>> Contractile_initForce;
    
    std::vector<OpenMM::CustomExternalForce*> ext_force;
    //OpenMM::CustomExternalForce* ext_force;
    //bool force_update = false;
    //int force_stepnum = 100000;
  
    void Kelvin_Nominal_length_calc()
    {
        if(Kelvin_Voigt)
        {
            std::vector<double> Nominal_Length_InNm ;
            
            const int Num_Bonds = Kelvin_VoigtBond->getNumBonds();
            int atom1, atom2 ;
            double length, stiffness;
            
            for(int i=0; i<Num_Bonds ; ++i)
            {
                Kelvin_VoigtBond->getBondParameters(i, atom1, atom2, length, stiffness);
                Nominal_Length_InNm.push_back(length);
            }
            Kelvin_Voigt_initNominal_length_InNm = Nominal_Length_InNm;
        }
    }
    
    
    
    void Kelvin_dist_calc(MyAtomInfo atoms[])
    {
        if(Kelvin_Voigt)
        {
            std::vector<double> distInNm ;
            const int Num_Bonds = Kelvin_VoigtBond->getNumBonds();
            int atom1, atom2 ;
            double length, stiffness , dist;
            
            for(int i=0; i<Num_Bonds; i++)
            {
                Kelvin_VoigtBond->getBondParameters(i, atom1, atom2, length, stiffness);
                dist =sqrt ( ( atoms[atom1].posInNm[0] - atoms[atom2].posInNm[0]) * ( atoms[atom1].posInNm[0] - atoms[atom2].posInNm[0]) + ( atoms[atom1].posInNm[1] - atoms[atom2].posInNm[1]) * ( atoms[atom1].posInNm[1] - atoms[atom2].posInNm[1]) + ( atoms[atom1].posInNm[2] - atoms[atom2].posInNm[2]) * ( atoms[atom1].posInNm[2] - atoms[atom2].posInNm[2]) ) ;
                
                distInNm.push_back(dist);
            }
            
            
            Kelvin_Voigt_distInNm.push_back(distInNm);
        }
    }
    
    
    
    
    void COM_calculator(MyAtomInfo atoms[])
    {
        double total_mass = 0;
        for(int i=0; atoms[i].type != -1 ; i++)
            {
                if(atoms[i].class_label.compare("ECM") != 0)
                {
                    COM[0] = COM[0] + atoms[i].mass * atoms[i].posInNm[0];
                    COM[1] = COM[1] + atoms[i].mass * atoms[i].posInNm[1];
                    COM[2] = COM[2] + atoms[i].mass * atoms[i].posInNm[2];
                    total_mass = total_mass + atoms[i].mass ;
                }
                
            }
        
        COM[0] = COM[0] / total_mass ;
        COM[1] = COM[1] / total_mass ;
        COM[2] = COM[2] / total_mass ;
            
            
    }
    
    
    
    
    
    
    
    
    
    void hill_dist_calc(MyAtomInfo atoms[])
    {
        if(HillForce)
        {
            std::vector<std::vector<double>> all_distInNm;
            for(int j=0; j<Hill_force.size() ; j++)
            {
            std::vector<double> distInNm ;
            distInNm.clear();
            const int Num_Bonds = Hill_force[j]->getNumBonds();
            int atom1, atom2 ;
            double dist;
            std::vector<double> parameters;
            
            for(int i=0; i<Num_Bonds; i++)
            {
                Hill_force[j]->getBondParameters(i, atom1, atom2, parameters);
                dist =sqrt ( ( atoms[atom1].posInNm[0] - atoms[atom2].posInNm[0]) * ( atoms[atom1].posInNm[0] - atoms[atom2].posInNm[0]) + ( atoms[atom1].posInNm[1] - atoms[atom2].posInNm[1]) * ( atoms[atom1].posInNm[1] - atoms[atom2].posInNm[1]) + ( atoms[atom1].posInNm[2] - atoms[atom2].posInNm[2]) * ( atoms[atom1].posInNm[2] - atoms[atom2].posInNm[2]) ) ;
                
                distInNm.push_back(dist);
            }
                all_distInNm.push_back(distInNm);
            
            }
            
            hill_distInNm.push_back(all_distInNm);
            
        }
    }
    
    
    
};




#endif /* OpenMM_structs_h */

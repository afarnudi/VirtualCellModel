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
    char* pdb;
    char symbol;
    double initPosInAng[3];
    double posInAng[3];
    double velocityInAngperPs[3];
    double mass;
    double radius;
    double sigma_LJ_12_6;
    double epsilon_LJ_12_6;
    double force[3];
    double energy;
    double stretching_energy;
    int ext_force_model;
    double ext_force_constants[3];
};

struct Bonds{
    int type;
    int atoms[2];
    std::string class_label;
    double nominalLengthInAngstroms, stiffnessInKcalPerAngstrom2, stiffnessInKcalPerAngstrom4;
    double dampInKcalPsPerAngstrom2;
    double FENE_lmax, FENE_lmin, FENE_le0, FENE_le1;
    bool   canConstrain;
};

struct Dihedrals{
    int type;
    std::string class_label;
    std::vector<int> atoms;
    double bending_stiffness_value;
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
    OpenMM::Integrator*     integrator;
    OpenMM::Context*  context;

    OpenMM::HarmonicBondForce* harmonic;
    std::vector<OpenMM::CustomBondForce*> x4harmonic;
    std::vector<OpenMM::CustomCompoundBondForce*> Dihedral;

    std::vector<OpenMM::CustomNonbondedForce*> EV;
};


struct TimeDependantData {
    OpenMM::HarmonicBondForce*  Kelvin_VoigtBond;
    bool Kelvin_Voigt = false;
    int Kelvin_stepnum = 100;
    std::vector<double> Kelvin_Voigt_damp;
    std::vector<std::vector<double>> Kelvin_Voigt_distInAng;
    std::vector<double> Kelvin_Voigt_initNominal_length_InNm;
    
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
            std::vector<double> distInAng ;
            const int Num_Bonds = Kelvin_VoigtBond->getNumBonds();
            int atom1, atom2 ;
            double length, stiffness , dist;
            
            for(int i=0; i<Num_Bonds; i++)
            {
                Kelvin_VoigtBond->getBondParameters(i, atom1, atom2, length, stiffness);
                dist =sqrt ( ( atoms[atom1].posInAng[0] - atoms[atom2].posInAng[0]) * ( atoms[atom1].posInAng[0] - atoms[atom2].posInAng[0]) + ( atoms[atom1].posInAng[1] - atoms[atom2].posInAng[1]) * ( atoms[atom1].posInAng[1] - atoms[atom2].posInAng[1]) + ( atoms[atom1].posInAng[2] - atoms[atom2].posInAng[2]) * ( atoms[atom1].posInAng[2] - atoms[atom2].posInAng[2]) ) ;
                
                distInAng.push_back(dist);
            }
            
            
            Kelvin_Voigt_distInAng.push_back(distInAng);
        }
    }
    
    
};




#endif /* OpenMM_structs_h */

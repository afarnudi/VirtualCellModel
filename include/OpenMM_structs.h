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
    double initPosInAng[3];
    double posInAng[3];
    double mass;
    double radius;
    double sigma_LJ_12_6;
    double epsilon_LJ_12_6;
};

struct Bonds{
    int type;
    int atoms[2];
    std::string class_label;
    double nominalLengthInAngstroms, stiffnessInKcalPerAngstrom2, stiffnessInKcalPerAngstrom4;
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
 * with the OpenMM API. This is a clean approach for interfacing with any MD
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
    

};




#endif /* OpenMM_structs_h */

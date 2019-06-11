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
    const char* pdb;
    double initPosInAng[3];
    double posInAng[3];
    double mass;
};

struct Bonds{
    int type;
    int atoms[2];
    double nominalLengthInAngstroms, stiffnessInKcalPerAngstrom2;
    bool   canConstrain;
};

struct Dihedrals{
    int type;
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
};


/** This function and an opaque structure are used to interface our main
 * programme with OpenMM without the main programme having any direct interaction
 * with the OpenMM API.
 * This function initilises the openmm system + contex + forces.
 */
MyOpenMMData* myInitializeOpenMM(const MyAtomInfo       atoms[],
                                        double          stepSizeInFs,
                                        std::string&    platformName,
                                        Bonds*          bonds,
                                        Dihedrals*      dihedrals);


using OpenMM::Vec3;
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

#endif /* OpenMM_structs_h */

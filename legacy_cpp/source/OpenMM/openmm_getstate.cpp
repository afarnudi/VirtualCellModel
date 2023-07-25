#include "Membrane.h"
#include "General_functions.hpp"
#include "Global_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include <vector>
#include <stdlib.h>
#include <math.h>


using OpenMM::Vec3;
void myGetOpenMMState(OpenMM::Context* context,
                      double& timeInPs,
                      double& energyInKJ,
                      double& potential_energyInKJ,
                      MyAtomInfo atoms[])
{

    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    infoMask += OpenMM::State::Velocities;  // for kinetic energy (cheapm)
    if (generalParameters.WantEnergy) {
        infoMask += OpenMM::State::Energy;     // for pot. energy (expensive)
    }
    // Forces are also available (and cheap).
    const OpenMM::State state = context->getState(infoMask);
//    const OpenMM::State state = omm->context->getState(infoMask,GenConst::Periodic_box);
    
    timeInPs = state.getTime(); // OpenMM time is in ps already
    
    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    const std::vector<Vec3>& velInNmperPs  = state.getVelocities();
    
//    if (generalParameters.Periodic_condtion_status) {
//        Vec3 Lboxx, Lboxy, Lboxz;
//        state.getPeriodicBoxVectors(Lboxx, Lboxy, Lboxz);
//        
//        
//        GenConst::Lboxdims.clear();
//        GenConst::Lboxdims.resize(3);
//        for (int i=0; i<3; i++) {
//            GenConst::Lboxdims[i].resize(3,0);
//        }
//        for (int j=0; j<3; j++) {
//            GenConst::Lboxdims[0][j]=Lboxx[j];
//            GenConst::Lboxdims[1][j]=Lboxy[j];
//            GenConst::Lboxdims[2][j]=Lboxz[j];
//        }
//    }
    
    for (int i=0; i < (int)positionsInNm.size(); ++i){
        for (int j=0; j < 3; ++j){
            atoms[i].posInNm[j] = positionsInNm[i][j];
            atoms[i].velocityInNmperPs[j] = velInNmperPs[i][j];
            
        }
    }
//    if (GenConst::WantForce) {
//        const std::vector<Vec3>& Forces = state.getForces();
//        for (int i=0; i < (int)positionsInNm.size(); ++i){
//            for (int j=0; j < 3; ++j){
//                atoms[i].posInNm[j] = positionsInNm[i][j];
//                atoms[i].velocityInNmperPs[j] = velInNmperPs[i][j];
//                if (GenConst::WantForce) {
//                    atoms[i].force[j]    = Forces[i][j];
//                }
//            }
//        }
//    }
    
    // If energy has been requested, obtain it in kJ/mol.
    energyInKJ = 0;

    if (generalParameters.WantEnergy){
        energyInKJ = state.getPotentialEnergy() + state.getKineticEnergy();
        potential_energyInKJ = state.getPotentialEnergy();
    }
    
}



void Cheap_GetOpenMMState(OpenMM::Context* context,
                          MyAtomInfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    const OpenMM::State state = context->getState(infoMask);
    
    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    
    for (int i=0; i < (int)positionsInNm.size(); ++i){
        for (int j=0; j < 3; ++j){
            atoms[i].posInNm[j] = positionsInNm[i][j];
        }
    }
}

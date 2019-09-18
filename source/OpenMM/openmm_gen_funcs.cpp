#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include <vector>

void myStepWithOpenMM(MyOpenMMData* omm, int numSteps) {
    omm->integrator->step(numSteps);
}

void myTerminateOpenMM(MyOpenMMData* omm) {
    delete omm;
}

using OpenMM::Vec3;
void myGetOpenMMState(MyOpenMMData* omm, bool wantEnergy,
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


using OpenMM::Vec3;
OpenMM::State getCurrentState(MyOpenMMData* omm, bool wantEnergy,
                 double& timeInPs,  bool wantforce)
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    if (wantEnergy) {
        infoMask += OpenMM::State::Velocities; // for kinetic energy (cheap)
        infoMask += OpenMM::State::Energy;     // for pot. energy (expensive)
    }
    // Forces are also available (and cheap).
    if (wantforce){
        infoMask+=OpenMM::State::Forces;
    }
    OpenMM::State currentstate = omm->context->getState(infoMask);
    timeInPs = currentstate.getTime(); // OpenMM time is in ps already
    return(currentstate);
}

void setNewState(MyOpenMMData* omm, bool wantEnergy,
                 double& energyInKcal,
                 MyAtomInfo atoms[], bool wantforce)
{   
    double timeInPs;
    const OpenMM::State newstate= getCurrentState(omm, wantEnergy,timeInPs,  wantforce);
        // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = newstate.getPositions();
    for (int i=0; i < (int)positionsInNm.size(); ++i){
        for (int j=0; j < 3; ++j){
    atoms[i].posInAng[j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;}}
    
    // If energy has been requested, obtain it and convert from kJ to kcal.
    energyInKcal = 0;
    if (wantEnergy)
        energyInKcal = (newstate.getPotentialEnergy() + newstate.getKineticEnergy())
        * OpenMM::KcalPerKJ;
}
using OpenMM::Vec3;
void          myreinitializeOpenMMState(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals){
    bool preservestate=1;
    
    //OpenMM::HarmonicBondForce*   HarmonicBond = omm->system->getForce(0);
    vector<double> parameters;
    vector<double> bendingparameter ;
    vector<int>particles ={2 ,  1, 3, 0};
    parameters.resize(2);
    
    int bond_index=2;
    int particle1=1;
    int particle2=3;
    double r_rest= bonds[0].nominalLengthInAngstroms*OpenMM::NmPerAngstrom;
    double k_bond= bonds[0].stiffnessInKcalPerAngstrom4
                                                 * OpenMM::KJPerKcal
                                                 * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
    double K_bend=dihedrals[0].bending_stiffness_value* OpenMM::KJPerKcal
                                                      * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
    parameters[0]=r_rest;
    parameters[1]=k_bond;
    bendingparameter.push_back(K_bend);
    //omm->x4harmonic->setBondParameters(bond_index,particle1, particle2, parameters);
    omm->Dihedral->setBondParameters(0,particles, bendingparameter);
    double lenght= bonds[0].nominalLengthInAngstroms* OpenMM::NmPerAngstrom;
    double k= bonds[0].stiffnessInKcalPerAngstrom2* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
    omm->harmonic->setBondParameters(bond_index,particle1, particle2, lenght, k );
    //bond_index=2;
    //particle1=1;
    //particle2=2;
    //omm->harmonic->setBondParameters(bond_index,particle1, particle2, lenght, k );
    omm->context->reinitialize(preservestate);
    
}

//                               PDB FILE WRITER
// Given state data, output a single frame (pdb "model") of the trajectory.
void myWritePDBFrame(int frameNum, double timeInPs, double energyInKcal,
                const MyAtomInfo atoms[], std::string traj_name)
{
    int EndOfList=-1;
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"MODEL     %d\n", frameNum);
    fprintf(pFile,"REMARK 250 time=%.3f ps; energy=%.3f kcal/mole\n",
            timeInPs, energyInKcal);
    
    for (int n=0; atoms[n].type != EndOfList; ++n){
        
        fprintf(pFile,"ATOM  %5d %4s ETH     1    %8.3f%8.3f%8.3f  1.00  0.00\n",
                n+1, atoms[n].pdb,
                atoms[n].posInAng[0], atoms[n].posInAng[1], atoms[n].posInAng[2]);
    }
    fprintf(pFile,"ENDMDL\n");
    fclose (pFile);
}

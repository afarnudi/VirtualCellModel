#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"


void set_multithermos(MyOpenMMData* omm, NonBondInteractionMap  &interaction_map, double stepSizeInFs, vector<set<int> >      &membrane_set, const MyAtomInfo  atoms[]){
    
    OpenMM::CustomNonbondedForce* thermoforce;
    thermoforce = new OpenMM::CustomNonbondedForce("r");
//    thermoforce->addGlobalParameter("tempa", 1);
    thermoforce->addInteractionGroup(membrane_set[0], membrane_set[0]);
//    thermoforce->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
//    thermoforce->addEnergyParameterDerivative("tempa");
    thermoforce->setForceGroup(31);
    
    
    int EndOfList=-1;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        thermoforce->addParticle();
    }
    omm->system->addForce(thermoforce);
    
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionInPs;
    double kBT = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("a", exp(-0.5*friction*dt));
    omm->CustomIntegrator->addGlobalVariable("b", sqrt(1-exp(-friction*dt)));
    omm->CustomIntegrator->addGlobalVariable("c", (1-exp(-0.5*friction*dt))/friction );
    omm->CustomIntegrator->addGlobalVariable("kTa", kBT);
    omm->CustomIntegrator->addGlobalVariable("kTb", kBT*5);
//    omm->CustomIntegrator->addGlobalVariable("tempa", 1);
    omm->CustomIntegrator->addPerDofVariable("stat", 0);
    
    omm->CustomIntegrator->addUpdateContextState();
    omm->CustomIntegrator->addComputePerDof("stat", "f31");
//    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f/m + b*sqrt(kTa/m)*gaussian");
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f0/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + (1-delta(stat))*b*sqrt(kTb/m)*gaussian");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(stat)*b*sqrt(kTa/m)*gaussian");
//    omm->CustomIntegrator->addComputePerDof("v", "v + (1-deriv(energy31, tempa))*b*sqrt(kTb/m)*gaussian");
    omm->CustomIntegrator->addComputePerDof("x", "x + dt*v");
//    omm->CustomIntegrator->addComputePerDof("v", "v + stat*(a*v - v + c*f/m + b*sqrt(kTa/m)*gaussian)");
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f0/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + (1-delta(stat))*b*sqrt(kTb/m)*gaussian");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(stat)*b*sqrt(kTa/m)*gaussian");
//    omm->CustomIntegrator->addComputePerDof("v", "v + (1-deriv(energy31, tempa))*b*sqrt(kTb/m)*gaussian");
    
}

void set_customLangevin(MyOpenMMData* omm, double stepSizeInFs){
    
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionInPs;
    double kBT = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("a", exp(-0.5*friction*dt));
    omm->CustomIntegrator->addGlobalVariable("b", sqrt(1-exp(-friction*dt)));
    omm->CustomIntegrator->addGlobalVariable("c", (1-exp(-0.5*friction*dt))/friction );
    omm->CustomIntegrator->addGlobalVariable("kT", kBT);
    
    omm->CustomIntegrator->addUpdateContextState();
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f/m + b*sqrt(kT/m)*gaussian");
    omm->CustomIntegrator->addComputePerDof("x", "x + dt*v");
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f/m + b*sqrt(kT/m)*gaussian");
}

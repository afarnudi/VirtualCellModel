#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"


void set_multithermos(MyOpenMMData* omm, NonBondInteractionMap  &interaction_map, double stepSizeInFs, vector<set<int> >      &membrane_set, const MyAtomInfo  atoms[]){
    /**
                Leap-frog Langevin integrator. Reference:
                Jesús A. Izaguirre, Chris R. Sweet, and Vijay S. Pande. Multiscale dynamics of macromolecules using Normal Mode Langevin. Pacific Symposium on Biocomputing, 15:240–251, 2010.
     */
    set<int> memSetWithOneElement;
    int setint = *membrane_set[0].begin();
    memSetWithOneElement.insert(setint);
    OpenMM::CustomNonbondedForce* thermoforce;
    thermoforce = new OpenMM::CustomNonbondedForce("r");
    thermoforce->addInteractionGroup(membrane_set[0], memSetWithOneElement);
    thermoforce->setForceGroup(31);
        
    int EndOfList=-1;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        thermoforce->addParticle();
    }
    omm->system->addForce(thermoforce);
    
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionInPs;
    double kBT    = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
    double kBTmem = generalParameters.BoltzmannKJpermolkelvin*generalParameters.customtemperature;
    
    
    int groupcount = interaction_map.get_ForceGroupCount();
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("a", exp(-0.5*friction*dt));
    omm->CustomIntegrator->addGlobalVariable("b", sqrt(1-exp(-friction*dt)));
    omm->CustomIntegrator->addGlobalVariable("c", (1-exp(-0.5*friction*dt))/friction );
    omm->CustomIntegrator->addGlobalVariable("kTa", kBT);
    omm->CustomIntegrator->addGlobalVariable("kTb", kBTmem);
    omm->CustomIntegrator->addPerDofVariable("stat", 0);
    
    
    omm->CustomIntegrator->addUpdateContextState();
    omm->CustomIntegrator->addComputePerDof("stat", "f31");
    
    //integrate all bonded forces
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f0/m");
    //Add
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(stat)*c*f1/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + (1-delta(stat))*b*sqrt(kTb/m)*gaussian");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(stat)*b*sqrt(kTa/m)*gaussian");
//    omm->CustomIntegrator->addComputePerDof("v", "v + (1-deriv(energy31, tempa))*b*sqrt(kTb/m)*gaussian");
    omm->CustomIntegrator->addComputePerDof("x", "x + dt*v");
//    omm->CustomIntegrator->addComputePerDof("v", "v + stat*(a*v - v + c*f/m + b*sqrt(kTa/m)*gaussian)");
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f0/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(stat)*c*f1/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + (1-delta(stat))*b*sqrt(kTb/m)*gaussian");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(stat)*b*sqrt(kTa/m)*gaussian");
//    omm->CustomIntegrator->addComputePerDof("v", "v + (1-deriv(energy31, tempa))*b*sqrt(kTb/m)*gaussian");
    
}

void set_customLangevinforminimisation(MyOpenMMData* omm, double stepSizeInFs){
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionInPs;
    double kBT    = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
    
    omm->LangevinMinimisation = new OpenMM::CustomIntegrator(dt);
    
    
    omm->LangevinMinimisation->addGlobalVariable("a", exp(-0.5*friction*dt));
    omm->LangevinMinimisation->addGlobalVariable("b", sqrt(1-exp(-friction*dt)));
    omm->LangevinMinimisation->addGlobalVariable("c", (1-exp(-0.5*friction*dt))/friction );
    omm->LangevinMinimisation->addGlobalVariable("kT", kBT);
    
    omm->LangevinMinimisation->addUpdateContextState();
    
//    omm->LangevinMinimisation->addComputePerDof("v", "a*v + c*f/m + b*sqrt(kT/m)*gaussian");
    omm->LangevinMinimisation->addComputePerDof("v", " c*f/m + b*sqrt(kT/m)*gaussian");
    
    omm->LangevinMinimisation->addComputePerDof("x", "x + dt*v");

    //    omm->LangevinMinimisation->addComputePerDof("v", "a*v + c*f/m + b*sqrt(kT/m)*gaussian");
}

void set_Langevin(MyOpenMMData* omm, double stepSizeInFs){
//    double dt = stepSizeInFs* OpenMM::PsPerFs;
//    double friction = generalParameters.frictionInPs;
//    double kBT    = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
//
//    omm->CustomLangevin = new OpenMM::CustomIntegrator(dt);
//
//
//    omm->CustomLangevin->addGlobalVariable("a", exp(-0.5*friction*dt));
//    omm->CustomLangevin->addGlobalVariable("b", sqrt(1-exp(-friction*dt)));
//    omm->CustomLangevin->addGlobalVariable("c", (1-exp(-0.5*friction*dt))/friction );
//    omm->CustomLangevin->addGlobalVariable("kT", kBT);
//
//    omm->CustomLangevin->addUpdateContextState();
//    omm->CustomLangevin->addComputePerDof("v", "a*v + c*f/m + b*sqrt(kT/m)*gaussian");
//    omm->CustomLangevin->addComputePerDof("x", "x + dt*v");
//    omm->CustomLangevin->addComputePerDof("v", "a*v + c*f/m + b*sqrt(kT/m)*gaussian");
}

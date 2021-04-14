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


void set_multithermos_dropNewton3(MyOpenMMData* omm,
                                  double stepSizeInFs,
                                  vector<OpenMM::CustomCompoundBondForce*> &DihedralForces,
                                  vector<OpenMM::CustomNonbondedForce*>    &WCAs,
                                  const MyAtomInfo  atoms[]){
    /**
                Leap-frog Langevin integrator. Reference:
                Jesús A. Izaguirre, Chris R. Sweet, and Vijay S. Pande. Multiscale dynamics of macromolecules using Normal Mode Langevin. Pacific Symposium on Biocomputing, 15:240–251, 2010.
     */
    
    DihedralForces[0]->setForceGroup(31);
    WCAs[0]->setForceGroup(30);
    
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionInPs;
    double kBT    = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
    double kBTmem = generalParameters.BoltzmannKJpermolkelvin*generalParameters.customtemperature;
    
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("a", exp(-0.5*friction*dt));
    omm->CustomIntegrator->addGlobalVariable("b", sqrt(1-exp(-friction*dt)));
    omm->CustomIntegrator->addGlobalVariable("c", (1-exp(-0.5*friction*dt))/friction );
    omm->CustomIntegrator->addGlobalVariable("kTa", kBT);
    omm->CustomIntegrator->addGlobalVariable("kTb", kBTmem);
    omm->CustomIntegrator->addPerDofVariable("stat", 0);
    
    
    omm->CustomIntegrator->addUpdateContextState();
    omm->CustomIntegrator->addComputePerDof("stat", "f31");
    
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*(f0)/m + b*( (1-delta(stat))*(sqrt(kTb/m) + delta(stat)*sqrt(kTa/m) )*gaussian)");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(stat)*c*(f30)/m");
//    omm->CustomIntegrator->addComputePerDof("v", "v + b*( (1-delta(stat))*(sqrt(kTb/m) + delta(stat)*sqrt(kTa/m) )*gaussian)");

//    omm->CustomIntegrator->addComputePerDof("v", "v + b*sqrt(kTa/m)*gaussian");

    omm->CustomIntegrator->addComputePerDof("x", "x + dt*v");
    
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*(f0)/m + b*( (1-delta(stat))*(sqrt(kTb/m) + delta(stat)*sqrt(kTa/m) )*gaussian)");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(stat)*c*(f30)/m");
//    omm->CustomIntegrator->addComputePerDof("v", "v + b*( (1-delta(stat))*(sqrt(kTb/m) + delta(stat)*sqrt(kTa/m) )*gaussian)");
    
}





void set_customLangevinforminimisation(MyOpenMMData* omm, double stepSizeInFs, double restraint){
    
    cout<<"Using the modified Langevin thermostat for minimisation with position restriction with tolerance "<<TFILE<<restraint<<TRESET<<endl;
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionInPs;
    double kBT    = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
    
    omm->LangevinMinimisation = new OpenMM::CustomIntegrator(dt);
    
    
//    omm->LangevinMinimisation->addGlobalVariable("a", exp(-0.5*friction*dt));
//    omm->LangevinMinimisation->addGlobalVariable("b", sqrt(1-exp(-friction*dt)));
//    omm->LangevinMinimisation->addGlobalVariable("c", (1-exp(-0.5*friction*dt))/friction );
//    omm->LangevinMinimisation->addGlobalVariable("kT", kBT);
    omm->LangevinMinimisation->addGlobalVariable("tol", restraint);
//    omm->LangevinMinimisation->addGlobalVariable("scale", 0.000001);
//    omm->LangevinMinimisation->addPerDofVariable("fmin", 1);
    
    
    omm->LangevinMinimisation->addUpdateContextState();
    

    
    
//    omm->LangevinMinimisation->addComputePerDof("v", "c*f/m + b*sqrt(kT/m)*gaussian");
//    omm->LangevinMinimisation->addComputePerDof("v", "f/m");
    
//    omm->LangevinMinimisation->addComputePerDof("fmin", "f");
//    omm->LangevinMinimisation->addComputePerDof("fmin", "sqrt(dot(f,f))");
//    omm->LangevinMinimisation->beginIfBlock("fmin > tol ");
//    omm->LangevinMinimisation->addComputePerDof("fmin", "scale*f");
//    omm->LangevinMinimisation->endBlock();

    
    
//    omm->LangevinMinimisation->addComputePerDof("v", " c*f/m + b*sqrt(kT/m)*gaussian");
//    omm->LangevinMinimisation->beginIfBlock("v > tol ");
//    omm->LangevinMinimisation->addComputePerDof("v", "tol");
//    omm->LangevinMinimisation->endBlock();
    omm->LangevinMinimisation->addComputePerDof("x", "x + max(-tol,min(0.5*dt*dt*f/m,tol))");

//    omm->LangevinMinimisation->addComputePerDof("v", "a*v + c*f/m + b*sqrt(kT/m)*gaussian");
    omm->LangevinMinimisation->addComputePerDof("v", "0");
    
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

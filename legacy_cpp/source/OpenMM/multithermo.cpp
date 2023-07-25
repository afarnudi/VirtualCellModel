#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"



void set_customLangevinforminimisation(MyOpenMMData* omm, double stepSizeInFs, double restraint){
    
    cout<<"Using the modified Langevin thermostat for minimisation with position restriction with tolerance "<<TFILE<<restraint<<TRESET<<endl;
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionIninvertPs;
    double kBT    = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    
    omm->CustomIntegrator->addGlobalVariable("tol", restraint);
    
    omm->CustomIntegrator->addUpdateContextState();

    omm->CustomIntegrator->addComputePerDof("x", "x + max(-tol,min(0.5*dt*dt*f/m,tol))");

    omm->CustomIntegrator->addComputePerDof("v", "0");
    
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



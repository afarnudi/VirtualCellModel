#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

using namespace std;

void set_LFLangevin_multithermos(MyOpenMMData* omm,
                                 double stepSizeInPs,
                                 double friction_invertPs,
                                 double kBT,
                                 double kBTmem,
                                 vector<OpenMM::CustomCompoundBondForce*> &DihedralForces,
                                 vector<OpenMM::CustomNonbondedForce*>    &WCAs
                                 ){
    /**
                Leap-frog Langevin integrator. Reference:
                Jesús A. Izaguirre, Chris R. Sweet, and Vijay S. Pande. Multiscale dynamics of macromolecules using Normal Mode Langevin. Pacific Symposium on Biocomputing, 15:240–251, 2010.
     */
    DihedralForces[0]->setForceGroup(31);
    WCAs[0]->setForceGroup(30);
    
    
    double dt = stepSizeInPs;

    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("a", exp(-0.5*friction_invertPs*dt));
    omm->CustomIntegrator->addGlobalVariable("b", sqrt(1-exp(-friction_invertPs*dt)));
    omm->CustomIntegrator->addGlobalVariable("c", (1-exp(-0.5*friction_invertPs*dt))/friction_invertPs );
    omm->CustomIntegrator->addGlobalVariable("sqkTa", sqrt(kBT));
    omm->CustomIntegrator->addGlobalVariable("sqkTb", sqrt(kBTmem));
    
    omm->CustomIntegrator->addPerDofVariable("f_n31", 0.0);
    
    //integrate all bonded forces
    omm->CustomIntegrator->addUpdateContextState();
    omm->CustomIntegrator->addComputePerDof("f_n31", "f31");
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f0/m + b*gaussian*( (1-delta(f_n31))*sqkTb + delta(f_n31)*sqkTa )/sqrt(m)");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f_n31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f30)/m");

    omm->CustomIntegrator->addComputePerDof("x", "x + dt*v");
    
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f0/m + b*gaussian*( (1-delta(f_n31))*sqkTb + delta(f_n31)*sqkTa )/sqrt(m)");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f30)/m");
    
}

void set_LFLangevin_multithermos_dropNewton3(MyOpenMMData* omm,
                                             double stepSizeInPs,
                                             double friction_invertPs,
                                             double kBT,
                                             double kBTmem,
                                             vector<OpenMM::CustomCompoundBondForce*> &DihedralForces,
                                             vector<OpenMM::CustomNonbondedForce*>    &WCAs
                                             ){
    /**
                Leap-frog Langevin integrator. Reference:
                Jesús A. Izaguirre, Chris R. Sweet, and Vijay S. Pande. Multiscale dynamics of macromolecules using Normal Mode Langevin. Pacific Symposium on Biocomputing, 15:240–251, 2010.
     */
    DihedralForces[0]->setForceGroup(31);
    WCAs[0]->setForceGroup(30);
    
    
    double dt = stepSizeInPs;

    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("a", exp(-0.5*friction_invertPs*dt));
    omm->CustomIntegrator->addGlobalVariable("b", sqrt(1-exp(-friction_invertPs*dt)));
    omm->CustomIntegrator->addGlobalVariable("c", (1-exp(-0.5*friction_invertPs*dt))/friction_invertPs );
    omm->CustomIntegrator->addGlobalVariable("sqkTa", sqrt(kBT));
    omm->CustomIntegrator->addGlobalVariable("sqkTb", sqrt(kBTmem));
    
    omm->CustomIntegrator->addPerDofVariable("f_n31", 0.0);
    
    //integrate all bonded forces
    omm->CustomIntegrator->addUpdateContextState();
    omm->CustomIntegrator->addComputePerDof("f_n31", "f31");
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f0/m + b*gaussian*( (1-delta(f_n31))*sqkTb + delta(f_n31)*sqkTa )/sqrt(m)");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f_n31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(f_n31)*c*(f30)/m");

    omm->CustomIntegrator->addComputePerDof("x", "x + dt*v");
    
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*f0/m + b*gaussian*( (1-delta(f_n31))*sqkTb + delta(f_n31)*sqkTa )/sqrt(m)");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(f_n31)*c*(f30)/m");
    
}

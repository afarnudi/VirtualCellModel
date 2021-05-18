#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

using namespace std;

void set_GJF(MyOpenMMData* omm,
             double stepSizeInPs,
             double friction_invertPs,
             double kBT){
    /**
     A simple and effective Verlet-type algorithm for simulating Langevin dynamics
     
     Niels Grønbech-Jensen  & Oded Farago
     Accepted author version posted online: 09 Jan 2013.Published online: 14 Feb 2013.
     DOI:10.1080/00268976.2012.760055
     */
    
   
    double dt = stepSizeInPs;
    double sigma  = sqrt(2*friction_invertPs*kBT*dt);
    double b_gjf = 1./(1+(dt*friction_invertPs)/2.);
    double a_gjf = (1-(dt*friction_invertPs)/2.)/(1+(dt*friction_invertPs)/2.);
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("g1", b_gjf*dt );
    omm->CustomIntegrator->addGlobalVariable("g2", b_gjf*dt*dt*0.5 );
    omm->CustomIntegrator->addGlobalVariable("g3", b_gjf*dt*sigma*0.5 );
    omm->CustomIntegrator->addGlobalVariable("a_gjf", (1-(dt*friction_invertPs)/2.)/(1+(dt*friction_invertPs)/2.) );
    omm->CustomIntegrator->addGlobalVariable("g4", dt*a_gjf*0.5 );
    omm->CustomIntegrator->addGlobalVariable("g5", b_gjf*sigma );
    
    omm->CustomIntegrator->addPerDofVariable("f_n", 0.0);
    omm->CustomIntegrator->addPerDofVariable("beta_n_1", 0.0);
    

    omm->CustomIntegrator->addUpdateContextState();
    
    omm->CustomIntegrator->addComputePerDof("beta_n_1", "gaussian");
    omm->CustomIntegrator->addComputePerDof("f_n", "f");
    
    omm->CustomIntegrator->addComputePerDof("x", "x + g1*v + g2*f_n/m + g3*beta_n_1/sqrt(m)");
    omm->CustomIntegrator->addComputePerDof("v", "a_gjf*v  + g4*f_n/m + g5*beta_n_1/sqrt(m) + dt*f/(2*m)");
    
}

void set_GJF_multithermos(MyOpenMMData* omm,
                          double stepSizeInPs,
                          double friction_invertPs,
                          double kBT,
                          double kBTmem,
                          vector<OpenMM::CustomCompoundBondForce*> &DihedralForces,
                          vector<OpenMM::CustomNonbondedForce*>    &WCAs){
    /**
     A simple and effective Verlet-type algorithm for simulating Langevin dynamics
     
     Niels Grønbech-Jensen  & Oded Farago
     Accepted author version posted online: 09 Jan 2013.Published online: 14 Feb 2013.
     DOI:10.1080/00268976.2012.760055
     */
    
    DihedralForces[0]->setForceGroup(31);
    WCAs[0]->setForceGroup(30);
    
    double dt = stepSizeInPs;
    
    double sigmaA  = sqrt(2*friction_invertPs*kBT*dt);
    double sigmaB  = sqrt(2*friction_invertPs*kBT*dt);
    double b_gjf = 1./(1+(dt*friction_invertPs)/2.);
    double a_gjf = (1-(dt*friction_invertPs)/2.)/(1+(dt*friction_invertPs)/2.);
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("g1", b_gjf*dt );
    omm->CustomIntegrator->addGlobalVariable("g2", b_gjf*dt*dt*0.5 );
    omm->CustomIntegrator->addGlobalVariable("g3A", b_gjf*dt*sigmaA*0.5 );
    omm->CustomIntegrator->addGlobalVariable("g3B", b_gjf*dt*sigmaB*0.5 );
    omm->CustomIntegrator->addGlobalVariable("a_gjf", (1-(dt*friction_invertPs)/2.)/(1+(dt*friction_invertPs)/2.) );
    omm->CustomIntegrator->addGlobalVariable("g4", dt*a_gjf*0.5 );
    omm->CustomIntegrator->addGlobalVariable("g5A", b_gjf*sigmaA );
    omm->CustomIntegrator->addGlobalVariable("g5B", b_gjf*sigmaB );
    
    omm->CustomIntegrator->addPerDofVariable("f_n0", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n30", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n31", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n1_0", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n1_30", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n1_31", 0.0);
    omm->CustomIntegrator->addPerDofVariable("beta_n_1", 0.0);
    

    omm->CustomIntegrator->addUpdateContextState();
    
    omm->CustomIntegrator->addComputePerDof("beta_n_1", "gaussian");
    omm->CustomIntegrator->addComputePerDof("f_n0", "f0");
    omm->CustomIntegrator->addComputePerDof("f_n30", "f30");
    omm->CustomIntegrator->addComputePerDof("f_n31", "f31");
    
    omm->CustomIntegrator->addComputePerDof("x", "x + g1*v + g2*(f_n0 + f_n30 + f_n31)/m + (g3A*delta(f_n31) +  (1-delta(f_n31))*g3B )*beta_n_1/sqrt(m)");
    
    omm->CustomIntegrator->addComputePerDof("f_n1_0", "f0");
    omm->CustomIntegrator->addComputePerDof("f_n1_30", "f30");
    omm->CustomIntegrator->addComputePerDof("f_n1_31", "f31");
    
    omm->CustomIntegrator->addComputePerDof("v", "a_gjf*v  + g4*(f_n0 + f_n30 + f_n31)/m + (g5A*delta(f_n31) +  (1-delta(f_n31))*g5B )*beta_n_1/sqrt(m) + dt*(f_n1_0 + f_n1_30 + f_n1_31)/(2*m)");
    
}

void set_GJF_multithermos_DropNewton3(MyOpenMMData* omm,
                                      double stepSizeInPs,
                                      double friction_invertPs,
                                      double kBT,
                                      double kBTmem,
                                      vector<OpenMM::CustomCompoundBondForce*> &DihedralForces,
                                      vector<OpenMM::CustomNonbondedForce*>    &WCAs){
    /**
     A simple and effective Verlet-type algorithm for simulating Langevin dynamics
     
     Niels Grønbech-Jensen  & Oded Farago
     Accepted author version posted online: 09 Jan 2013.Published online: 14 Feb 2013.
     DOI:10.1080/00268976.2012.760055
     */
    
    DihedralForces[0]->setForceGroup(31);
    WCAs[0]->setForceGroup(30);
    
    double dt = stepSizeInPs;
    
    double sigmaA  = sqrt(2*friction_invertPs*kBT*dt);
    double sigmaB  = sqrt(2*friction_invertPs*kBT*dt);
    double b_gjf = 1./(1+(dt*friction_invertPs)/2.);
    double a_gjf = (1-(dt*friction_invertPs)/2.)/(1+(dt*friction_invertPs)/2.);
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("g1", b_gjf*dt );
    omm->CustomIntegrator->addGlobalVariable("g2", b_gjf*dt*dt*0.5 );
    omm->CustomIntegrator->addGlobalVariable("g3A", b_gjf*dt*sigmaA*0.5 );
    omm->CustomIntegrator->addGlobalVariable("g3B", b_gjf*dt*sigmaB*0.5 );
    omm->CustomIntegrator->addGlobalVariable("a_gjf", (1-(dt*friction_invertPs)/2.)/(1+(dt*friction_invertPs)/2.) );
    omm->CustomIntegrator->addGlobalVariable("g4", dt*a_gjf*0.5 );
    omm->CustomIntegrator->addGlobalVariable("g5A", b_gjf*sigmaA );
    omm->CustomIntegrator->addGlobalVariable("g5B", b_gjf*sigmaB );
    
    omm->CustomIntegrator->addPerDofVariable("f_n0", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n30", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n31", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n1_0", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n1_30", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n1_31", 0.0);
    omm->CustomIntegrator->addPerDofVariable("beta_n_1", 0.0);
    

    omm->CustomIntegrator->addUpdateContextState();
    
    omm->CustomIntegrator->addComputePerDof("beta_n_1", "gaussian");
    omm->CustomIntegrator->addComputePerDof("f_n0", "f0");
    omm->CustomIntegrator->addComputePerDof("f_n30", "f30");
    omm->CustomIntegrator->addComputePerDof("f_n31", "f31");
    
    omm->CustomIntegrator->addComputePerDof("x", "x + g1*v + g2*(f_n0 + delta(f_n31)*f_n30 + f_n31)/m + (g3A*delta(f_n31) +  (1-delta(f_n31))*g3B )*beta_n_1/sqrt(m)");
    
    omm->CustomIntegrator->addComputePerDof("f_n1_0", "f0");
    omm->CustomIntegrator->addComputePerDof("f_n1_30", "f30");
    omm->CustomIntegrator->addComputePerDof("f_n1_31", "f31");
    
    omm->CustomIntegrator->addComputePerDof("v", "a_gjf*v  + g4*(f_n0 + delta(f_n31)*f_n30 + f_n31)/m + (g5A*delta(f_n31) +  (1-delta(f_n31))*g5B )*beta_n_1/sqrt(m) + dt*(f_n1_0 + delta(f_n31)*f_n1_30 + f_n1_31)/(2*m)");
    
}

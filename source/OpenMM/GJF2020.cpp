#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

using namespace std;

void set_GJF2020(MyOpenMMData* omm,
                 double stepSizeInPs,
                 double friction_invertPs,
                 double kBT,
                 string GJFcase){
    /**
     Defining velocities for accurate kinetic statistics in the GJF thermostat
     Niels Grønbech-Jensen and Oded Farago
     DOI: 10.1103/PhysRevE.101.022123
     */
    
    
       
    double dt = stepSizeInPs;
    double sigma_gjf = sqrt(2*friction_invertPs*kBT*dt);

    double b_gjf = 1./(1+(dt*friction_invertPs)/2.);
    double a_gjf = (1-(dt*friction_invertPs)/2.)/(1+(dt*friction_invertPs)/2.);

    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    if (GJFcase=="A") {
        omm->CustomIntegrator->addGlobalVariable("a1", sqrt(b_gjf) );
        omm->CustomIntegrator->addGlobalVariable("a2", sqrt(b_gjf)*sigma_gjf*0.5 );
        omm->CustomIntegrator->addGlobalVariable("a3", sqrt(b_gjf)*dt*0.5 );
        omm->CustomIntegrator->addGlobalVariable("a4", sqrt(b_gjf)*dt );
        omm->CustomIntegrator->addGlobalVariable("a5", a_gjf/sqrt(b_gjf) );
        omm->CustomIntegrator->addGlobalVariable("a6", sigma_gjf*0.5 );
    } else if (GJFcase=="B") {
        omm->CustomIntegrator->addGlobalVariable("a_gjf", a_gjf );
        omm->CustomIntegrator->addGlobalVariable("b1", b_gjf*dt );
        omm->CustomIntegrator->addGlobalVariable("b2", b_gjf*dt*sigma_gjf*0.5 );
        omm->CustomIntegrator->addGlobalVariable("b3", b_gjf*sigma_gjf );
    } else if (GJFcase=="C") {
        omm->CustomIntegrator->addGlobalVariable("b_gjf", b_gjf );
        omm->CustomIntegrator->addGlobalVariable("c1", sqrt( b_gjf*(1+b_gjf) )*sigma_gjf );
        omm->CustomIntegrator->addGlobalVariable("c2", b_gjf*dt*0.5 );
        omm->CustomIntegrator->addGlobalVariable("c3", ( b_gjf-sqrt( b_gjf*(1+b_gjf) ) )*dt*0.5*sigma_gjf );
        omm->CustomIntegrator->addGlobalVariable("c4", a_gjf/b_gjf );
        omm->CustomIntegrator->addGlobalVariable("c5", sigma_gjf*( 2*b_gjf*b_gjf - a_gjf*sqrt( b_gjf*(1+b_gjf) ) )/(2*b_gjf) );
    }
    
    omm->CustomIntegrator->addPerDofVariable("u_n", 0.0);
    omm->CustomIntegrator->addPerDofVariable("beta_n_1", 0.0);
    

    omm->CustomIntegrator->addUpdateContextState();
    
    omm->CustomIntegrator->addComputePerDof("beta_n_1", "gaussian");
    if (GJFcase=="A") {
        omm->CustomIntegrator->addComputePerDof("u_n", "a1*v + a2*beta_n_1/sqrt(m) + a3*f/m");
        omm->CustomIntegrator->addComputePerDof("x"  , "x + a4*u_n");
        omm->CustomIntegrator->addComputePerDof("v"  , "a5*u_n + a6*beta_n_1/sqrt(m) + 0.5*dt*f/m");
    } else if (GJFcase=="B") {
        omm->CustomIntegrator->addComputePerDof("u_n", "v + dt*f/(2*m)");
        omm->CustomIntegrator->addComputePerDof("x"  , "x + b1*u_n + b2*beta_n_1/sqrt(m)");
        omm->CustomIntegrator->addComputePerDof("v"  , "a_gjf*u_n + b3*beta_n_1/sqrt(m) + dt*f/(2*m)");
    } else if (GJFcase=="C") {
        omm->CustomIntegrator->addComputePerDof("u_n", "b_gjf*v + c1*beta_n_1 + c2*beta_n_1/sqrt(m)");
        omm->CustomIntegrator->addComputePerDof("x"  , "x + dt*u_n + c3*beta_n_1/sqrt(m) ");
        omm->CustomIntegrator->addComputePerDof("v"  , "c4*u_n + c5*beta_n_1/sqrt(m) + dt*f/(2*m)");
    }

    
}

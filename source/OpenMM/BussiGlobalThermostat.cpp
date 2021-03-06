#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

using namespace std;

void set_Bussi_Global_thermostat(MyOpenMMData* omm,
                                 double stepSizeInPs,
                                 double friction_invertPs,
                                 double kBT,
                                 int number_of_atoms,
                                 bool CMMotionRemover){
    /**
     Giovanni Bussi and Michele Parrinello.
     Stochastics thermostats : comparison of local and global schemes.
     http://dx.doi.org/10.1016/j.cpc.2008.01.006
     */
    
   
    double dt = stepSizeInPs;
    int dof = 3*number_of_atoms -int(CMMotionRemover)*3;
    double c = exp(-2.0 * dt / friction_invertPs);
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("c_global", c );
    omm->CustomIntegrator->addGlobalVariable("dofs_global", dof);
//    omm->CustomIntegrator->addGlobalVariable("mke", dof * 0.5 * kBT );
    omm->CustomIntegrator->addGlobalVariable("ke", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("alpha", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("sign_alpha", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("R_global", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("S_global", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("R2", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("sum_R2", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("global_b1", ( 1-c )* 0.5 * kBT     );
    omm->CustomIntegrator->addGlobalVariable("global_b2", 2.0 * sqrt( c*(1-c)*0.5 * kBT ) );
    omm->CustomIntegrator->addGlobalVariable("global_b3", sqrt( c/( (1-c)* 0.5 * kBT) ) );
    
    // global thermostat coefficients
    omm->CustomIntegrator->addComputeSum("ke", "0.5*m*v*v");
    //"S_global" is an approximate since we assume "dofs_global" goes to inf
    //for "dofs_global" > 1000 (at least) this is a good approximation
    omm->CustomIntegrator->addComputeGlobal("R_global", "gaussian");
    
    omm->CustomIntegrator->beginIfBlock("dofs_global>1000");
    omm->CustomIntegrator->addComputeGlobal("S_global", "sqrt(2*(dofs_global-1))*gaussian+(dofs_global-1)");
    omm->CustomIntegrator->endBlock();
    
    omm->CustomIntegrator->beginIfBlock("dofs_global<=1000");
    omm->CustomIntegrator->addComputePerDof("R2", "gaussian^2");
    omm->CustomIntegrator->addComputeSum("sum_R2", "R2");
    omm->CustomIntegrator->addComputeGlobal("S_global", "sum_R2-gaussian^2");
    omm->CustomIntegrator->endBlock();
    
    omm->CustomIntegrator->addComputeGlobal("alpha","sqrt(c_global + global_b1*(S_global+R_global*R_global)/ke + global_b2*R_global*sqrt(1.0/ke) )");
    omm->CustomIntegrator->addComputeGlobal("sign_alpha","2*step(R_global + global_b3*sqrt(ke) )-1");
    //momenta
    omm->CustomIntegrator->addComputePerDof("v", "sign_alpha*alpha*v");
    //velocity-Verlet
    omm->CustomIntegrator->addUpdateContextState();
    omm->CustomIntegrator->addComputePerDof("v", "v + 0.5*dt*f/m");
    omm->CustomIntegrator->addComputePerDof("x", "x+dt*v");
    omm->CustomIntegrator->addComputePerDof("v", "v + 0.5*dt*f/m");
}

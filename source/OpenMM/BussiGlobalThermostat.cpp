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

void set_Bussi_Global_thermostat_multithermos_DropNewton3(MyOpenMMData* omm,
                                                          double stepSizeInPs,
                                                          double friction_invertPs,
                                                          double kBT,
                                                          int number_of_atoms,
                                                          bool CMMotionRemover,
                                                          int number_of_mem_atoms,
                                                          double kBTmem,
                                                          vector<OpenMM::CustomCompoundBondForce*> &DihedralForces,
                                                          vector<OpenMM::CustomNonbondedForce*>    &WCAs){
    /**
     Giovanni Bussi and Michele Parrinello.
     Stochastics thermostats : comparison of local and global schemes.
     http://dx.doi.org/10.1016/j.cpc.2008.01.006
     */
    
   
    double dt = stepSizeInPs;
    int dofCh = 3* (number_of_atoms -number_of_mem_atoms -int(CMMotionRemover));
    int dofMem = 3* (number_of_mem_atoms -int(CMMotionRemover));
    double c = exp(-2.0 * dt / friction_invertPs);
    
    DihedralForces[0]->setForceGroup(31);
    WCAs[0]->setForceGroup(30);
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
    
    omm->CustomIntegrator->addGlobalVariable("c_global", c );
    
    omm->CustomIntegrator->addGlobalVariable("dofs_globalCh", dofCh);
    omm->CustomIntegrator->addGlobalVariable("dofs_globalMem", dofMem);

    omm->CustomIntegrator->addGlobalVariable("keCh", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("keMem", 0.0 );
    
    omm->CustomIntegrator->addGlobalVariable("alphaCh", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("alphaMem", 0.0 );
    
    omm->CustomIntegrator->addGlobalVariable("sign_alphaCh", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("sign_alphaMem", 0.0 );
    
    omm->CustomIntegrator->addGlobalVariable("R_globalCh", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("R_globalMem", 0.0 );
    
    omm->CustomIntegrator->addGlobalVariable("S_globalCh", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("S_globalMem", 0.0 );
    
    omm->CustomIntegrator->addGlobalVariable("R2Ch", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("R2Mem", 0.0 );
    
    omm->CustomIntegrator->addGlobalVariable("sum_R2Ch", 0.0 );
    omm->CustomIntegrator->addGlobalVariable("sum_R2Mem", 0.0 );
    
    omm->CustomIntegrator->addGlobalVariable("global_b1Ch", ( 1-c )* 0.5 * kBT     );
    omm->CustomIntegrator->addGlobalVariable("global_b1Mem", ( 1-c )* 0.5 * kBTmem     );
    
    omm->CustomIntegrator->addGlobalVariable("global_b2Ch", 2.0 * sqrt( c*(1-c)*0.5 * kBT ) );
    omm->CustomIntegrator->addGlobalVariable("global_b2Mem", 2.0 * sqrt( c*(1-c)*0.5 * kBTmem ) );
    
    omm->CustomIntegrator->addGlobalVariable("global_b3Ch", sqrt( c/( (1-c)* 0.5 * kBT) ) );
    omm->CustomIntegrator->addGlobalVariable("global_b3Mem", sqrt( c/( (1-c)* 0.5 * kBTmem) ) );
    
    omm->CustomIntegrator->addPerDofVariable("f_n31", 0.0);
    
    
    omm->CustomIntegrator->addComputePerDof("f_n31", "f31");
    // global thermostat coefficients
//    omm->CustomIntegrator->addComputeSum("ke", "0.5*m*v*v");
    omm->CustomIntegrator->addComputeSum("keCh", "0.5*m*v*v *delta(f_n31)");
    omm->CustomIntegrator->addComputeSum("keMem","0.5*m*v*v *(1-delta(f_n31))");
    //"S_global" is an approximate since we assume "dofs_global" goes to inf
    //for "dofs_global" > 1000 (at least) this is a good approximation
    omm->CustomIntegrator->addComputeGlobal("R_globalCh", "gaussian");
    omm->CustomIntegrator->addComputeGlobal("R_globalMem", "gaussian");
    
    omm->CustomIntegrator->beginIfBlock("dofs_globalCh>1000");
    omm->CustomIntegrator->addComputeGlobal("S_globalCh", "sqrt(2*(dofs_globalCh-1))*gaussian+(dofs_globalCh-1)*delta(f_n31)");
    omm->CustomIntegrator->endBlock();
    omm->CustomIntegrator->beginIfBlock("dofs_globalMem>1000");
    omm->CustomIntegrator->addComputeGlobal("S_globalMem", "sqrt(2*(dofs_globalMem-1))*gaussian+(dofs_globalMem-1)*(1-delta(f_n31))");
    omm->CustomIntegrator->endBlock();
    
    
    omm->CustomIntegrator->beginIfBlock("dofs_globalCh<=1000");
    omm->CustomIntegrator->addComputePerDof("R2Ch", "gaussian^2");
    omm->CustomIntegrator->addComputeSum("sum_R2Ch", "R2Ch*delta(f_n31)");
    omm->CustomIntegrator->addComputeGlobal("S_globalCh", "sum_R2Ch-gaussian^2");
    omm->CustomIntegrator->endBlock();
    omm->CustomIntegrator->beginIfBlock("dofs_globalMem<=1000");
    omm->CustomIntegrator->addComputePerDof("R2Mem", "gaussian^2");
    omm->CustomIntegrator->addComputeSum("sum_R2Mem", "R2Mem*(1-delta(f_n31))");
    omm->CustomIntegrator->addComputeGlobal("S_globalMem", "sum_R2Mem-gaussian^2");
    omm->CustomIntegrator->endBlock();
    
    
    omm->CustomIntegrator->addComputeGlobal("alphaCh","sqrt(c_global + global_b1Ch*(S_globalCh+R_globalCh*R_globalCh)/keCh + global_b2Ch*R_globalCh*sqrt(1.0/keCh) )");
    omm->CustomIntegrator->addComputeGlobal("alphaMem","sqrt(c_global + global_b1Mem*(S_globalMem+R_globalMem*R_globalMem)/keMem + global_b2Mem*R_globalMem*sqrt(1.0/keMem) )");
    
    omm->CustomIntegrator->addComputeGlobal("sign_alphaCh","2*step(R_globalCh + global_b3Ch*sqrt(keCh) )-1");
    omm->CustomIntegrator->addComputeGlobal("sign_alphaMem","2*step(R_globalMem + global_b3Mem*sqrt(keMem) )-1");
    //momenta
    omm->CustomIntegrator->addComputePerDof("v", "sign_alphaCh*alphaCh*v*delta(f_n31) + sign_alphaMem*alphaMem*v*(1-delta(f_n31))");
    //velocity-Verlet
    omm->CustomIntegrator->addUpdateContextState();
    
    
    
    omm->CustomIntegrator->addComputePerDof("v", "v + 0.5*dt*(f0+f_n31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + 0.5*dt*f30*(1-delta(f_n31))/m");
    
    omm->CustomIntegrator->addComputePerDof("x", "x+dt*v");
    
    omm->CustomIntegrator->addComputePerDof("v", "v + 0.5*dt*(f0+f_n31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + 0.5*dt*f30*(1-delta(f_n31))/m");
}

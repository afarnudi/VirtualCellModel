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


void set_multithermos_dropNewton3_Langevin(MyOpenMMData* omm,
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
    omm->CustomIntegrator->addPerDofVariable("f_n31", 0);
    
    
    omm->CustomIntegrator->addUpdateContextState();
    omm->CustomIntegrator->addComputePerDof("f_n31", "f31");
    
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*(f0)/m + b*( (1-delta(f_n31))*(sqrt(kTb/m) + delta(f_n31)*sqrt(kTa/m) )*gaussian)");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f_n31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(f_n31)*c*(f30)/m");

    omm->CustomIntegrator->addComputePerDof("x", "x + dt*v");
    
    omm->CustomIntegrator->addComputePerDof("v", "a*v + c*(f0)/m + b*( (1-delta(f_n31))*(sqrt(kTb/m) + delta(f_n31)*sqrt(kTa/m) )*gaussian)");
    omm->CustomIntegrator->addComputePerDof("v", "v + c*(f31)/m");
    omm->CustomIntegrator->addComputePerDof("v", "v + delta(f_n31)*c*(f30)/m");

    
}

void set_multithermos_dropNewton3_GJF(MyOpenMMData* omm,
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
    int    num_of_atoms = 0;
    for (int i=0; atoms[i].type!=-1; i++) {
        num_of_atoms++;
    }
    
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);

    omm->CustomIntegrator->addGlobalVariable("a_gjf", (1-(dt*friction)/2.)/(1+(dt*friction)/2.) );
    omm->CustomIntegrator->addGlobalVariable("b_gjf", 1./(1+(dt*friction)/2.) );
    omm->CustomIntegrator->addGlobalVariable("x_sigma_gjfa", sqrt(0.5* friction*kBT*dt *dt*dt) );
    omm->CustomIntegrator->addGlobalVariable("v_sigma_gjfa", sqrt(  2* friction*kBT*dt ) );
    omm->CustomIntegrator->addGlobalVariable("x_sigma_gjfb", sqrt(0.5* friction*kBTmem*dt *dt*dt) );
    omm->CustomIntegrator->addGlobalVariable("v_sigma_gjfb", sqrt(  2* friction*kBTmem*dt ) );
    
//    omm->CustomIntegrator->addPerDofVariable("x_n", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n0", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n30", 0.0);
    omm->CustomIntegrator->addPerDofVariable("f_n31", 0.0);
    omm->CustomIntegrator->addPerDofVariable("beta_n_1", 0.0);
    

    omm->CustomIntegrator->addUpdateContextState();
    
    
    omm->CustomIntegrator->addComputePerDof("beta_n_1", "gaussian");
    omm->CustomIntegrator->addComputePerDof("f_n0", "f0");
    omm->CustomIntegrator->addComputePerDof("f_n30", "f30");
    omm->CustomIntegrator->addComputePerDof("f_n31", "f31");
    
    omm->CustomIntegrator->addComputePerDof("x", "x + b_gjf*dt*v + b_gjf*dt*dt*(f_n0+f_n31+delta(f_n31)*f_n30)/(2*m) + b_gjf*((1-delta(f_n31))*x_sigma_gjfb + delta(f_n31)*x_sigma_gjfa )*beta_n_1/sqrt(m)");
    
    
    
    omm->CustomIntegrator->addComputePerDof("v", "a_gjf*v + (dt/(2*m))*(a_gjf*(f_n0+f_n31+delta(f_n31)*f_n30) + f0) + b_gjf*((1-delta(f_n31))*v_sigma_gjfb + delta(f_n31)*v_sigma_gjfa )*beta_n_1/sqrt(m) ");
    omm->CustomIntegrator->addComputePerDof("v", "v+ (dt/(2*m))*(f31)");
    omm->CustomIntegrator->addComputePerDof("v", "v+ (dt/(2*m))*(delta(f_n31)*f30)");

}

void set_multithermos_GJF(MyOpenMMData* omm,
                          double stepSizeInFs,
                          vector<OpenMM::CustomCompoundBondForce*> &DihedralForces,
                          vector<OpenMM::CustomNonbondedForce*>    &WCAs,
                          const MyAtomInfo  atoms[]){
    /**
                Leap-frog Langevin integrator. Reference:
                Jesús A. Izaguirre, Chris R. Sweet, and Vijay S. Pande. Multiscale dynamics of macromolecules using Normal Mode Langevin. Pacific Symposium on Biocomputing, 15:240–251, 2010.
     */
    
   
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionInPs;
    double kBT    = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
//    double kBTmem = generalParameters.BoltzmannKJpermolkelvin*generalParameters.customtemperature;
    int    num_of_atoms = 0;
    for (int i=0; atoms[i].type!=-1; i++) {
        num_of_atoms++;
    }
    
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);

    omm->CustomIntegrator->addGlobalVariable("a_gjf", (1-(dt*friction)/2.)/(1+(dt*friction)/2.) );
    omm->CustomIntegrator->addGlobalVariable("b_gjf", 1./(1+(dt*friction)/2.) );
    omm->CustomIntegrator->addGlobalVariable("x_sigma_gjf", sqrt(0.5* friction*kBT*dt *dt*dt) );
    omm->CustomIntegrator->addGlobalVariable("v_sigma_gjf", sqrt(  2* friction*kBT*dt ) );
    
    omm->CustomIntegrator->addPerDofVariable("f_n", 0.0);
    omm->CustomIntegrator->addPerDofVariable("beta_n_1", 0.0);
    

    omm->CustomIntegrator->addUpdateContextState();
    
    omm->CustomIntegrator->addComputePerDof("beta_n_1", "gaussian");
    omm->CustomIntegrator->addComputePerDof("f_n", "f");
    
    omm->CustomIntegrator->addComputePerDof("x", "x + b_gjf*dt*v + b_gjf*dt*dt*f_n/(2*m) + b_gjf*x_sigma_gjf*beta_n_1/sqrt(m)");
    omm->CustomIntegrator->addComputePerDof("v", "a_gjf*v + (dt/(2*m))*(a_gjf*f_n + f) + b_gjf*v_sigma_gjf*beta_n_1/sqrt(m) ");

}




void set_customLangevinforminimisation(MyOpenMMData* omm, double stepSizeInFs, double restraint){
    
    cout<<"Using the modified Langevin thermostat for minimisation with position restriction with tolerance "<<TFILE<<restraint<<TRESET<<endl;
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionInPs;
    double kBT    = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
    
    omm->LangevinMinimisation = new OpenMM::CustomIntegrator(dt);
    
    
    omm->LangevinMinimisation->addGlobalVariable("tol", restraint);
    
    omm->LangevinMinimisation->addUpdateContextState();

    omm->LangevinMinimisation->addComputePerDof("x", "x + max(-tol,min(0.5*dt*dt*f/m,tol))");

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

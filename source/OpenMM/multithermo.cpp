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
     A simple and effective Verlet-type algorithm for simulating Langevin dynamics
     
     Niels Grønbech-Jensen  & Oded Farago
     Accepted author version posted online: 09 Jan 2013.Published online: 14 Feb 2013.
     DOI:10.1080/00268976.2012.760055
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


void set_multithermos_GJF2020(MyOpenMMData* omm,
                              double stepSizeInFs,
                              vector<OpenMM::CustomCompoundBondForce*> &DihedralForces,
                              vector<OpenMM::CustomNonbondedForce*>    &WCAs,
                              const MyAtomInfo  atoms[],
                              string GJFcase){
    /**
     Defining velocities for accurate kinetic statistics in the GJF thermostat
     Niels Grønbech-Jensen and Oded Farago
     DOI: 10.1103/PhysRevE.101.022123
     */
    
    
       
    double dt = stepSizeInFs* OpenMM::PsPerFs;
    double friction = generalParameters.frictionInPs;
    double kBT    = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
//    double kBTmem = generalParameters.BoltzmannKJpermolkelvin*generalParameters.customtemperature;
    int    num_of_atoms = 0;
    for (int i=0; atoms[i].type!=-1; i++) {
        num_of_atoms++;
    }
    double gamma1, gamma5;
    double b_gjf = 1./(1+(dt*friction)/2.);
    double a_gjf = (1-(dt*friction)/2.)/(1+(dt*friction)/2.);
    double Gamma4, Gamma5;
    
    if (GJFcase=="A") {
        gamma1 = 1/b_gjf;
        gamma5 = -0.5;
        Gamma4 = 2*b_gjf;
        Gamma5 = 0;
    } else if (GJFcase=="B") {
        gamma1 = 1/b_gjf;
        gamma5 = -0.5;
        Gamma4 = 2*b_gjf;
        Gamma5 = 0;
    } else if (GJFcase=="C") {
        gamma1 = 1;
        gamma5 = -0.5*(b_gjf-sqrt(b_gjf*(b_gjf+1)));
        Gamma4 = 2*b_gjf*b_gjf - a_gjf*sqrt( b_gjf*(1+b_gjf) );
        Gamma5 = sqrt( b_gjf*(1+b_gjf) );
    }
    
//    double Gamma4 = b_gjf*gamma1 - 2*a_gjf*gamma5;
//    double Gamma5 = b_gjf*gamma1 + 2*gamma5;
    
    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);

    omm->CustomIntegrator->addGlobalVariable("a_gjf", a_gjf );
    omm->CustomIntegrator->addGlobalVariable("b_gjf", b_gjf );
    omm->CustomIntegrator->addGlobalVariable("sigma_gjf", sqrt(2*friction*kBT) );
//    omm->CustomIntegrator->addGlobalVariable("g1", gamma1 );
//    omm->CustomIntegrator->addGlobalVariable("g5", gamma5 );
//    omm->CustomIntegrator->addGlobalVariable("h4", Gamma4 );
//    omm->CustomIntegrator->addGlobalVariable("h5", Gamma5 );
    
    omm->CustomIntegrator->addPerDofVariable("u_n", 0.0);
    omm->CustomIntegrator->addPerDofVariable("beta_n_1", 0.0);
    

    omm->CustomIntegrator->addUpdateContextState();
    
    omm->CustomIntegrator->addComputePerDof("beta_n_1", "gaussian");
    if (GJFcase=="A") {
        omm->CustomIntegrator->addComputePerDof("u_n", "sqrt(b_gjf)*(v + (sigma_gjf*beta_n_1 + dt*f)/(2*m) )");
        omm->CustomIntegrator->addComputePerDof("x"  , "x + sqrt(b_gjf)*dt*u_n");
        omm->CustomIntegrator->addComputePerDof("v"  , "a_gjf*u_n/sqrt(b_gjf) + (sigma_gjf*beta_n_1 + dt*f)/(2*m)");
    } else if (GJFcase=="B") {
        omm->CustomIntegrator->addComputePerDof("u_n", "v + dt*f/(2*m)");
        omm->CustomIntegrator->addComputePerDof("x"  , "x + b_gjf*dt*u_n + b_gjf*dt*sigma_gjf*beta_n_1/(2*m)");
        omm->CustomIntegrator->addComputePerDof("v"  , "a_gjf*u_n + b_gjf*sigma_gjf*beta_n_1/m + dt*f/(2*m)");
    } else if (GJFcase=="C") {
        omm->CustomIntegrator->addComputePerDof("u_n", "b_gjf*( v + sqrt( (1+b_gjf)/b_gjf )*sigma_gjf*beta_n_1 + dt*sigma_gjf*beta_n_1/m )");
        omm->CustomIntegrator->addComputePerDof("x"  , "x + dt*u_n + (b_gjf-sqrt( b_gjf*(1+b_gjf) ) )*dt*beta_n_1*sigma_gjf/m ");
        omm->CustomIntegrator->addComputePerDof("v"  , "a_gjf*u_n/b_gjf + sigma_gjf*beta_n_1*(2*b_gjf - a_gjf*sqrt( (1+b_gjf)/b_gjf ) )/(2*m) + dt*f/(2*m)");
    }

}

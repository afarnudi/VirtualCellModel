#include "Membrane.h"
#include "General_functions.hpp"
#include "Global_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include <vector>
#include <stdlib.h>
#include <math.h>

void myStepWithOpenMM(MyOpenMMData* omm,
                      TimeDependantData* time_dependant_data,
                      MyAtomInfo atoms[],
                      int numSteps ,
                      int& total_step) {
    
    OpenMM::Vec3 z(0,0,0);
     std::vector<OpenMM::Vec3> a ;
     int num = omm->system->getNumParticles();
     
     for(int i=0; i<num ; i++)
     {
         a.push_back(z);
     }
    
    if((time_dependant_data->Kelvin_Voigt) || (time_dependant_data->HillForce))
    {
        for (int i=0; i<numSteps; i++)
        {
            if(time_dependant_data->Kelvin_Voigt)
            {
                if (total_step % time_dependant_data->Kelvin_stepnum ==0)
                {
                    Cheap_GetOpenMMState(omm,atoms);
                    time_dependant_data->Kelvin_dist_calc(atoms);
                    if(time_dependant_data->Kelvin_Voigt_distInNm.size()>1)
                    {
                        //mys_state_update
                        Kelvin_Voigt_update(omm,time_dependant_data);
                        time_dependant_data->Kelvin_Voigt_distInNm.erase(time_dependant_data->Kelvin_Voigt_distInNm.begin());
                    }
                }
            }
            
            if(time_dependant_data->HillForce)
            {
                //for test
                if((total_step<34000000) && (total_step>170000) && (total_step % 1000 ==1))
                //if (total_step % time_dependant_data->hill_stepnum ==0)
                {
                    Cheap_GetOpenMMState(omm,atoms);
                    time_dependant_data->hill_dist_calc(atoms);
                    time_dependant_data->COM_calculator(atoms);
                    
                    //time_dependant_data->COM_calculator(atoms);
                    cout<<"x_com = " <<time_dependant_data->COM[0]<<'\n';
                    
                    if(time_dependant_data->hill_distInNm.size()>1)
                    {
                        //mys_state_update
                        hill_update(omm,time_dependant_data,atoms);
                    time_dependant_data->hill_distInNm.erase(time_dependant_data->hill_distInNm.begin());
                    }
                }
            }

            
            if ( omm->integrator != NULL ) {
                
                omm->integrator->step(1);
                total_step++;
            } else {
                
                omm->Lintegrator->step(1);
                total_step++;
            }
            
            //relax membrane
              if((total_step<163000) && (total_step>88000) && (total_step % 18000 ==1) )
                {
                    omm->context->setVelocities(a);
               }
        }
        
    }
    
    
    else
    {
        if ( omm->integrator != NULL ) {
            
            omm->integrator->step(numSteps);
            //        omm->context->computeVirtualSites();
                    total_step += numSteps;
        } else {
            
            omm->Lintegrator->step(numSteps);
            //        omm->context->computeVirtualSites();
                    total_step += numSteps;
        }
        
    }
    
}

void myTerminateOpenMM(MyOpenMMData* omm,
                       TimeDependantData* time_dependant_data) {
    delete omm;
    delete time_dependant_data;
}

using OpenMM::Vec3;

void Kelvin_Voigt_update(MyOpenMMData* omm,
                         TimeDependantData* time_dependant_data)
{
    const int Num_Bonds = time_dependant_data->Kelvin_VoigtBond->getNumBonds();
    int atom1, atom2 ;
    double length, stiffness;
    
    for(int i=0; i<Num_Bonds ; i++)
    {
        time_dependant_data->Kelvin_VoigtBond->getBondParameters(i, atom1, atom2, length, stiffness);
        
        length = time_dependant_data->Kelvin_Voigt_initNominal_length_InNm[i] - (time_dependant_data->Kelvin_Voigt_distInNm[1][i] - time_dependant_data->Kelvin_Voigt_distInNm[0][i]) * (time_dependant_data->Kelvin_Voigt_damp[i] * OpenMM::FsPerPs / stiffness)/(time_dependant_data->Kelvin_stepnum * generalParameters.Step_Size_In_Fs) ;
        
        time_dependant_data->Kelvin_VoigtBond->setBondParameters(i, atom1, atom2, length, stiffness);
    }
    time_dependant_data->Kelvin_VoigtBond->updateParametersInContext(*omm->context);
    
}






void hill_update(MyOpenMMData* omm,
                         TimeDependantData* time_dependant_data,  MyAtomInfo atoms[] )
{
    for(int j=0; j<time_dependant_data->Hill_force.size() ; j++)
    {
    int Num_Bonds = time_dependant_data->Hill_force[j]->getNumBonds();
        
    int atom1, atom2 ;
    //double length, stiffness;
    std::vector<double> parameters;
        
    
    for(int i=0; i<Num_Bonds ; i++)
    {
        time_dependant_data->Hill_force[j]->getBondParameters(i, atom1, atom2, parameters);
        
        
        if( 0.5*(atoms[atom1].posInNm[0] + atoms[atom2].posInNm[0]) < time_dependant_data->COM[0] )
        {
            if(0.5*(atoms[atom1].posInNm[1] + atoms[atom2].posInNm[1]) < time_dependant_data->COM[1] )
            {
               // parameters[1] = 350 ;
                parameters[1] = parameters[1]+25;
                double a = parameters[1];
                parameters[1] = fmin(80, a);
            }
            else
            {
                
//                parameters[1] = (-1)*parameters[1]+40;
//                double a = parameters[1];
//                parameters[1] = (-1) * fmin(150, a);
                
                
                parameters[1] = parameters[1]+6;
                double a = parameters[1];
                parameters[1] = fmin(20, a);
                
                
            }
        }
        else
        {
            if(0.5*(atoms[atom1].posInNm[1] + atoms[atom2].posInNm[1]) < time_dependant_data->COM[1] )
            {
              //  parameters[1] = 50 ;
                parameters[1] = parameters[1]+6;
                double a = parameters[1];
                parameters[1] = fmin(20, a);
                
//                parameters[1] = parameters[1]+50;
//                               double a = parameters[1];
//                               parameters[1] = fmin(320, a);
            }
            else
            {
                
//                parameters[1] = (-1)*parameters[1]+25;
//                double a = parameters[1];
//                parameters[1] = (-1) * fmin(60, a);
                
                parameters[1] = parameters[1]+25;
                               double a = parameters[1];
                               parameters[1] = fmin(80, a);
                
            }
        }
        
        
        
//         if( 0.5*(atoms[atom1].posInNm[0] + atoms[atom2].posInNm[0]) < time_dependant_data->COM[0] )
//                {
//                    if(0.5*(atoms[atom1].posInNm[1] + atoms[atom2].posInNm[1]) < time_dependant_data->COM[1] )
//                    {
//                       //parameters[1] = -140 ;
//
//                                        parameters[1] = (-1)*parameters[1]+40;
//                                        double a = parameters[1];
//                                        parameters[1] = (-1) * fmin(230, a);
//
//                    }
//                    else
//                    {
//                        //parameters[1] = -70 ;
//                        parameters[1] = (-1)*parameters[1]+20;
//                        double a = parameters[1];
//                        parameters[1] = (-1) * fmin(100, a);
//
//
////                        parameters[1] = parameters[1]+4;
////                        double a = parameters[1];
////                        parameters[1] = fmin(15, a);
//
//
//                    }
//                }
//                else
//                {
//                    if(0.5*(atoms[atom1].posInNm[1] + atoms[atom2].posInNm[1]) < time_dependant_data->COM[1] )
//                    {
//                      // parameters[1] = -140 ;
//
//                        parameters[1] = (-1)*parameters[1]+40;
//                                                         double a = parameters[1];
//                                                         parameters[1] = (-1) * fmin(230, a);
//
//
//        //                parameters[1] = parameters[1]+50;
//        //                               double a = parameters[1];
//        //                               parameters[1] = fmin(320, a);
//                    }
//                    else
//                    {
//                       // parameters[1] = -70 ;
//        parameters[1] = (-1)*parameters[1]+20;
//                             double a = parameters[1];
//                             parameters[1] = (-1) * fmin(100, a);
//
////                        parameters[1] = parameters[1]+15;
////                                       double a = parameters[1];
////                                       parameters[1] = fmin(60, a);
//
//                    }
//                }
        
        
        
        
        //real hill model
        //parameters[1] = (time_dependant_data->hill_const_force[i]) /(1 + abs(( time_dependant_data->hill_distInAng[1][j][i] - time_dependant_data->hill_distInAng[0][j][i]) / (time_dependant_data->hill_stepnum * GenConst::Step_Size_In_Fs  * parameters[4] * OpenMM::PsPerFs))  );
        
        
        time_dependant_data->Hill_force[j]->setBondParameters(i, atom1, atom2, parameters);
    }
        
        //std::cout<<parameters[1]<<'\n';
    time_dependant_data->Hill_force[j]->updateParametersInContext(*omm->context);
    }
}


void customLangevinIntegrator(MyOpenMMData* omm, double stepSizeInFs){
//    double dt = stepSizeInFs* OpenMM::PsPerFs;
//    double friction = GenConst::frictionInPs;
//    double temperature = GenConst::temperature;
//    double kB = GenConst::BoltzmannKJpermolkelvin;
//    omm->CustomIntegrator = new OpenMM::CustomIntegrator(dt);
//
//    omm->CustomIntegrator->addGlobalVariable("a", exp(-friction*dt));
//    omm->CustomIntegrator->addGlobalVariable("b", sqrt(1-exp(-2*friction*dt)));
//    omm->CustomIntegrator->addGlobalVariable("c", (1-exp(-friction*dt))/friction);
//    omm->CustomIntegrator->addGlobalVariable("kT", kB*temperature);
//    omm->CustomIntegrator->addUpdateContextState();
//
//    omm->CustomIntegrator->addComputePerDof("v", "v*a + c*f/m + b*sqrt(kT/m)*gaussian");
//
//    omm->CustomIntegrator->addComputePerDof("x", "x + dt*v");
//    omm->CustomIntegrator->addConstrainVelocities();
//    omm->CustomIntegrator->addConstrainPositions();
    
    
    omm->CustomIntegrator->addPerDofVariable("x0", 0);
    omm->CustomIntegrator->addUpdateContextState();
    omm->CustomIntegrator->addComputePerDof("x0", "x");
    omm->CustomIntegrator->addComputePerDof("v", "v+dt*f/m");
    omm->CustomIntegrator->addComputePerDof("x", "x+dt*v");
    omm->CustomIntegrator->addConstrainPositions();
    omm->CustomIntegrator->addComputePerDof("v", "(x-x0)/dt");
}

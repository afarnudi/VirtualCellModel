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
                      int& numSteps ,
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
                    Cheap_GetOpenMMState(omm->context,atoms);
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
                    Cheap_GetOpenMMState(omm->context,atoms);
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
            
            
            if ( generalParameters.Integrator_type=="Verlet" ) {
                omm->VerletIntegrator->step(1);
                total_step++;
            } else if ( generalParameters.Integrator_type=="Brownian" ) {
                omm->BrownianIntegrator->step(1);
                total_step++;
            } else if ( generalParameters.Integrator_type=="Langevin" ) {
                omm->LangevinIntegrator->step(1);
                total_step++;
            } else if ( generalParameters.Integrator_type=="LangevinMiddle" ) {
                omm->LangevinMiddleIntegrator->step(1);
                total_step++;
            } else if ( generalParameters.Integrator_type=="CustomLangevinDropNewton3" ) {
                omm->CustomIntegrator->step(1);
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
        if (generalParameters.expSampling) {
            CalculateSilentSteps(numSteps, total_step, generalParameters.expSamplingExponent);
        }
        if ( generalParameters.Integrator_type=="Verlet" ) {
            omm->VerletIntegrator->step(numSteps);
        } else if ( generalParameters.Integrator_type=="Brownian" ) {
            omm->BrownianIntegrator->step(numSteps);
        } else if ( generalParameters.Integrator_type=="Langevin" ) {
            omm->LangevinIntegrator->step(numSteps);
        } else if ( generalParameters.Integrator_type=="LangevinMiddle" ) {
            omm->LangevinMiddleIntegrator->step(numSteps);
        } else {
            omm->CustomIntegrator->step(numSteps);
        }
        total_step+=numSteps;
        //        else if ( generalParameters.Integrator_type=="LFLangevinMulti-thermos" || generalParameters.Integrator_type=="LFLangevinMulti-thermosDropNewton3" || generalParameters.Integrator_type=="GJF" || generalParameters.Integrator_type=="GJF2013Multi-thermos" || generalParameters.Integrator_type=="GJF2013Multi-thermosDropNewton3" || generalParameters.Integrator_type=="GJF2020" || generalParameters.Integrator_type=="LangevinMinimise" || generalParameters.Integrator_type=="Bussi2008") {
        //            omm->CustomIntegrator->step(numSteps);
        //            total_step+=numSteps;
        //        }
        //        else if ( generalParameters.Integrator_type=="LangevinMinimise" ) {
        //            omm->LangevinMinimisation->step(numSteps);
        //            total_step+=numSteps;
        //        }
        
        
    }
    
}

void CalculateSilentSteps(int& numSteps ,
                          int total_step,
                          double exponent){
    //    cout<<total_step<<endl;
    if (total_step==0) {
        numSteps = int( exp(exponent)  );
    } else {
        int m = round(  log(total_step)/exponent );
        numSteps = int( pow(M_E, (m+1)*exponent) - total_step  );
    }
    
    int power=2;
    while (numSteps==0) {
        int m = round(  log(total_step)/exponent );
        numSteps = int( pow(M_E, (m+power)*exponent) - total_step  );
        power++;
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


void set_pbcvectors(OpenMM::System &system){
    std::vector<Vec3> pbcxyz;
    if (generalParameters.Periodic_condtion_status) {
        pbcxyz.resize(3);
        
        pbcxyz[0][0]=generalParameters.Simulation_box_length;
        pbcxyz[0][1]=0;
        pbcxyz[0][2]=0;
        
        pbcxyz[1][0]=0;
        pbcxyz[1][1]=generalParameters.Simulation_box_length;
        pbcxyz[1][2]=0;
        
        pbcxyz[2][0]=0;
        pbcxyz[2][1]=0;
        pbcxyz[2][2]=generalParameters.Simulation_box_length;
        
        system.setDefaultPeriodicBoxVectors(pbcxyz[0], pbcxyz[1], pbcxyz[2]);
        
        pbcxyz[0][0]=0;
        pbcxyz[0][1]=0;
        pbcxyz[0][2]=0;
        
        pbcxyz[1][0]=0;
        pbcxyz[1][1]=0;
        pbcxyz[1][2]=0;
        
        pbcxyz[2][0]=0;
        pbcxyz[2][1]=0;
        pbcxyz[2][2]=0;
        
        system.getDefaultPeriodicBoxVectors(pbcxyz[0], pbcxyz[1], pbcxyz[2]);
        cout<<"Periodic Boundry Condition "<<TON<<"On"<<TRESET<<endl;
        cout<<"Periodic vectors (X,Y,Z):\n";
        cout<<TGRAY;
        for (int i=0; i<3; i++) {
            cout<<i<<": "<<pbcxyz[i][0]<<"\t"<<pbcxyz[i][1]<<"\t"<<pbcxyz[i][2]<<"\n";
        }
        cout<<TRESET;
    } else {
        cout<<"Periodic Boundry Condition "<<TOFF<<"Off"<<TRESET<<endl;
    }
}


string generate_parameter_name(string    parameter_name,
                               int       set_1_index,
                               int       set_2_index,
                               string    set_1_name,
                               string    set_2_name){
    string parameter   = parameter_name + set_1_name + std::to_string(set_1_index) + set_2_name + std::to_string(set_2_index) ;
    return parameter;
}
string generate_parameter_name(string    parameter_name,
                               string    set_1_name,
                               string    set_2_name){
    string parameter   = parameter_name + set_1_name + set_2_name ;
    return parameter;
}

set<int> get_flat_set(vector<vector<set<int> > > vec_vec_set,
                      int                        set_1_index){
    set<int> flat_set;
    for (int i=0; i<vec_vec_set[set_1_index].size(); i++) {
        flat_set.insert(vec_vec_set[set_1_index][i].begin(),vec_vec_set[set_1_index][i].end());
    }
    return flat_set;
}

set<int> get_flat_set(vector<set<int> >  vec_vec_set,
                      int                set_1_index){
    return vec_vec_set[set_1_index];
}


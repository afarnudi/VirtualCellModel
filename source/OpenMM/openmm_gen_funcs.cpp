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
                    if(time_dependant_data->Kelvin_Voigt_distInAng.size()>1)
                    {
                        //mys_state_update
                        Kelvin_Voigt_update(omm,time_dependant_data);
                        time_dependant_data->Kelvin_Voigt_distInAng.erase(time_dependant_data->Kelvin_Voigt_distInAng.begin());
                    }
                }
            }
            
            if(time_dependant_data->HillForce)
            {
                //for test
                if((total_step<340000) && (total_step>170000) && (total_step % 1000 ==1))
                //if (total_step % time_dependant_data->hill_stepnum ==0)
                {
                    Cheap_GetOpenMMState(omm,atoms);
                    time_dependant_data->hill_dist_calc(atoms);
                    time_dependant_data->COM_calculator(atoms);
                    
                    //time_dependant_data->COM_calculator(atoms);
                    cout<<"x_com = " <<time_dependant_data->COM[0]<<'\n';
                    
                    if(time_dependant_data->hill_distInAng.size()>1)
                    {
                        //mys_state_update
                        hill_update(omm,time_dependant_data,atoms);
                    time_dependant_data->hill_distInAng.erase(time_dependant_data->hill_distInAng.begin());
                    }
                }
            }

            
            omm->integrator->step(1);
            total_step++;
            
            //relax membrane
              if((total_step<163000) && (total_step>88000) && (total_step % 18000 ==1) )
                {
                    omm->context->setVelocities(a);
               }
        }
        
    }
    
    
    else
    {
        omm->integrator->step(numSteps);
        total_step += numSteps;
    }
    
}

void myTerminateOpenMM(MyOpenMMData* omm,
                       TimeDependantData* time_dependant_data) {
    delete omm;
    delete time_dependant_data;
}

using OpenMM::Vec3;
void myGetOpenMMState(MyOpenMMData* omm,
                      double& timeInPs,
                      double& energyInKcal,
                      double& potential_energyInKcal,
                      MyAtomInfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    infoMask += OpenMM::State::Velocities;  // for kinetic energy (cheapm)
    if (GenConst::WantEnergy) {
        infoMask += OpenMM::State::Energy;     // for pot. energy (expensive)
    }
    if (GenConst::WantForce) {
        infoMask += OpenMM::State::Forces;
    }
    // Forces are also available (and cheap).

    const OpenMM::State state = omm->context->getState(infoMask);
    timeInPs = state.getTime(); // OpenMM time is in ps already
    
    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    const std::vector<Vec3>& velInNmperPs  = state.getVelocities();
    const std::vector<Vec3>& Forces        = state.getForces();
    for (int i=0; i < (int)positionsInNm.size(); ++i){
        for (int j=0; j < 3; ++j){
            atoms[i].posInAng[j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;
            atoms[i].velocityInAngperPs[j] = velInNmperPs[i][j]  * OpenMM::AngstromsPerNm;
            if (GenConst::WantForce) {
                atoms[i].force[j]    = Forces[i][j] * OpenMM::KcalPerKJ * OpenMM::NmPerAngstrom;
            }
        }
    }
    
    
    // If energy has been requested, obtain it and convert from kJ to kcal.
    energyInKcal = 0;

    if (GenConst::WantEnergy){
    energyInKcal = (state.getPotentialEnergy() + state.getKineticEnergy())
     * OpenMM::KcalPerKJ;
        potential_energyInKcal = (state.getPotentialEnergy())
        * OpenMM::KcalPerKJ;
    }

}

using OpenMM::Vec3;

void Cheap_GetOpenMMState(MyOpenMMData* omm,
                          MyAtomInfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    const OpenMM::State state = omm->context->getState(infoMask);
    
    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    
    for (int i=0; i < (int)positionsInNm.size(); ++i){
        for (int j=0; j < 3; ++j){
            atoms[i].posInAng[j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;
        }
    }
}

//                               PDB FILE WRITER
// Given state data, output a single frame (pdb "model") of the trajectory.
void myWritePDBFrame(int frameNum,
                     double timeInPs,
                     double energyInKcal,
                     double potential_energy,
                     const MyAtomInfo atoms[],
                     const Bonds bonds[],
                     std::string traj_name,
                     std::string force_name)
{
    int EndOfList=-1;
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"MODEL     %d\n", frameNum);
    fprintf(pFile,"REMARK 250 time=%.3f ps; energy=%.3f kcal/mole\n",
            timeInPs, energyInKcal);
    int index=0;
    string hist = atoms[0].pdb;
    if (atoms[0].class_label == "Chromatin") {
        hist.pop_back();
    }
    char chain[]={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
    
    //    double occ=1;
    for (int n=0; atoms[n].type != EndOfList; ++n){
        string new_label = atoms[n].pdb;
        if (atoms[n].class_label == "Chromatin") {
            new_label.pop_back();
        }
        if (atoms[n].class_label == "ECM") {
            new_label.pop_back();
        }

        if (hist != new_label) {
            index++;
            hist = new_label;
        }
        fprintf(pFile,"ATOM  %5d %4s ETH %c   %4.0f %8.3f%8.3f%8.3f%6.2f%6.1f          %c\n",
//        fprintf(pFile,"ATOM  %5d %4s ETH %c%4.0f %8.3f%8.3f%8.3f%6.2f%6.1f\n",
                n+1,
                atoms[n].pdb,
                chain[index],
                double(index),
                atoms[n].posInAng[0],
                atoms[n].posInAng[1],
                atoms[n].posInAng[2],
                atoms[n].stretching_energy,
                atoms[n].energy,
                atoms[n].symbol);
    }
    
    // visualize bonds in pdb file
    if(frameNum==1)
    {
    if (GenConst::write_bonds_to_PDB) {
        for (int n=0; bonds[n].type != EndOfList; ++n){
            if(bonds[n].atoms[0] < bonds[n].atoms[1])
            {
                fprintf(pFile, "CONECT%5d%5d\n",bonds[n].atoms[0]+1,bonds[n].atoms[1]+1);
            }
            if(bonds[n].atoms[0] > bonds[n].atoms[1])
            {
                fprintf(pFile, "CONECT%5d%5d\n",bonds[n].atoms[1]+1,bonds[n].atoms[0]+1);
            }
        }
    }
    }
    
    fprintf(pFile,"ENDMDL\n");
    fclose (pFile);




    // write forces and potential energy on Force.pdb
    FILE* sFile;
    sFile = fopen (force_name.c_str(),"a");
    
    fprintf(sFile,"MODEL     %d\n", frameNum);
    fprintf(sFile,"REMARK 250 time=%.3f ps; energy=%.3f kcal/mole\n",
            timeInPs, potential_energy);

    for (int n=0; atoms[n].type != EndOfList; ++n){
            string new_label = atoms[n].pdb;
            if (atoms[n].class_label == "Chromatin") {
                new_label.pop_back();
            }
            if (atoms[n].class_label == "ECM") {
                new_label.pop_back();
            }

            if (hist != new_label) {
                index++;
                hist = new_label;
            }
            fprintf(sFile,"ATOM  %5d %4s ETH %c   %4.0f %8.1f%8.1f%8.1f%6.2f%6.2f          %c\n",
    //        fprintf(pFile,"ATOM  %5d %4s ETH %c%4.0f %8.3f%8.3f%8.3f%6.2f%6.1f\n",
                    n+1,
                    atoms[n].pdb,
                    chain[index],
                    double(index),
                    atoms[n].force[0],
                    atoms[n].force[1],
                    atoms[n].force[2],
                    atoms[n].stretching_energy,
                    atoms[n].energy,
                    atoms[n].symbol);
        }
    
    fprintf(sFile,"ENDMDL\n");
    fclose (sFile);


}


void Kelvin_Voigt_update(MyOpenMMData* omm,
                         TimeDependantData* time_dependant_data)
{
    const int Num_Bonds = time_dependant_data->Kelvin_VoigtBond->getNumBonds();
    int atom1, atom2 ;
    double length, stiffness;
    
    for(int i=0; i<Num_Bonds ; i++)
    {
        time_dependant_data->Kelvin_VoigtBond->getBondParameters(i, atom1, atom2, length, stiffness);
        
        length = time_dependant_data->Kelvin_Voigt_initNominal_length_InNm[i] - (time_dependant_data->Kelvin_Voigt_distInAng[1][i] - time_dependant_data->Kelvin_Voigt_distInAng[0][i])*(OpenMM::NmPerAngstrom) * (time_dependant_data->Kelvin_Voigt_damp[i] / stiffness)/(time_dependant_data->Kelvin_stepnum * GenConst::Step_Size_In_Fs) ;
        
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
        
        double tot_forward_force = 0;
        double tot_backward_force = 0;
    
    for(int i=0; i<Num_Bonds ; i++)
    {
        time_dependant_data->Hill_force[j]->getBondParameters(i, atom1, atom2, parameters);
        
        
        if( 0.5*(atoms[atom1].posInAng[0] + atoms[atom2].posInAng[0]) < time_dependant_data->COM[0] )
        {
            if(0.5*(atoms[atom1].posInAng[1] + atoms[atom2].posInAng[1]) < time_dependant_data->COM[1] )
            {
               // parameters[1] = 350 ;
                parameters[1] = parameters[1]+50;
                double a = parameters[1];
                parameters[1] = fmin(300, a);
            }
            else
            {
                // parameters[1] = 350 ;
//                parameters[1] = parameters[1]-20;
//                double a = parameters[1];
//                parameters[1] = fmax(-150, a);
                
//                parameters[1] = (-1)*parameters[1]+40;
//                double a = parameters[1];
//                parameters[1] = (-1) * fmin(150, a);
                
                
                parameters[1] = parameters[1]+20;
                double a = parameters[1];
                parameters[1] = fmin(80, a);
                
                
            }
        }
        else
        {
            if(0.5*(atoms[atom1].posInAng[1] + atoms[atom2].posInAng[1]) < time_dependant_data->COM[1] )
            {
              //  parameters[1] = 50 ;
                parameters[1] = parameters[1]+20;
                double a = parameters[1];
                parameters[1] = fmin(80, a);
                
//                parameters[1] = parameters[1]+50;
//                               double a = parameters[1];
//                               parameters[1] = fmin(320, a);
            }
            else
            {
              // parameters[1] = 50 ;
//                parameters[1] = parameters[1]-25;
//                double a = parameters[1];
//                parameters[1] = fmax(-60, a);
                
//                parameters[1] = (-1)*parameters[1]+25;
//                double a = parameters[1];
//                parameters[1] = (-1) * fmin(60, a);
                
                parameters[1] = parameters[1]+50;
                               double a = parameters[1];
                               parameters[1] = fmin(300, a);
                
//                parameters[1] = parameters[1]+25;
//                                double a = parameters[1];
//                                parameters[1] = fmin(60, a);
            }
        }
        //parameters[1] = (time_dependant_data->hill_const_force[i]) /(1 + abs(( time_dependant_data->hill_distInAng[1][j][i] - time_dependant_data->hill_distInAng[0][j][i]) / (time_dependant_data->hill_stepnum * GenConst::Step_Size_In_Fs  * parameters[4] * OpenMM::PsPerFs))  );
        
       // cout<<"between x= "<<atoms[atom1].posInAng[0] <<" and x= "<< atoms[atom2].posInAng[0] <<" Force is "<<parameters[1]<<'\n';
        
//        if( (atoms[atom1].posInAng[0] + atoms[atom2].posInAng[0]) >0 )
//        {
//            tot_forward_force+=parameters[1];
//        }
//        else
//        {
//            tot_backward_force+=parameters[1];
//        }
        
//        if(i==1)
//        {
//            cout<<"Force= "<<parameters[1];
//        }
        
        time_dependant_data->Hill_force[j]->setBondParameters(i, atom1, atom2, parameters);
    }
        
//        cout<<"tot_forward_force = "<<tot_forward_force<<'\n';
//        cout<<"tot_backward_force = "<<tot_backward_force<<'\n';
        
        //std::cout<<parameters[1]<<'\n';
    time_dependant_data->Hill_force[j]->updateParametersInContext(*omm->context);
    }
}

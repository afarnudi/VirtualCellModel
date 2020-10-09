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
void myGetOpenMMState(MyOpenMMData* omm,
                      double& timeInPs,
                      double& energyInKJ,
                      double& potential_energyInKJ,
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
    const OpenMM::State state = omm->context->getState(infoMask,GenConst::Periodic_box);
    
    timeInPs = state.getTime(); // OpenMM time is in ps already
    
    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    const std::vector<Vec3>& velInNmperPs  = state.getVelocities();
    
    if (GenConst::Periodic_box) {
        Vec3 Lboxx, Lboxy, Lboxz;
        state.getPeriodicBoxVectors(Lboxx, Lboxy, Lboxz);
        
        
        GenConst::Lboxdims.clear();
        GenConst::Lboxdims.resize(3);
        for (int i=0; i<3; i++) {
            GenConst::Lboxdims[i].resize(3,0);
        }
        for (int j=0; j<3; j++) {
            GenConst::Lboxdims[0][j]=Lboxx[j];
            GenConst::Lboxdims[1][j]=Lboxy[j];
            GenConst::Lboxdims[2][j]=Lboxz[j];
        }
    }
    
    for (int i=0; i < (int)positionsInNm.size(); ++i){
        for (int j=0; j < 3; ++j){
            atoms[i].posInNm[j] = positionsInNm[i][j];
            atoms[i].velocityInNmperPs[j] = velInNmperPs[i][j];
            
        }
    }
    if (GenConst::WantForce) {
        const std::vector<Vec3>& Forces = state.getForces();
        for (int i=0; i < (int)positionsInNm.size(); ++i){
            for (int j=0; j < 3; ++j){
                atoms[i].posInNm[j] = positionsInNm[i][j];
                atoms[i].velocityInNmperPs[j] = velInNmperPs[i][j];
                if (GenConst::WantForce) {
                    atoms[i].force[j]    = Forces[i][j];
                }
            }
        }
    }
    
    // If energy has been requested, obtain it in kJ/mol.
    energyInKJ = 0;

    if (GenConst::WantEnergy){
        energyInKJ = state.getPotentialEnergy() + state.getKineticEnergy();
        potential_energyInKJ = state.getPotentialEnergy();
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
            atoms[i].posInNm[j] = positionsInNm[i][j];
        }
    }
}

//                               PDB FILE WRITER
// Given state data, output a single frame (pdb "model") of the trajectory.
void myWritePDBFrame(int frameNum,
                     double timeInPs,
                     double energyInKJ,
                     double potential_energyInKJ,
                     const MyAtomInfo atoms[],
                     const Bonds bonds[])
{
    int EndOfList=-1;
    
    string traj_name= GenConst::trajectory_file_name+".pdb";
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"MODEL     %d\n", frameNum);
    fprintf(pFile,"REMARK 250 time=%.3f ps; energy=%6.6f potential energy=%.3f KJ/mole\n",
            timeInPs,
            energyInKJ,
            potential_energyInKJ);
//    cout<<endl;
//    cout<<"kbend*angle "<<potential_energyInKJ<<endl;
//    cout<<"thetha "<<potential_energyInKJ*180/M_PI<<endl;
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
//        fprintf(pFile,"ATOM  %5d %4s ETH %c   %4.0f %8.3f%8.3f%8.3f%6.2f%6.1f          %c\n",
        fprintf(pFile,"ATOM  %5d %4s ETH %c%4.0f    %8.3f%8.3f%8.3f%6.2f%6.1f\n",
                n+1,
                atoms[n].pdb,
                chain[index],
                double(index),
                atoms[n].posInNm[0],
                atoms[n].posInNm[1],
                atoms[n].posInNm[2],
                atoms[n].stretching_energy,
                atoms[n].energyInKJ);//,
//                atoms[n].symbol);
    }
    
    // visualize bonds in pdb file
    if(frameNum==0 && GenConst::write_bonds_to_PDB)
    {
    
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
    
    fprintf(pFile,"ENDMDL\n");
    fclose (pFile);
}

void writeXYZFrame  (int atom_count,
                     const MyAtomInfo atoms[],
                     std::string traj_name)
{
    int EndOfList=-1;
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"%d\ncomment\n", atom_count);

    for (int n=0; atoms[n].type != EndOfList; ++n){
        fprintf(pFile,"%4s\t%8.3f\t%8.3f\t%8.3f\n",
                atoms[n].pdb,
                atoms[n].posInNm[0],
                atoms[n].posInNm[1],
                atoms[n].posInNm[2]);
    }
    fclose (pFile);
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
        
        length = time_dependant_data->Kelvin_Voigt_initNominal_length_InNm[i] - (time_dependant_data->Kelvin_Voigt_distInNm[1][i] - time_dependant_data->Kelvin_Voigt_distInNm[0][i]) * (time_dependant_data->Kelvin_Voigt_damp[i] * OpenMM::FsPerPs / stiffness)/(time_dependant_data->Kelvin_stepnum * GenConst::Step_Size_In_Fs) ;
        
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





void myGetOpenMMState2(MyOpenMMData* omm,
                      double& timeInPs,
                      double& energyInKJ,
                      double& potential_energyInKJ,
                      vector<Membrane> mems,
                      MyAtomInfo atoms[])
{
    
    int infoMask = 0;
    infoMask += OpenMM::State::Forces;
    
    const OpenMM::State state1 = omm->context->getState(infoMask,GenConst::Periodic_box, 0b1000'0000);
    timeInPs = state1.getTime(); // OpenMM time is in ps already
    const std::vector<Vec3>& Forces1 = state1.getForces();
    
    mems[0].calculate_surface_area_with_voronoi();
    
    double pressure=0;
    for (int n=0; n< mems[0].get_num_of_nodes(); ++n){
        double Force[3]={Forces1[n][0],Forces1[n][1],Forces1[n][2]};
        
        
        vector<double>  norm = mems[0].get_node_normal_vec(n);
        double Normal[3]={norm[0],norm[1],norm[2]};
        
        
        
        double normalForce = abs( innerproduct(Force, Normal) );
        pressure += normalForce/mems[0].get_node_voronoi_area(n);
        
        
    }
    string traj_name= GenConst::trajectory_file_name+"ch0mem0_pressure.txt";
    std::ofstream write_pressure;
    write_pressure.open(traj_name.c_str(),std::ios_base::app);
    write_pressure<<timeInPs<<"\t"<<pressure<<endl;
}

#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

void myStepWithOpenMM(MyOpenMMData* omm, MyAtomInfo atoms[], int numSteps) {
    if((omm->Kelvin_Voigt || omm->Custom_Kelvin_Voigt))
    {
        for (int i=0; i<numSteps; i++)
        {
            Cheap_GetOpenMMState(omm,atoms);
    
                omm->Kelvin_Voigt_distInAng.push_back(dist_calc(omm, atoms, 0));
                if(omm->Kelvin_Voigt_distInAng.size()>1)
                {
                    //mys_state_update
                    my_system_update(omm, 0, omm->Kelvin_Voigt);
                    omm->Kelvin_Voigt_distInAng.erase(omm->Kelvin_Voigt_distInAng.begin());
                }
            
            
                omm->Custom_Kelvin_Voigt_distInAng.push_back(dist_calc(omm, atoms, 1));
                if(omm->Custom_Kelvin_Voigt_distInAng.size()>1)
                {
                    //my_system_update
                    my_system_update(omm, 1, omm->Custom_Kelvin_Voigt);
                    omm->Custom_Kelvin_Voigt_distInAng.erase(omm->Custom_Kelvin_Voigt_distInAng.begin());
                }
            
            omm->integrator->step(1);
            
        }
            
    }
    
    else
    {
        omm->integrator->step(numSteps);
    }
    
}

void myTerminateOpenMM(MyOpenMMData* omm) {
    delete omm;
}

using OpenMM::Vec3;
void myGetOpenMMState(MyOpenMMData* omm,
                      bool wantEnergy,
                      bool wantForce,
                      double& timeInPs,
                      double& energyInKcal,
                      MyAtomInfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    infoMask += OpenMM::State::Velocities;  // for kinetic energy (cheap)
    if (wantEnergy) {
        infoMask += OpenMM::State::Energy;     // for pot. energy (expensive)
    }
    if (wantForce) {
        infoMask += OpenMM::State::Forces;
    }
    // Forces are also available (and cheap).
    
    const OpenMM::State state = omm->context->getState(infoMask);
    timeInPs = state.getTime(); // OpenMM time is in ps already
    
    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    const std::vector<Vec3>& velInAngperPs  = state.getVelocities();
    const std::vector<Vec3>& Forces        = state.getForces();
    for (int i=0; i < (int)positionsInNm.size(); ++i){
        for (int j=0; j < 3; ++j){
            atoms[i].posInAng[j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;
            atoms[i].velocityInAngperPs[j] = velInAngperPs[i][j];
            if (wantForce) {
                atoms[i].force[j]    = Forces[i][j] * OpenMM::KcalPerKJ * OpenMM::NmPerAngstrom;
            }
        }
    }
    // If energy has been requested, obtain it and convert from kJ to kcal.
    energyInKcal = 0;
    if (wantEnergy)
        energyInKcal = (state.getPotentialEnergy() + state.getKineticEnergy())
        * OpenMM::KcalPerKJ;
}


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
                     bool wantforce,
                     double timeInPs,
                     double energyInKcal,
                     const MyAtomInfo atoms[],
                     std::string traj_name)
{
    int EndOfList=-1;
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"MODEL     %d\n", frameNum);
    fprintf(pFile,"REMARK 250 time=%.3f ps; energy=%.3f kcal/mole\n",
            timeInPs, energyInKcal);
    int index=1;
    string hist = atoms[0].pdb;
    
    char chain[]={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
    
    double occ=1;
    for (int n=0; atoms[n].type != EndOfList; ++n){
        if (hist != atoms[n].pdb) {
            index++;
            hist = atoms[n].pdb;
        }
        fprintf(pFile,"ATOM  %5d %4s ETH %c   %4.0f %8.3f%8.3f%8.3f%6.2f%6.1f          %c\n",
                n+1,
                atoms[n].pdb,
                chain[index],
                double(index),
                atoms[n].posInAng[0],
                atoms[n].posInAng[1],
                atoms[n].posInAng[2],
                occ,
                atoms[n].energy,
                atoms[n].symbol);
    }
    fprintf(pFile,"ENDMDL\n");
    fclose (pFile);
}


void calc_energy(vector<Membrane>    mem,
                MyAtomInfo          atoms[]){
    
    int node_A, node_B, node_C, node_D;
    
    double points[3][3];
    
    int mem_count=0;
    
    for (int i=0; i<mem.size(); i++) {
        for (int j=0; j<mem[i].get_num_of_nodes(); j++) {
            atoms[j+mem_count].energy = 0;
        }
        for (int j=0; j<mem[i].Triangle_Pair_Nodes.size(); j++) {
            
            double A, B, C, E, F, G;
            
            node_C=mem[i].Triangle_Pair_Nodes[j][0];
            node_A=mem[i].Triangle_Pair_Nodes[j][1];
            node_B=mem[i].Triangle_Pair_Nodes[j][2];
            node_D=mem[i].Triangle_Pair_Nodes[j][3];
 
            for (int k=0; k<3; k++) {
                points[0][k]=atoms[node_A+mem_count].posInAng[k];
                points[1][k]=atoms[node_B+mem_count].posInAng[k];
                points[2][k]=atoms[node_C+mem_count].posInAng[k];
            }
            calc_surface_coefficeints(points, A, B, C);
            
            for (int k=0; k<3; k++) {
                points[2][k]=atoms[node_D+mem_count].posInAng[k];
            }
            calc_surface_coefficeints(points, E, F, G);
            double cosine=0;
            double denominator = sqrt(A*A + B*B + C*C) * sqrt(E*E + F*F + G*G);
            if (denominator > 0.001) {
                cosine = ( A*E + B*F + C*G )/denominator;
            }
            
            
            //shift the value between 0 and 2 only for representation purposes.
            cosine += 1;
            double scale=5;
            atoms[node_A+mem_count].energy += scale*cosine/4;
            atoms[node_B+mem_count].energy += scale*cosine/4;
            atoms[node_C+mem_count].energy += scale*cosine/4;
            atoms[node_D+mem_count].energy += scale*cosine/4;
            
        }
        mem_count += mem[i].get_num_of_nodes();
    }
}

void calc_energy_2(vector<Membrane>    mem,
                   MyAtomInfo          atoms[]){
    
    int node_A, node_B, node_C, node_D;
    
    double points[3][3];
    
    int mem_count=0;
    
    for (int i=0; i<mem.size(); i++) {
        for (int j=0; j<mem[i].get_num_of_nodes(); j++) {
            atoms[j+mem_count].energy = 0;
        }
        for (int j=0; j<mem[i].Triangle_Pair_Nodes.size(); j++) {
            
            double AB[3], BA[3], AC[3], BD[3], ABxAC[3], BAxBD[3];
            
            node_C=mem[i].Triangle_Pair_Nodes[j][0];
            node_A=mem[i].Triangle_Pair_Nodes[j][1];
            node_B=mem[i].Triangle_Pair_Nodes[j][2];
            node_D=mem[i].Triangle_Pair_Nodes[j][3];
            
            for (int k=0; k<3; k++) {
                AB[k]=atoms[node_B+mem_count].posInAng[k]-atoms[node_A+mem_count].posInAng[k];
                BA[k]=-AB[k];
                AC[k]=atoms[node_C+mem_count].posInAng[k]-atoms[node_A+mem_count].posInAng[k];
                BD[k]=atoms[node_D+mem_count].posInAng[k]-atoms[node_B+mem_count].posInAng[k];
            }

            crossvector(ABxAC, AB, AC);
            crossvector(BAxBD, BA, BD);
            
            double dot =innerproduct(ABxAC, BAxBD), lABxACl=vector_length(ABxAC), lBAxBDl=vector_length(BAxBD);
            
            double cosine=(dot)/(lABxACl*lBAxBDl);
            double sine = (1.00001 - cosine)*0.5;
            
            sine = sqrt(sine);
            
            double scale=5;
            atoms[node_A+mem_count].energy += scale*sine/4;
            atoms[node_B+mem_count].energy += scale*sine/4;
            atoms[node_C+mem_count].energy += scale*sine/4;
            atoms[node_D+mem_count].energy += scale*sine/4;
            
        }
        mem_count += mem[i].get_num_of_nodes();
    }
}


void my_system_update(MyOpenMMData* omm,
                      int bondtype,
                      bool WantUpdate)
{
    if (WantUpdate)
    {
        switch (bondtype) {
            //Kelvin-Voigt
            case 0:
            {
                const int Num_Bonds = omm->Kelvin_VoigtBond->getNumBonds();
                int atom1, atom2 ;
                double length, stiffness;
                
                for(int i=0; i<Num_Bonds ; ++i)
                {
                    omm->Kelvin_VoigtBond->getBondParameters(i, atom1, atom2, length, stiffness);
                    
                    length = omm->Kelvin_Voigt_initNominal_length_InNm[i] - (omm->Kelvin_Voigt_distInAng[1][i] - omm->Kelvin_Voigt_distInAng[0][i])*(OpenMM::NmPerAngstrom) * (omm->damp * OpenMM::FsPerPs / stiffness)/(GenConst::Step_Size_In_Fs) ;
                    
                    omm->Kelvin_VoigtBond->setBondParameters(i, atom1, atom2, length, stiffness);
                }
                omm->Kelvin_VoigtBond->updateParametersInContext(*omm->context);

            }
            
                break;
                
                //Custom Kelvin-Voigt
            case 1:
            {
                for(int j=0; j< omm->Custom_Kelvin_VoigtBond.size(); j++)
                {
                    
                    const int Num_Bonds = omm->Custom_Kelvin_VoigtBond[j]->getNumBonds();
                    int atom1, atom2 ;
                    // parameters{nominal length, stiffness, damp}
                    std::vector<double> parameters;
                    
                    for(int i=0; i<Num_Bonds ; ++i)
                    {
                        omm->Custom_Kelvin_VoigtBond[j]->getBondParameters(i, atom1, atom2, parameters);
                        
                        parameters[0] = omm->Custom_Kelvin_Voigt_initNominal_length_InNm[i] - (omm->Custom_Kelvin_Voigt_distInAng[1][i] - omm->Custom_Kelvin_Voigt_distInAng[0][i])*(OpenMM::NmPerAngstrom) * (parameters[2] * OpenMM::FsPerPs / parameters[1])/(GenConst::Step_Size_In_Fs) ;
                        
                        omm->Custom_Kelvin_VoigtBond[j]->setBondParameters(i, atom1, atom2, parameters);
                    }
                    omm->Custom_Kelvin_VoigtBond[j]->updateParametersInContext(*omm->context);
                    
                }
                
            }
                
                break;
                
            default:
                break;
        }
        
    }
}


// calculate distance InAng
std::vector<double> dist_calc(MyOpenMMData* omm,
                              MyAtomInfo atoms[],
                              int bondtype)
{
    std::vector<double> distInAng ;
    
    switch (bondtype) {
        //Kelvin-Voigt
        case 0:
        {
            if(omm->Kelvin_Voigt)
            {
                const int Num_Bonds = omm->Kelvin_VoigtBond->getNumBonds();
                int atom1, atom2 ;
                double length, stiffness , dist;
                
                for(int i=0; i<Num_Bonds; i++)
                {
                    omm->Kelvin_VoigtBond->getBondParameters(i, atom1, atom2, length, stiffness);
                    dist =sqrt ( ( atoms[atom1].posInAng[0] - atoms[atom2].posInAng[0]) * ( atoms[atom1].posInAng[0] - atoms[atom2].posInAng[0]) + ( atoms[atom1].posInAng[1] - atoms[atom2].posInAng[1]) * ( atoms[atom1].posInAng[1] - atoms[atom2].posInAng[1]) + ( atoms[atom1].posInAng[2] - atoms[atom2].posInAng[2]) * ( atoms[atom1].posInAng[2] - atoms[atom2].posInAng[2]) ) ;
                    
                    distInAng.push_back(dist);
                }
            }
        }
            break;
            
            //Custom Kelvin-Voigt
        case 1:
        {
            if(omm->Custom_Kelvin_Voigt)
            {
                for(int j=0; j< omm->Custom_Kelvin_VoigtBond.size(); j++)
                {
                    const int Num_Bonds = omm->Custom_Kelvin_VoigtBond[j]->getNumBonds();
                    int atom1, atom2 ;
                    double dist;
                    // parameters{nominal length, stiffness, damp}
                    std::vector<double> parameters;
                    
                    for(int i=0; i<Num_Bonds; i++)
                    {
                        omm->Custom_Kelvin_VoigtBond[j]->getBondParameters(i, atom1, atom2, parameters);
                        dist =sqrt ( ( atoms[atom1].posInAng[0] - atoms[atom2].posInAng[0]) * ( atoms[atom1].posInAng[0] - atoms[atom2].posInAng[0]) + ( atoms[atom1].posInAng[1] - atoms[atom2].posInAng[1]) * ( atoms[atom1].posInAng[1] - atoms[atom2].posInAng[1]) + ( atoms[atom1].posInAng[2] - atoms[atom2].posInAng[2]) * ( atoms[atom1].posInAng[2] - atoms[atom2].posInAng[2]) ) ;
                        
                        distInAng.push_back(dist);
                    }
                }
            }
        }
            break;
            
        default:
            break;
    }
    
    
    return distInAng;
}
// Nominal length in nm
std::vector<double> Nominal_length_calc(MyOpenMMData* omm, int bondtype)
{
    std::vector<double> Nominal_Length_InNm ;
    
    switch (bondtype) {
        //Kelvin-Voigt
        case 0:
        {
            const int Num_Bonds = omm->Kelvin_VoigtBond->getNumBonds();
            int atom1, atom2 ;
            double length, stiffness;
            
            for(int i=0; i<Num_Bonds ; ++i)
            {
                omm->Kelvin_VoigtBond->getBondParameters(i, atom1, atom2, length, stiffness);
                Nominal_Length_InNm.push_back(length);
            }
            
        }
            
            break;
            
        //Custom Kelvin-Voigt
        case 1:
        {
            for(int j=0; j< omm->Custom_Kelvin_VoigtBond.size(); j++)
            {
                const int Num_Bonds = omm->Custom_Kelvin_VoigtBond[j]->getNumBonds();
                int atom1, atom2 ;
                // parameters{nominal length, stiffness, damp}
                std::vector<double> parameters;
                
                for(int i=0; i<Num_Bonds; i++)
                {
                    omm->Custom_Kelvin_VoigtBond[j]->getBondParameters(i, atom1, atom2, parameters);
                    Nominal_Length_InNm.push_back(parameters[0]);
                }
            }
        }
            break;
            
        default:
            break;
    }
    
    return Nominal_Length_InNm;
}

#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

double Membrane::calculating_the_bond_energy(int index, bool initial_or_final, MyAtomInfo  atoms[],int number_of_privious_mem_nodes){
    vector<double> bond;
    double bond_length;
    double bond_energy=0;
    int A,B;
    bond.resize(3);
    //In Triangle_Pair_Nodes nodes of the common bond are sorted in a way that locate between the uncommon nodes of triangles (indices 1 and 2) 
    if (initial_or_final==0){// initial, common bond 
        A=Triangle_Pair_Nodes[index][1];
        B=Triangle_Pair_Nodes[index][2];
        }
    else{
        A=Triangle_Pair_Nodes[index][0];
        B=Triangle_Pair_Nodes[index][3];
    }
    for(int i=0; i<3; i++){    
    bond[i]= (atoms[A+number_of_privious_mem_nodes].posInAng[i] -atoms[B+number_of_privious_mem_nodes].posInAng[i]) *OpenMM::NmPerAngstrom;
    //cout<<"Atom A  "<<atoms[A+number_of_privious_mem_nodes].posInAng[i]<<endl;
    //cout<<"Atom B  "<<atoms[B+number_of_privious_mem_nodes].posInAng[i]<<endl;
    }
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);
    
    //cout<<"bond between nodes  "<<A<<"  and   "<<B<<"  lenght   "<<bond_length<<endl;
    switch (spring_model){
        case 1: //Fene
        {
            cout<<"FEne is under construction"<<endl;
            
        }
        case 2: //harmonic
        {
            bond_energy= 0.5*Spring_coefficient
                                      //  * OpenMM::KJPerKcal
                                        * OpenMM::AngstromsPerNm
                                        * OpenMM::AngstromsPerNm
                                        *(bond_length-(Average_node_pair_length* OpenMM::NmPerAngstrom))
                                        *(bond_length-(Average_node_pair_length* OpenMM::NmPerAngstrom));

        }
        case 5: //realharmonic :)) (x4harmonic is the name but its potential is really 1/2*k*x^2)
        {
            bond_energy= 0.5*Spring_coefficient
                                      //  * OpenMM::KJPerKcal
                                        * OpenMM::AngstromsPerNm
                                        * OpenMM::AngstromsPerNm
                                        *(bond_length-(Average_node_pair_length* OpenMM::NmPerAngstrom))
                                        *(bond_length-(Average_node_pair_length* OpenMM::NmPerAngstrom));

        }
    }
    return(bond_energy);
}




double Membrane::calculating_the_bond_energy_check(int p1, int p2, MyAtomInfo atoms[]){
    
    vector<double> bond;
    double Avg=get_avg_node_dist();
    double bond_length;
    double bond_energy=0;
    bond.resize(3);    
    for(int i=0; i<3; i++){    
    bond[i]= (atoms[p1].posInAng[i] -atoms[p2].posInAng[i]) *OpenMM::NmPerAngstrom;
   }
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);
   // cout<<"bond_length  "<<bond_length<<endl;
    //cout<<"bond_length-Avg  "<< bond_length-(Avg* OpenMM::NmPerAngstrom)<<endl;
    
    switch (spring_model){
//        case 1: //Fene
//        {
//            cout<<"FEne is under construction"<<endl;
//            
//        }
        case 2: //harmonic
        {
            bond_energy= 0.5*Spring_coefficient
                                      //  * OpenMM::KJPerKcal
                                        * OpenMM::AngstromsPerNm
                                        * OpenMM::AngstromsPerNm
                                        *(bond_length-(Avg* OpenMM::NmPerAngstrom))
                                        *(bond_length-(Avg* OpenMM::NmPerAngstrom));
            

        }
//    case 3: //realharmonic :)) (x4harmonic is the name but its potential is really 1/2*k*x^2)
//        {
//            bond_energy= 0.5*Spring_coefficient
//                                      //  * OpenMM::KJPerKcal
//                                        * OpenMM::AngstromsPerNm
//                                        * OpenMM::AngstromsPerNm
//                                        *(bond_length-(Avg* OpenMM::NmPerAngstrom))
//                                        *(bond_length-(Avg* OpenMM::NmPerAngstrom));
//
//        }
    }
    return(bond_energy);
    
}



double Membrane::calculating_the_bond_length_check(int p1, int p2, MyAtomInfo atoms[]){
    
    vector<double> bond;
//    double Avg=get_avg_node_dist();
    double bond_length;
    bond.resize(3);    
    for(int i=0; i<3; i++){    
    bond[i]= (atoms[p1].posInAng[i] -atoms[p2].posInAng[i]) *OpenMM::NmPerAngstrom;
   }
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);

    return(bond_length);
    
}
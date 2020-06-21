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
    bond[i]= (atoms[A+number_of_privious_mem_nodes].posInNm[i] -atoms[B+number_of_privious_mem_nodes].posInNm[i]);
    
    }
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);
    
    switch (spring_model){
        case 1: //Fene
        {
            cout<<"FEne is under construction"<<endl;
            
        }
        case 2: //harmonic
        {
            bond_energy= 0.5 * Spring_coefficient
                             * (bond_length - Average_node_pair_length)
                             * (bond_length - Average_node_pair_length);

        }
        case 5: //realharmonic :)) (x4harmonic is the name but its potential is really 1/2*k*x^2)
        {
            bond_energy= 0.5 * Spring_coefficient
                             * (bond_length - Average_node_pair_length)
                             * (bond_length - Average_node_pair_length);

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
    bond[i]= (atoms[p1].posInNm[i] -atoms[p2].posInNm[i]);
   }
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);
    
    switch (spring_model){

        case 2: //harmonic
        {
            bond_energy= 0.5 * Spring_coefficient
                             * (bond_length-Avg)
                             * (bond_length-Avg);
            

        }
    }
    return(bond_energy);
    
}



double Membrane::calculating_the_bond_length_check(int p1, int p2, MyAtomInfo atoms[]){
    
    vector<double> bond;
//    double Avg=get_avg_node_dist();
    double bond_length;
    bond.resize(3);    
    for(int i=0; i<3; i++){    
    bond[i]= (atoms[p1].posInNm[i] -atoms[p2].posInNm[i]);
   }
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);

    return(bond_length);
    
}

#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

void Membrane::monte_carlo_flip(){
    bool accept=0;
    //Randomly choosing a Dihedral and finding the neighbours.

    int triangle_A, triangle_B;
    int initial_pair;
    vector<int> A_neighbours_dihedral_index;
    vector<int> B_neighbours_dihedral_index;
    initial_pair=rand()%(Triangle_pair_list.size()-1);
    triangle_A=Triangle_pair_list[initial_pair][0];
    triangle_B=Triangle_pair_list[initial_pair][1];
    for (int i=0; i<Triangle_pair_list.size() ; i++){
        if ((Triangle_pair_list[i][0]==triangle_A and Triangle_pair_list[i][1]!=triangle_B)  || (Triangle_pair_list[i][1]==triangle_A and Triangle_pair_list[i][0]!=triangle_B)){
            A_neighbours_dihedral_index.push_back(i);}
        if ((Triangle_pair_list[i][0]==triangle_B and Triangle_pair_list[i][1]!=triangle_A)  || (Triangle_pair_list[i][1]==triangle_B and Triangle_pair_list[i][0]!=triangle_A)){
            B_neighbours_dihedral_index.push_back(i);}
        
    }
    if (A_neighbours_dihedral_index.size()==2 and B_neighbours_dihedral_index.size()==2){
        accept=1;
    }
    //if (accept==0){}
    
    //calculating the initial state energy
    //calculating the bond energy
    bool initial_or_final=0;
    double initial_bond_energy=calculating_the_bond_energy(initial_pair, initial_or_final);
    //calculating the bend energy
    
 
    
}

double Membrane::calculating_the_bond_energy(int index, bool initial_or_final){
    vector<double> bond;
    double bond_length;
    double bond_energy=0;
    int A,B;
    bond.resize(3);
    //In Triangle_pair_nodes nodes of the common bond are sorted in a way that locate between the uncommon nodes of triangles (indices 1 and 2) 
    if (initial_or_final==0){// initial, common bond 
        A=1;
        B=2;
        }
    else{
        A=0;
        B=3;
        }
    bond[0]= get_node_position(Triangle_Pair_Nodes[index][A], 0) - get_node_position(Triangle_Pair_Nodes[index][B], 0);
    bond[1]= get_node_position(Triangle_Pair_Nodes[index][A], 1) - get_node_position(Triangle_Pair_Nodes[index][B], 1);
    bond[2]= get_node_position(Triangle_Pair_Nodes[index][A], 2) - get_node_position(Triangle_Pair_Nodes[index][B], 2);    
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);
    
    switch (spring_model){
        case 1: //Fene
        {
            cout<<"FEne under construction"<<endl;
            
        }
        case 2: //harmonic
        {
            bond_energy= 0.5*Spring_coefficient
                                        * OpenMM::KJPerKcal
                                        * OpenMM::AngstromsPerNm
                                        * OpenMM::AngstromsPerNm
                                        *(bond_length-Average_node_pair_length)
                                        *(bond_length-Average_node_pair_length);
        }
    }
    return(bond_energy);
}

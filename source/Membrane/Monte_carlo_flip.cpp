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

double Membrane::calculating_the_bend_energy(int uncommon1, int common2, int common3, int uncommon4){
    double p1[3],p2[3],p3[3],p4[3];
    double p3p4[3],p2p4[3],p3p1[3],p2p1[3];
    double outward_vector[3];
    double N1[3], N2[3];
    double N1_length,N2_length;
    double cosine=1;
    double bending_energy;
    for (int index=0; index<3; index++){
        p1[index]=Node_Position[uncommon1][index];
        p2[index]=Node_Position[common2][index];
        p3[index]=Node_Position[common3][index];
        p4[index]=Node_Position[uncommon4][index];
        
        p3p4[index]=p4[index]-p3[index];
        p2p4[index]=p4[index]-p2[index];
        p3p1[index]=p1[index]-p3[index];
        p2p1[index]=p1[index]-p2[index];
        
        outward_vector[index]=p1[index]- COM_position[index];
    }
    crossvector(N1, p3p1, p2p1);
    crossvector(N2, p3p4, p2p4);
    if (innerproduct(N1,outward_vector)<0){
        N1[0]=-N1[0];
        N1[1]=-N1[1];
        N1[2]=-N1[2];
    }
    if (innerproduct(N2,outward_vector)<0){
        N2[0]=-N2[0];
        N2[1]=-N2[1];
        N2[2]=-N2[2];
    }
    N1_length=vector_length(N1);
    N2_length=vector_length(N2);
    if(N1_length !=0 and N2_length !=0){
    cosine= innerproduct(N1, N2)/(vector_length(N1)*vector_length(N2));}
    else{
        cout<<"warning! The trianglur mesh messed up. its impossible to calculate the bending energy." <<endl;
        EXIT_FAILURE;
    }
    bending_energy= Bending_coefficient * OpenMM::KJPerKcal
                                        *(1.00001-cosine);
    return(bending_energy);
}
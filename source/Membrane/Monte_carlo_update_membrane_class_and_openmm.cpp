#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

void Membrane::update_Membrane_class_and_openmm(int initial_pair,int triangle_A,int triangle_B, int new_neighbour_dihedrals[4][6], MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals){
   // check_monte_carlo=1;
    
    int u1=Triangle_Pair_Nodes[initial_pair][0];
    int c2=Triangle_Pair_Nodes[initial_pair][1];
    int c3= Triangle_Pair_Nodes[initial_pair][2];
    int u4=Triangle_Pair_Nodes[initial_pair][3];
    //cout<<"u1 "<<u1<<"  c2  "<<c2<<"  c3  "<<c3<<"  u4  "<<u4<<endl;
    
        
    //updating Node_bond_list
    int bond_index;
    Update_Membrane(initial_pair,triangle_A,triangle_B, new_neighbour_dihedrals, bond_index);

    
    //updating the openmm
    ////updating the struct bonds and setbondparameter of spring
    vector<double> parameters;
    bonds[bond_index].atoms[0]=u1;
    bonds[bond_index].atoms[1]=u4;
    
    for (int i=0; i< omm->EV.size() ; i++)
    {
        omm->EV[i]->setExclusionParticles(bond_index, u1, u4);
    }
    
    for (int i=0; i< omm->LJ.size() ; i++)
    {
        omm->LJ[i]->setExclusionParticles(bond_index, u1, u4);
    }
    
    switch (spring_model){
        case 1: //Fene
            { 
                    cout<<"FEne is under construction"<<endl; //openmm setparameter
            }
        case 2: //harmonic
            {
                parameters.resize(2);
                parameters[0]=bonds[bond_index].nominalLengthInNm; //I must check the bond index in case of having other objects
                parameters[1]=bonds[bond_index].stiffnessInKJPerNm2;
                omm->harmonic->setBondParameters(bond_index,u1,u4, parameters[0], parameters[1] );
                    
            }
    }
    
    ////updating dihedral structs and setBondParameters of Dihedral Forces.
    vector<double> bendingparameter;
    bendingparameter.resize(1);
    bendingparameter[0]=Bending_coefficient* OpenMM::KJPerKcal;
    
    dihedrals[initial_pair].atoms[0]=c2;
    dihedrals[initial_pair].atoms[1]=u1;
    dihedrals[initial_pair].atoms[2]=u4;
    dihedrals[initial_pair].atoms[3]=c3;
    omm->Dihedral[0]->setBondParameters(initial_pair,Triangle_Pair_Nodes[initial_pair], bendingparameter);
     
    for(int i=0; i<4; i++){
        int index= new_neighbour_dihedrals[i][4];
        omm->Dihedral[0]->setBondParameters(index,Triangle_Pair_Nodes[index], bendingparameter);
    }
  

 /* for(int k=0; k<Num_of_Triangles; k++){
        cout<<"newww Triangle list  "<<k<<"   "<<  Triangle_list[k][0]<<"  "<<Triangle_list[k][1]<<"   "<<Triangle_list[k][2]<<endl;}

    for(int p=0; p<Num_of_Triangle_Pairs; p++){
      cout<<"new neighbours list   "<<Triangle_pair_list[p][0]<<"  "<<Triangle_pair_list[p][1]<<endl; 
      cout<<"new neighbour nodes   "<<Triangle_Pair_Nodes[p][0]<<"  "<<Triangle_Pair_Nodes[p][1]<<"  "<<Triangle_Pair_Nodes[p][2]<<"  "<<Triangle_Pair_Nodes[p][3]<<endl;
      
    }*/
     
}



void Membrane::Update_Membrane(int initial_pair,int triangle_A,int triangle_B, int new_neighbour_dihedrals[4][6], int& bond_index){
    int u1=Triangle_Pair_Nodes[initial_pair][0];
    int c2=Triangle_Pair_Nodes[initial_pair][1];
    int c3= Triangle_Pair_Nodes[initial_pair][2];
    int u4=Triangle_Pair_Nodes[initial_pair][3];
    //cout<<"u1 "<<u1<<"  c2  "<<c2<<"  c3  "<<c3<<"  u4  "<<u4<<endl;
    
        
    //updating Node_bond_list
    
    for(int i=0; i<Num_of_Node_Pairs; i++){
        if((Node_Bond_list[i][0]==c2 and Node_Bond_list[i][1]==c3) ||(Node_Bond_list[i][1]==c2 and Node_Bond_list[i][0]==c3)){
            bond_index=i;
            Node_Bond_list[i][0]=u1;
            Node_Bond_list[i][1]=u4;
        }
    }
    
    
    //updating triangle list
    Triangle_list[triangle_A][0]=u1;
    Triangle_list[triangle_A][1]=c2;
    Triangle_list[triangle_A][2]=u4;
    
    
    Triangle_list[triangle_B][0]=u1;
    Triangle_list[triangle_B][1]=c3;
    Triangle_list[triangle_B][2]=u4;
    
    
    //updating Triangle_Pair_Nodes
    Triangle_Pair_Nodes[initial_pair][0]=c2;
    Triangle_Pair_Nodes[initial_pair][1]=u1;
    Triangle_Pair_Nodes[initial_pair][2]=u4;
    Triangle_Pair_Nodes[initial_pair][3]=c3;
     
    for(int i=0; i<4; i++){
        int index= new_neighbour_dihedrals[i][4];
        for (int j=0; j<4; j++){
            Triangle_Pair_Nodes[index][j]=new_neighbour_dihedrals[i][j];}
            
    //updating Triangle_Pair_list
        if(new_neighbour_dihedrals[i][5]==1){
            if(Triangle_pair_list[index][0]==triangle_A){
                Triangle_pair_list[index][0]=triangle_B;
            }else if(Triangle_pair_list[index][1]==triangle_A){
                Triangle_pair_list[index][1]=triangle_B;
            }else if(Triangle_pair_list[index][0]==triangle_B){
                Triangle_pair_list[index][0]=triangle_A;
            }else if(Triangle_pair_list[index][1]==triangle_B){
                Triangle_pair_list[index][1]=triangle_A;
            }
        }
    }
}

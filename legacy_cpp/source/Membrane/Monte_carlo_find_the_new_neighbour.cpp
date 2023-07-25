#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

void Membrane::find_the_new_neighbour(int neighbour_id[6], int previous_dihedral_index , int initial_pair, bool A_or_B){
    //cout<<"find_the_new_neighbours"<<endl;
//    int Uncommon_node;
    int uncommon1=Triangle_Pair_Nodes[initial_pair][0];
    int common2=Triangle_Pair_Nodes[initial_pair][1];
    int common3=Triangle_Pair_Nodes[initial_pair][2];
    int uncommon4=Triangle_Pair_Nodes[initial_pair][3];
    neighbour_id[4]=previous_dihedral_index;
    if(A_or_B==0){
        if(Triangle_Pair_Nodes[previous_dihedral_index][0]==common2){
          neighbour_id[0]=uncommon4;
          for (int i=1; i<4; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=1;}
        
        if(Triangle_Pair_Nodes[previous_dihedral_index][3]==common2){
          neighbour_id[3]=uncommon4;
          for (int i=0; i<3; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=1;}
        
        if(Triangle_Pair_Nodes[previous_dihedral_index][0]==common3){
          neighbour_id[0]=uncommon4;
          for (int i=1; i<4; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=0;}
          
        if(Triangle_Pair_Nodes[previous_dihedral_index][3]==common3){
          neighbour_id[3]=uncommon4;
          for (int i=0; i<3; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=0;}
        
    }
    else{
        if(Triangle_Pair_Nodes[previous_dihedral_index][0]==common2){
          neighbour_id[0]=uncommon1;
          for (int i=1; i<4; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=0;}
        
        if(Triangle_Pair_Nodes[previous_dihedral_index][3]==common2){
          neighbour_id[3]=uncommon1;
          for (int i=0; i<3; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=0;}
        
        if(Triangle_Pair_Nodes[previous_dihedral_index][0]==common3){
          neighbour_id[0]=uncommon1;
          for (int i=1; i<4; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=1;}
          
        if(Triangle_Pair_Nodes[previous_dihedral_index][3]==common3){
          neighbour_id[3]=uncommon1;
          for (int i=0; i<3; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=1;}
        
    }
    
   
        
}
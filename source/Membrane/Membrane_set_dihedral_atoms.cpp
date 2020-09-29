//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 24/09/2020.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"


void Membrane::set_dihedral_atoms(void){
    dihedral_atoms.clear();
    dihedral_atoms.resize(Num_of_Triangle_Pairs);
    for (auto &nodes: dihedral_atoms){
        nodes.resize(4,-1);
    }
    
    
    update_COM_position();
//    double COM[3]={0,0,0};
    int pos1,pos2,pos3,pos4;
    double  b1[3], b2[3], b3[3];
    double N1[3], N2[3], Nout[3];
    
    for (int i=0; i<Num_of_Triangle_Pairs; i++)
    {
        pos1=Triangle_Pair_Nodes[i][0];
        pos2=Triangle_Pair_Nodes[i][1];
        pos3=Triangle_Pair_Nodes[i][2];
        pos4=Triangle_Pair_Nodes[i][3];
        
        for (int index=0; index<3; index++) {
            b1[index]  = Node_Position[pos2][index]-Node_Position[pos1][index];
            b2[index]  = Node_Position[pos3][index]-Node_Position[pos2][index];
            b3[index]  = Node_Position[pos4][index]-Node_Position[pos3][index];
            Nout[index]= Node_Position[pos2][index]-COM_position[index];
//            Nout[index]= Node_Position[pos2][index]-COM[index];
        }
        
        crossvector(N1, b1, b2);
        crossvector(N2, b2, b3);
        
        
        if (innerproduct(N1, Nout)<0) {
            
            dihedral_atoms[i][0]=pos1;
            dihedral_atoms[i][1]=pos2;
            dihedral_atoms[i][2]=pos3;
            dihedral_atoms[i][3]=pos4;
        } else {
            
            dihedral_atoms[i][0]=pos1;
            dihedral_atoms[i][1]=pos3;
            dihedral_atoms[i][2]=pos2;
            dihedral_atoms[i][3]=pos4;
        }
    }
    
    
}

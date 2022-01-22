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
    
    for (int i=0; i<Num_of_Triangle_Pairs; i++)
    {
        vector<int> dihedral_atom_indices = check_dihedral_direction(Triangle_Pair_Nodes[i]);
        dihedral_atoms[i]=dihedral_atom_indices;
    }
}

vector<int> Membrane::check_dihedral_direction(vector<int> dihedral_atom_indices){
    int pos1,pos2,pos3,pos4;
    double  b1[3], b2[3], b3[3];
    double N1[3], N2[3], Nout[3];
    
    vector<int> correct_arrengment(4,-1);
    
    pos1=dihedral_atom_indices[0];
    pos2=dihedral_atom_indices[1];
    pos3=dihedral_atom_indices[2];
    pos4=dihedral_atom_indices[3];
    
    for (int index=0; index<3; index++) {
        b1[index]  = Node_Position[pos2][index]-Node_Position[pos1][index];
        b2[index]  = Node_Position[pos3][index]-Node_Position[pos2][index];
        b3[index]  = Node_Position[pos4][index]-Node_Position[pos3][index];
        if (UseXYZinMembrane) {
            Nout[index]= Node_Position[pos2][index]-XYZinMembrane[index];
        } else {
            Nout[index]= Node_Position[pos2][index]-COM_position[index];
        }
    }
    
    crossvector(N1, b1, b2);
    crossvector(N2, b2, b3);
    
    if (innerproduct(N1, Nout)<0) {
        correct_arrengment[0]=pos1;
        correct_arrengment[1]=pos2;
        correct_arrengment[2]=pos3;
        correct_arrengment[3]=pos4;
    } else {
        correct_arrengment[0]=pos1;
        correct_arrengment[1]=pos3;
        correct_arrengment[2]=pos2;
        correct_arrengment[3]=pos4;
    }
    return correct_arrengment;;
}

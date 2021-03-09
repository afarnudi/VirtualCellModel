//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//
#include <iostream>
#include "Global_functions.hpp"

void print_statistics(int num_of_atoms,
                      int num_of_bonds,
                      int num_of_dihedrals,
                      int num_of_angles,
                      vector<Membrane> &Membranes,
                      vector<Chromatin> &Chromatins){
    
    for (int i=0; i<Membranes.size(); i++) {
        if (Membranes[i].get_spring_model() == potentialModelIndex.Model["None"]) {
            num_of_bonds-=Membranes[i].get_num_of_node_pairs();
        }
        if (Membranes[i].get_bending_model() == potentialModelIndex.Model["None"]) {
            num_of_dihedrals-=Membranes[i].get_num_of_triangle_pairs();
        }
    }
    for (int i=0; i<generalParameters.Num_of_Chromatins; i++) {
        if (Chromatins[i].get_spring_model() == potentialModelIndex.Model["None"]) {
            num_of_bonds-=Chromatins[i].get_num_of_bonds();
        }
    }
    cout<<num_of_atoms<<" Atoms"<<endl;
    cout<<num_of_bonds<<" Bonds"<<endl;
    cout<<num_of_dihedrals<<" Dihedrals"<<endl;
    cout<<num_of_angles<<" Angles"<<endl;
}

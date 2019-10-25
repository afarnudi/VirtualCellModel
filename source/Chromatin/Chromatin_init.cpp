//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Chromatin.h"
#include <algorithm>
//void build_random_chain(void);

using std::cout;
using std::endl;

void Chromatin::initialise(void){
    vector<vector<double> > zero_vec(Num_of_Nodes, vector<double> (Num_of_Nodes));
    Contact_Matrix=zero_vec;
    
    ABC_index.resize(Num_of_Nodes);
    generate(ABC_index.begin(), ABC_index.end(), [&r=num_of_node_types]() {
        return rand() % r;
    });
    cout<<"\nInitialising the Chromatin Class..."<<endl;
    build_random_chain();
    shift_node_positions();
    
    while (label.length()>3) {
        label.pop_back();
    }
    while (label.length()<3) {
        label += "0";
    }
    
    if (index>=100) {
        cout<<"Trouble with PDB labeling:\nPDB name: class_label + calss_index + class_node_type. \nToo many classes (100>) to make a PDB name with the correct format.\n";
        std::exit(EXIT_FAILURE);
    } else if (index>=10){
        if(num_of_node_types>=10){
            cout<<"Trouble with PDB labeling:\nPDB name: class_label + calss_index + class_node_type. \nToo many node types (# of classes >10 and # of node types >10) to make a PDB name with the correct format.\n";
            std::exit(EXIT_FAILURE);
        } else {
            label.pop_back();
            label.pop_back();
            label += std::to_string(index);
        }
    } else {
        if(num_of_node_types>=10){
            label.pop_back();
            label.pop_back();
            label += std::to_string(index);
        } else {
            label.pop_back();
            label += std::to_string(index);
        }
    }
    cout<<"Chromatin class initiated.\n\n";
    
}

void Chromatin::initialise(double min_radius){
    vector<vector<double> > zero_vec(Num_of_Nodes, vector<double> (Num_of_Nodes));
    Contact_Matrix=zero_vec;
    
    ABC_index.resize(Num_of_Nodes);
    generate(ABC_index.begin(), ABC_index.end(), [&r=num_of_node_types]() {
        return rand() % r;
    });
    
    cout<<"\nInitialising the Chromatin Class..."<<endl;
    build_random_chain();
    Pack(min_radius-2);
    shift_node_positions();
    cout<<"Chromatin class initiated.\n";
    
}

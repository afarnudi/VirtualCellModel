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
    
    pdb_label_check();
    
    Contact_Matrix.resize(Num_of_Nodes);
    for (int i=0; i<Num_of_Nodes; i++) {
        Contact_Matrix[i].resize(Num_of_Nodes, 0);
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

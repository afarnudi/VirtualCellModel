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
//    vector<vector<double> > zero_vec(Num_of_Nodes, vector<double> (Num_of_Nodes));
//    Contact_Matrix=zero_vec;
    
    ABC_index.resize(Num_of_Nodes);
    generate(ABC_index.begin(), ABC_index.end(), [&r=num_of_node_types]() {
        return rand() % r;
    });
    cout<<TPINK<<"\nInitialising the Chromatin Class..."<<TRESET<<endl;
    build_random_chain();
    cout<<"Generated a self avoiding random walk chain."<<endl;
    shift_node_positions();
    cout<<"Shifted positions. to "<<Shift_position_xyzVector[0]<< " "<<Shift_position_xyzVector[1]<< " "<<Shift_position_xyzVector[2]<< " " <<endl;
    shift_node_velocities();
    cout<<"Shifted velocities."<<endl;
    
    pdb_label_check();
    
    if (ExportGeneratedCoordinates) {
        export_coordinates();
    }
//    Contact_Matrix.resize(Num_of_Nodes);
//    for (int i=0; i<Num_of_Nodes; i++) {
//        Contact_Matrix[i].resize(Num_of_Nodes, 0);
//    }
    
    cout<<TSUCCESS<<"Chromatin class initiated.\n"<<TRESET<<
                  "**************************\n"<<endl;
    
}

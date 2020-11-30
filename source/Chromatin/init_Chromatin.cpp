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
    cout<<TPINK<<"\nInitialising the Chromatin Class..."<<TRESET<<endl;
    
    if (ImportCoordinates) {
        import_coordinates(import_file_name);
    } else {
        build_random_chain();
    }
    
    if (rescale_factor!=1 && rescale_factor!=0) {
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Position[i][0]*=rescale_factor;
            Node_Position[i][1]*=rescale_factor;
            Node_Position[i][2]*=rescale_factor;
        }
        
    }
    
    ABC_index.resize(Num_of_Nodes);
    generate(ABC_index.begin(), ABC_index.end(), [&r=num_of_node_types]() {
        return rand() % r;
    });
    
    set_bond_nominal_lengths();
    set_node_radius();
    
    if (spring_model == GenConst::potential.Model["None"]) {
        cout<<TWARN<<"\nChromatin spring model is set to 'None'."<<TRESET<<endl;
    }
    
    shift_node_positions();
    cout<<"Shifted positions. to "<<Shift_position_xyzVector[0]<< " "<<Shift_position_xyzVector[1]<< " "<<Shift_position_xyzVector[2]<< " " <<endl;
    shift_node_velocities();
//    cout<<"Shifted velocities."<<endl;
    
    pdb_label_check();
    
    
//    Contact_Matrix.resize(Num_of_Nodes);
//    for (int i=0; i<Num_of_Nodes; i++) {
//        Contact_Matrix[i].resize(Num_of_Nodes, 0);
//    }
    cout<<TSUCCESS<<"Chromatin class initiated.\n"<<TRESET<<
                  "**************************\n"<<endl;
    
}

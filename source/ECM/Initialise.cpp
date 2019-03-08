//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.h"


void ECM::initialise(string Mesh_file_name, int dimension){
//    T_Kinetic_Energy.resize(100);
    cout<<"Initialising the ECM Class..."<<endl;
    if (dimension == 2) {
        read_gmesh_file_2D(Mesh_file_name);
    }
    
    output_file_neme=Mesh_file_name;
    cout<<"# of Nodes = "<<Num_of_Nodes<<endl;
    cout<<"# of Triangles = "<<Num_of_Triangles<<endl;
    Node_Bond_identifier();
    Node_neighbour_list_constructor();
    
//    Normal_direction_Identifier();
//    Triangle_pair_counter();
//    cout<<"# of triangle pairs="<<Num_of_Triangle_Pairs<<endl;
//    if (Num_of_Triangle_Pairs != 3*(Triangle_list.size())/2){
//        cout<<"Warning! some triangles have less or more neighbour than 3"<<endl;
//        
//    }
//    //        Triangle_Pair_and_Node_Bonds_Identifier();
    
//    Triangle_pair_identifier();

//    
    cout<<"\nECM class initiated.\n******************************\n\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}

//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"


void Actin::initialise(string Mesh_file_name){
//    T_Kinetic_Energy.resize(100);
    cout<<"Initialising the Actin Class..."<<endl;
    read_gmesh_file(Mesh_file_name);
//    output_file_neme=Mesh_file_name;
//    Radius= sqrt((Node_Position[0][0]-X_in)*(Node_Position[0][0]-X_in) + (Node_Position[0][1]-Y_in)*(Node_Position[0][1]-Y_in) + (Node_Position[0][2]-Z_in)*(Node_Position[0][2]-Z_in));
//    cout<<"\nRadius="<<Radius<<endl;
//    cout<<"# of Nodes="<<Num_of_Nodes<<endl;
//    cout<<"# of triangles="<<Num_of_Triangles<<endl;
//    Normal_direction_Identifier();
//    Triangle_pair_counter();
//    cout<<"# of triangle pairs="<<Num_of_Triangle_Pairs<<endl;
//    if (Num_of_Triangle_Pairs != 3*(Triangle_list.size())/2){
//        cout<<"Warning! some triangles have less or more neighbour than 3"<<endl;
//        
//    }
//    //        Triangle_Pair_and_Node_Bonds_Identifier();
//    Node_Bonds_identifier();
//    Node_neighbour_list_constructor();
//    Triangle_pair_identifier();
//    check();
//    
    cout<<"\nMembrane class initiated.\n******************************\n\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}

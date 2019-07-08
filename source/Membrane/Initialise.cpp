//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::initialise(std::string Mesh_file_name){
//    T_Kinetic_Energy.resize(100);
    cout<<"Initialising the Membrane Class..."<<endl;
    if (mesh_format==1){
        read_gmesh_file(Mesh_file_name);
    }else if (mesh_format==2){
        read_ply_file(Mesh_file_name);
    }
    output_file_neme=Mesh_file_name;
//    if (rescale==true){
//        Rescale(rescale_factor);
//        cout<<"mesh rescaled by factor"<< rescale_factor <<endl;
//    }
    
//    if(Shift==true){
//    shift_position ( Shift_in_X_direction, Shift_in_Y_direction, Shift_in_Z_direction);
////    cout<<"The Whole Membrane is Shifted"<<endl;
////    cout<<"Shift in X direction = "<< Shift_in_X_direction<<endl;
////    cout<<"Shift in Y direction = "<< Shift_in_Y_direction<<endl;
////    cout<<"Shift in Z direction = "<< Shift_in_Z_direction<<endl;
//    }
    
    Radius= sqrt((Node_Position[0][0]-X_in)*(Node_Position[0][0]-X_in) + (Node_Position[0][1]-Y_in)*(Node_Position[0][1]-Y_in) + (Node_Position[0][2]-Z_in)*(Node_Position[0][2]-Z_in));
    cout<<"\nRadius="<<Radius<<endl;
    cout<<"# of Nodes="<<Num_of_Nodes<<endl;
    cout<<"# of triangles="<<Num_of_Triangles<<endl;
    Normal_direction_Identifier();
    Triangle_pair_counter();
    cout<<"# of triangle pairs="<<Num_of_Triangle_Pairs<<endl;
    if (Num_of_Triangle_Pairs != 3*(Triangle_list.size())/2){
        cout<<"Warning! some triangles have less or more neighbour than 3"<<endl;
    }
    //        Triangle_Pair_and_Node_Bonds_Identifier();
    Node_Bonds_identifier();
    Node_neighbour_list_constructor();
    Triangle_pair_identifier();
    //DamperCheck.resize(Num_of_Node_Pairs);
    //SinusCheck.resize(Num_of_Node_Pairs);
    //check();
    ECM_Node_neighbour_list.resize(Num_of_Nodes);
    shift_position(Shift_in_X_direction, Shift_in_Y_direction, Shift_in_Z_direction);
    cout<<"\nMembrane class initiated.\n******************************\n\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}

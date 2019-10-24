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
    
    Radius= sqrt((Node_Position[0][0]-X_in)*(Node_Position[0][0]-X_in) + (Node_Position[0][1]-Y_in)*(Node_Position[0][1]-Y_in) + (Node_Position[0][2]-Z_in)*(Node_Position[0][2]-Z_in));
    cout<<"\nRadius="<<Radius<<endl;
    cout<<"# of Nodes="<<Num_of_Nodes<<endl;
    cout<<"# of triangles="<<Num_of_Triangles<<endl;
   

    //this part is for testing  just one bond (be carful to pass a mesh with just 2 nodes to the programm in case of using it)
/*
    Node_Bond_list.resize(1);
    Node_Bond_list[0].resize(2);
    Node_Bond_list[0][0]=0;
    Node_Bond_list[0][1]=1;
    
    Num_of_Node_Pairs=1;
    Num_of_Triangle_Pairs=0;
    cout<<"\nMembrane class initiated.\n******************************\n\n";
*/
    Normal_direction_Identifier();
    Triangle_pair_counter();
    cout<<"# of triangle pairs="<<Num_of_Triangle_Pairs<<endl;
    if (Num_of_Triangle_Pairs != 3*(Triangle_list.size())/2){
        cout<<"Warning! some triangles have less or more neighbour than 3"<<endl;
    }
    Node_Bonds_identifier();
    Node_neighbour_list_constructor();
    Triangle_pair_identifier();
  

    ECM_Node_neighbour_list.resize(Num_of_Nodes);
    //shift_position(Shift_in_X_direction, Shift_in_Y_direction, Shift_in_Z_direction);
    check();

    check_radius_update_values();
    shift_velocity(x_speed, y_speed, z_speed);

    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;


}

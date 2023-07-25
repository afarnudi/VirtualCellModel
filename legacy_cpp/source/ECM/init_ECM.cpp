//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.h"


void ECM::initialise(string Mesh_file_name){
//    T_Kinetic_Energy.resize(100);
    cout<<TECM<<"Initialising the ECM Class..."<<TRESET<<endl;
    if (MeshType == "2d") {
        read_gmesh_file_2D(Mesh_file_name);
    } else if (MeshType == "3d"){
        read_gmesh_file_3D(Mesh_file_name);
    } else if (MeshType == "3dsquare"){
        read_gmesh_file_3D_square(Mesh_file_name);
    }
    
    output_file_neme=Mesh_file_name;
    cout<<"# of Nodes = "<<Num_of_Nodes<<endl;
    cout<<"# of Triangles = "<<Num_of_Triangles<<endl;
    
    if(MeshType=="2d"){
        Node_Bond_identifier_2D();
    }
    else if(MeshType=="3d"){
        Node_Bond_identifier_3D();
    }
    else if(MeshType=="3dsquare"){
        Node_Bond_identifier_3D_square();
    }
    cout<<"# of bonds = "<<Num_of_Node_Pairs<<endl;
    Node_neighbour_list_constructor();
    
    if (MeshType == "3d" ) {
       // normal_direction_Identifier();
    }
    
    //shift_node_positions();
   // shift_velocity(x_speed, y_speed, z_speed);
    shift_velocity(Shift_velocities_xyzVector[0], Shift_velocities_xyzVector[1], Shift_velocities_xyzVector[2]);
    shift_position(Shift_position_xyzVector[0], Shift_position_xyzVector[1], Shift_position_xyzVector[2]);
    check();
    
    while (label.length()>3) {
        label.pop_back();
    }
    while (label.length()<3) {
        label += "0";
    }
    if (index>=10){
        label.pop_back();
        label += std::to_string(index);
    } else {
        label += std::to_string(index);
    }
    
    cout<<TSUCCESS<<"\nECM class initiated.\n"<<TRESET<<
                    "********************\n"<<endl;
}

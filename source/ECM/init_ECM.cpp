//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.h"


void ECM::initialise(int dimension){
//    T_Kinetic_Energy.resize(100);
    cout<<TECM<<"Initialising the ECM Class..."<<TRESET<<endl;
    if (dimension == 2) {
        read_gmesh_file_2D(Mesh_file_name);
    } else if (dimension == 3){
        read_gmesh_file_3D(Mesh_file_name);
    } else if (dimension == 4){
        read_gmesh_file_3D_square(Mesh_file_name);
    }
    
    output_file_neme=Mesh_file_name;
    cout<<"# of Nodes = "<<Num_of_Nodes<<endl;
    cout<<"# of Triangles = "<<Num_of_Triangles<<endl;
    
    if(dimension==2){
        Node_Bond_identifier_2D();
    }
    else if(dimension==3){
        Node_Bond_identifier_3D();
    }
    else if(dimension==4){
        Node_Bond_identifier_3D_square();
    }
    cout<<"# of bonds = "<<Num_of_Node_Pairs<<endl;
    Node_neighbour_list_constructor();
    
    if (dimension == 3 ) {
       // normal_direction_Identifier();
    }
    
    shift_node_positions();
    shift_velocity(x_speed, y_speed, z_speed);
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

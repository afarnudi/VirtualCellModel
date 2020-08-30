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
    if (!GenConst::Testmode) {
        cout<<TBLUE<<"\nInitialising the Membrane Class..."<<TRESET<<endl;
    }
    
   

    if (mesh_format==1){
        read_gmesh_file(Mesh_file_name);
    }else if (mesh_format==2){
        read_ply_file(Mesh_file_name);
    }
    output_file_neme=Mesh_file_name;
    
    Radius= sqrt((Node_Position[0][0]-X_in)*(Node_Position[0][0]-X_in) + (Node_Position[0][1]-Y_in)*(Node_Position[0][1]-Y_in) + (Node_Position[0][2]-Z_in)*(Node_Position[0][2]-Z_in));
    if (!GenConst::Testmode) {
        
        cout<<"Number of ... :\n";
        cout<<"Nodes\t"<<Num_of_Nodes<<endl;
        cout<<"Triangles\t"<<Num_of_Triangles<<endl;
    }
    
    Normal_direction_Identifier();
    Triangle_pair_counter();
    if (!GenConst::Testmode) {
        cout<<"Triangle pairs\t"<<Num_of_Triangle_Pairs<<endl;
    }
    
    if (Num_of_Triangle_Pairs != 3*(Triangle_list.size())/2){
        cout<<TWARN<<"Warning"<<TRESET<<"! some triangles have less or more neighbour than 3"<<endl;
    }
    Node_Bonds_identifier();
    if (!GenConst::Testmode) {
        cout<<"Node pairs (Bonds) "<<Num_of_Node_Pairs<<endl;
    }
    
    cout<<"\nOther properties:\n";
    cout<<"\nRadius\t"<<Radius<<endl;
    Triangle_pair_identifier();
    
    Node_neighbour_list_constructor();
    Bond_triangle_neighbour_list_constructor();
    
    cout<<endl;
    check();

    check_radius_update_values();
    shift_velocity(x_speed, y_speed, z_speed);
    
    
    //limiting the labels to 4 charachters for use in the pdb writer
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
    
    if (spring_model == 1) {
        if (FENE_k == 0 || FENE_epsilon == 0 || FENE_max == 0 ) {
            cout<<TBLINK<<TWWARN<<"Warning"<<TRESET<<". Membrane spring model set to FENE but FENE parameters not set in the membrane configuration file. Please make sure you have set the following parameters: \nFENE_eps\nFENE_k\nFENE_min\nFENE_max (cannot be zero)\n";
            exit(EXIT_FAILURE);
        }
    }
    
    if (initial_random_rotation_coordinates){
        if (!GenConst::Testmode) {
            cout<<TWARN<<"Randomly rotating the mesh\n"<<TRESET;
        }
        
        srand (time(NULL));
        double scratch = rand();
        double phi   = ((double) rand() / (RAND_MAX))*2*M_PI;
        double theta = ((double) rand() / (RAND_MAX))*M_PI;
//        cout<<"theta: "<<theta<<"\n";
//        cout<<"phi  : "<<phi<<"\n";
        rotate_coordinates(theta, phi);
    }
    
//    if (GenConst::Wantvoronoi){
        node_voronoi_area.resize(Num_of_Nodes,0);
//    }
    shift_position(Shift_in_X_direction, Shift_in_Y_direction, Shift_in_Z_direction);
    if (!GenConst::Testmode) {
        cout<<"Bending energy = "<<calculate_bending_energy()<<endl;
        
        cout<<TSUCCESS<<"\nMembrane class initiated."<<TRESET<<"\n******************************\n\n";
    }
    

    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;


}

void Membrane::analysis_init(std::string Mesh_path){
    string mesh_extension = Mesh_path;
    mesh_extension.erase(mesh_extension.begin(),mesh_extension.end()-3 );
    if (mesh_extension=="ply"){
        read_ply_file(Mesh_path);
    } else {
        read_gmesh_file(Mesh_path);
    }
    
    Normal_direction_Identifier();
    Triangle_pair_counter();
    
    Node_Bonds_identifier();
    Triangle_pair_identifier();
    
    Node_neighbour_list_constructor();
    Bond_triangle_neighbour_list_constructor();
    
    node_voronoi_area.resize(Num_of_Nodes,0);
    
    
}

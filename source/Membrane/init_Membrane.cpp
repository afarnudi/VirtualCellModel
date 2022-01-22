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
    if (!generalParameters.Testmode) {
        cout<<TMEM<<"\nInitialising the Membrane Class..."<<TRESET<<endl;
    }
    
    read_mesh_file(Mesh_file_name);
    
    output_file_neme=Mesh_file_name;
    if (!generalParameters.Testmode) {
        
        cout<<"Number of ... :\n";
        cout<<"Nodes\t"<<Num_of_Nodes<<endl;
        cout<<"Triangles\t"<<Num_of_Triangles<<endl;
    }
    update_COM_position();
    Normal_direction_Identifier();
//    Triangle_pair_counter();
    Triangle_pair_identifier();
    if (!generalParameters.Testmode) {
        cout<<"Triangle pairs\t"<<Num_of_Triangle_Pairs<<endl;
    }
    
    if (Num_of_Triangle_Pairs != 3*(Triangle_list.size())/2){
        cout<<TWARN<<"Warning"<<TRESET<<"! some triangles have less or more neighbour than 3"<<endl;
    }
    Node_Bonds_identifier();
    if (!generalParameters.Testmode) {
        cout<<"Node pairs (Bonds) "<<Num_of_Node_Pairs<<endl;
    }
    
    cout<<"\nOther properties:\n";
    update_average_Membrane_radius();
    cout<<"\nRadius\t"<<Radius<<endl;
//    Triangle_pair_identifier();
    
    Node_neighbour_list_constructor();
    Bond_triangle_neighbour_list_constructor();
    
    assign_surface_volume_constraints();
    
    
    if (AddRandomModes) {
        randomundulationgenerator();
    }
    check();
    set_bond_nominal_length();
    set_dihedral_atoms();
    set_bending_nominal_angle();
    set_node_radius();
    calculate_volume_and_surface_area();
    
    if (bending_model==potentialModelIndex.Model["Julicher1996"]) {
        mean_curvature_init();
    }
    
    if (New_Radius!=-1) {
        check_radius_update_values();
    }
    
    shift_velocity(Shift_velocities_xyzVector[0], Shift_velocities_xyzVector[1], Shift_velocities_xyzVector[2]);
    
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
    
//    if (spring_model == potentialModelIndex.Model["FENE"]) {
//        if (FENE_k == 0 || FENE_epsilon == 0 || FENE_max == 0 ) {
//            string errorMessage = TWARN;
//            errorMessage+="Warning";
//            errorMessage+= TRESET;
//            errorMessage+= ". Membrane spring model set to FENE but FENE parameters not set in the membrane configuration file. Please make sure you have set the following parameters: \nFENE_eps\nFENE_k\nFENE_min\nFENE_max (cannot be zero)\n";
//            throw std::runtime_error(errorMessage);
//        }
//    } else
        if (spring_model == potentialModelIndex.Model["None"]) {
        cout<<TWARN<<"\nMembrnae spring model is set to 'None'."<<TRESET<<endl;
    }
    
    if (initial_random_rotation_coordinates){
        if (!generalParameters.Testmode) {
            cout<<TWARN<<"\nRandomly rotating the mesh\n"<<TRESET;
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
//        node_voronoi_area.resize(Num_of_Nodes,0);
//    }
    shift_position(Shift_position_xyzVector[0], Shift_position_xyzVector[1], Shift_position_xyzVector[2]);
    if (!generalParameters.Testmode) {
        cout<<TSUCCESS<<"\nMembrane class initiated."<<TRESET<<
                        "\n*************************\n"<<endl;
    }
}

void Membrane::analysis_init(std::string Mesh_path){
    string mesh_extension = Mesh_path;
    mesh_extension.erase(mesh_extension.begin(),mesh_extension.end()-3 );
    
    if (mesh_extension=="ply"){
        read_ply_file(Mesh_path);
    } else {
        read_gmesh_file(Mesh_path);
    }
    update_average_Membrane_radius();
    
    Normal_direction_Identifier();
    Triangle_pair_counter();
    
    Node_Bonds_identifier();
    Triangle_pair_identifier();
    
    Node_neighbour_list_constructor();
    Bond_triangle_neighbour_list_constructor();
}




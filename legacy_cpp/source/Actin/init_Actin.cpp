//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"

using std::cout;
using std::endl;

void Actin::initialise(string Mesh_file_name){
    //    T_Kinetic_Energy.resize(100);
    cout<<TACT<<"Initialising the Actin Class..."<<TRESET<<endl;
    if(mesh_format=="msh"){
        if (MeshType=="normal") {
            read_gmesh_file(Mesh_file_name);
        } else if (MeshType=="sajjad"){
            read_gmesh_file_2(Mesh_file_name);
        }
        
    } else if (mesh_format=="actin"){
       read_actin_file();
    }
    
    cout<<"# of Nodes="<<Num_of_Nodes<<endl;
    if(MeshType=="normal"){
        if (mesh_format=="actin") {
            cout<<"# of bonds = "<<Num_of_Node_Pairs<<endl;
        } else {
            Node_Bond_identifier_3();
            cout<<"# of bonds = "<<Num_of_Node_Pairs<<endl;
        }
    }
    else if(MeshType=="sajjad")
    {
        Node_Bond_identifier_2();
        cout<<"# of filaments = "<<Num_of_Node_Pairs<<endl;
        cout<<"# of abps = "<<Num_of_abp_Pairs<<endl;
        cout<<"# of MTs = "<<Num_of_MT_Pairs<<endl;
    }
    
    if(MeshType=="normal"){
        check();
    }
    
    Node_neighbour_list_constructor();
    
    set_bond_nominal_length();
    set_node_radius();
    
    
    if(MeshType=="sajjad")
    {
        initialise_node_bond_relaxed_length();
        initialise_abp_bond_relaxed_length();
        initialise_MT_bond_relaxed_length();
    }
    tau_Maxwell_relax = Dashpot_Viscosity/Spring_coefficient;
    exp_tau=exp(-generalParameters.Step_Size_In_Fs/tau_Maxwell_relax);
    
    
    if(MeshType=="sajjad")
    {
        check();
        check_2();
        check_3();
    }
    
    
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
    
    shift_velocity(Shift_velocities_xyzVector[0], Shift_velocities_xyzVector[1], Shift_velocities_xyzVector[2]);
    shift_position(Shift_position_xyzVector[0], Shift_position_xyzVector[1], Shift_position_xyzVector[2]);
    
    cout<<TSUCCESS<<"\nActin class initiated.\n"<<TRESET<<
    "************************\n"<<endl;
}

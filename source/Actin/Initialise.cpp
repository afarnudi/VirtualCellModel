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

void Actin::initialise(int type){
//    T_Kinetic_Energy.resize(100);
    cout<<"Initialising the Actin Class..."<<endl;
    if(type==1){
    read_gmesh_file(Mesh_file_name);
    }
    else if(type==2)
    {
       read_gmesh_file_2(Mesh_file_name);
    }
    cout<<"# of Nodes="<<Num_of_Nodes<<endl;
    if(type==1){
    Node_Bond_identifier();
        cout<<"# of bonds = "<<Num_of_Node_Pairs<<endl;
    }
    else if(type==2)
    {
        Node_Bond_identifier_2();
        cout<<"# of filaments = "<<Num_of_Node_Pairs<<endl;
        cout<<"# of abps = "<<Num_of_abp_Pairs<<endl;
        cout<<"# of MTs = "<<Num_of_MT_Pairs<<endl;
    }
    
    
    shift_node_positions();
    
    if(type==1){
    initialise_node_bond_relaxed_length();
    }
    else if(type==2)
    {
       initialise_node_bond_relaxed_length();
       initialise_abp_bond_relaxed_length();
        initialise_MT_bond_relaxed_length();
    }
    tau_Maxwell_relax=Dashpot_Viscosity/Spring_coefficient;
    exp_tau=exp(-GenConst::Step_Size_In_Fs/tau_Maxwell_relax);
    
    if(type==1){
    check();
    }
    else if(type==2)
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
    
    cout<<"\nActin class initiated.\n******************************\n\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}

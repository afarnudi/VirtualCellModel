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

void Actin::initialise(){
//    T_Kinetic_Energy.resize(100);
    cout<<"Initialising the Actin Class..."<<endl;
    read_gmesh_file(Mesh_file_name);
    cout<<"# of Nodes="<<Num_of_Nodes<<endl;
    Node_Bond_identifier();
    cout<<"# of bonds = "<<Num_of_Node_Pairs<<endl;
    
    shift_node_positions();
    
    initialise_node_bond_relaxed_length();
    tau_Maxwell_relax=Dashpot_Viscosity/Spring_coefficient;
    exp_tau=exp(-GenConst::Step_Size_In_Fs/tau_Maxwell_relax);
    
    cout<<"\nActin class initiated.\n******************************\n\n";
    //        cout<< "Average node distance is   "<<Average_Membrane_Node_Distance()<<endl;
}

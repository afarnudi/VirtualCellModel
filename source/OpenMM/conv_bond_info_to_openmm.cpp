#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"

Bonds* convert_membrane_bond_info_to_openmm(Membrane mem) {
    const int mem_num_bonds = mem.return_num_of_node_pairs();
    Bonds* bonds = new Bonds[mem_num_bonds];
    
    //used in openmm to specify different types of atoms. I don't know what the application is at the moment.
    int CC=0;
    
    for (int i=0; i<mem_num_bonds; i++) {
        bonds[i].type = CC;
        bonds[i].atoms[0]=mem.return_node_pair(i, 0);
        bonds[i].atoms[1]=mem.return_node_pair(i, 1);
        bonds[i].nominalLengthInAngstroms=mem.return_avg_node_dist();
        bonds[i].stiffnessInKcalPerAngstrom2=mem.return_spring_stiffness_coefficient();
    }
    //End of list
//    bonds[mem_num_bonds].type=-1;
    
    return bonds;
}

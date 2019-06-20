#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"

Bonds* convert_membrane_bond_info_to_openmm(Membrane mem) {
    const int mem_num_bonds = mem.get_num_of_node_pairs();
    Bonds* bonds = new Bonds[mem_num_bonds];
    for (int i=0; i<mem_num_bonds; i++) {
        bonds[i].type = mem.get_spring_model();
        bonds[i].atoms[0]=mem.get_node_pair(i, 0);
        bonds[i].atoms[1]=mem.get_node_pair(i, 1);
        switch (bonds[i].type) {
            //FENE
            case 1:
                mem.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
                bonds[i].stiffnessInKcalPerAngstrom2=mem.get_spring_stiffness_coefficient();
                break;
            //Harmonic
            case 2:
                bonds[i].nominalLengthInAngstroms=mem.get_avg_node_dist();
                bonds[i].stiffnessInKcalPerAngstrom2=mem.get_spring_stiffness_coefficient();
                break;
                
            
        }
       
        
    }
    
    return bonds;
}

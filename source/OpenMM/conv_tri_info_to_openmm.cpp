#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"

Dihedrals* convert_membrane_dihedral_info_to_openmm(Membrane mem) {
    const int mem_num_tris = mem.return_num_of_triangle();
    Dihedrals* diatoms = new Dihedrals[mem_num_tris+1];
    
    
    for (int i=0; i<mem_num_tris; i++) {
        vector<int> tri_pair_nodes(mem.return_traingle_pair_nodes(i));
        diatoms[i].atoms.resize(4);
        diatoms[i].type=0;
        diatoms[i].atoms[0]=tri_pair_nodes[0];
        diatoms[i].atoms[1]=tri_pair_nodes[1];
        diatoms[i].atoms[2]=tri_pair_nodes[2];
        diatoms[i].atoms[3]=tri_pair_nodes[3];
        diatoms[i].bending_stiffness_value = mem.return_spring_stiffness_coefficient() * 0.2 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
    }
    //End of list
    diatoms[mem_num_tris].type=-1;
    
    return diatoms;
}

#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"

Dihedrals* convert_membrane_dihedral_info_to_openmm(Membrane &mem) {
    const int mem_num_tri_pairs = mem.get_num_of_triangle_pairs();
//    cout<<"**++**++**++\n\t"<<mem_num_tri_pairs<<"\n**++**++**++\n";
    Dihedrals* diatoms = new Dihedrals[mem_num_tri_pairs];
    
    
    for (int i=0; i<mem_num_tri_pairs; i++) {
        vector<int> tri_pair_nodes(mem.get_traingle_pair_nodes_list(i));

        diatoms[i].atoms.resize(4);
        diatoms[i].type=0;
        diatoms[i].atoms[0]=tri_pair_nodes[0];
        diatoms[i].atoms[1]=tri_pair_nodes[1];
        diatoms[i].atoms[2]=tri_pair_nodes[2];
        diatoms[i].atoms[3]=tri_pair_nodes[3];

        diatoms[i].bendingStiffnessinKJ = mem.get_bending_stiffness_coefficient();
        diatoms[i].class_label = mem.get_label();
    }
    cout<<" harmonic"<<endl;
    cout<<"\tCoeficient (KJ . mol^-1 ) = "<<mem.get_bending_stiffness_coefficient() <<endl;
    
    return diatoms;
}

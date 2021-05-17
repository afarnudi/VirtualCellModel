#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"

Dihedrals* convert_membrane_dihedral_info_to_openmm(Membrane &mem) {
    
    
    bool dihedralPotential = false;
    bool noPotential = false;
    
    const int mem_num_tri_pairs = mem.get_num_of_dihedral_elements();
//    cout<<"**++**++**++\n\t"<<mem_num_tri_pairs<<"\n**++**++**++\n";
    Dihedrals* diatoms = new Dihedrals[mem_num_tri_pairs];
    
    
    for (int i=0; i<mem_num_tri_pairs; i++) {
//        vector<int> tri_pair_nodes(mem.get_traingle_pair_nodes_list(i));
        
        diatoms[i].type=mem.get_bending_model();
        if (diatoms[i].type == potentialModelIndex.Model["Dihedral"]) {
            dihedralPotential = true;
            vector<int> tri_pair_nodes(mem.get_dihedral_atoms_list(i));
            diatoms[i].atoms.resize(4);
            diatoms[i].atoms[0]=tri_pair_nodes[0];
            diatoms[i].atoms[1]=tri_pair_nodes[1];
            diatoms[i].atoms[2]=tri_pair_nodes[2];
            diatoms[i].atoms[3]=tri_pair_nodes[3];

            diatoms[i].bendingStiffnessinKJ = mem.get_bending_stiffness_coefficient();
            diatoms[i].class_label = mem.get_label();
            diatoms[i].spontaneousBendingAngleInRad = mem.get_spontaneous_angle_in_Rad(i);
        } else if (diatoms[i].type == potentialModelIndex.Model["None"]){
            noPotential = true;
        }
        
        
    }
    if (dihedralPotential) {
        cout<<" Cosine(dihedral)"<<endl;
        cout<<"\tCoeficient (KJ . mol^-1 ) = "<<mem.get_bending_stiffness_coefficient() <<endl;
    }
    
    if (noPotential || mem_num_tri_pairs==0) {
        cout<<" None"<<endl;
    }
    
    
    return diatoms;
}

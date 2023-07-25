#include "General_class_functions.h"


void Update_classes(std::vector<Membrane>   &membranes,
                    std::vector<Actin>      &actins,
                    std::vector<ECM>        &ecms,
                    std::vector<Chromatin>  &chromatins,
                    double                   time,
                    MyAtomInfo*              all_atoms)
{
    int atom_count=0;
    int bond_count=0;
    int dihe_count=0;
    
    for (int i=0; i<membranes.size(); i++) {
        membranes[i].export_for_resume(time/generalParameters.Step_Size_In_Fs, all_atoms, atom_count);
        atom_count += membranes[i].get_num_of_nodes();
        bond_count += membranes[i].get_num_of_node_pairs();
        dihe_count += membranes[i].get_num_of_triangle_pairs();
    }
    for (int i=0; i<actins.size(); i++) {
        actins[i].export_for_resume(time/generalParameters.Step_Size_In_Fs, all_atoms, atom_count);
        atom_count += actins[i].get_num_of_nodes();
        bond_count += actins[i].get_num_of_node_pairs();
    }
    for (int i=0; i<ecms.size(); i++) {
        ecms[i].export_for_resume(time/generalParameters.Step_Size_In_Fs, all_atoms, atom_count);
        atom_count += ecms[i].get_num_of_nodes();
        bond_count += ecms[i].get_num_of_node_pairs();
    }
    for (int i=0; i<chromatins.size(); i++) {
        chromatins[i].set_state(all_atoms, atom_count);
        chromatins[i].export_for_resume(time/generalParameters.Step_Size_In_Fs, all_atoms, atom_count);
        atom_count += chromatins[i].get_num_of_nodes();
        bond_count += chromatins[i].get_num_of_nodes()-1;
    }
    
    
}

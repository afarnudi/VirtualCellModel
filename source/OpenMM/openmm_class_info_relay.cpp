#include "OpenMM_funcs.hpp"

void OpenMM_membrane_info_relay (vector<Membrane>       membranes,
                                 vector<std::set<int> > &membrane_set,
                                 MyAtomInfo*            all_atoms,
                                 Bonds*                 all_bonds,
                                 Dihedrals*             all_dihedrals,
                                 int                    &atom_count,
                                 int                    &bond_count,
                                 int                    &dihe_count){
    
    for (int i=0; i<membranes.size(); i++) {
        
        //Create a set of the atom index to use for OpenMM's custom non bond interaction set.
        
        MyAtomInfo* atoms = convert_membrane_position_to_openmm(membranes[i]);
        for (int j=0;j<membranes[i].get_num_of_nodes(); j++) {
            all_atoms[j+atom_count]=atoms[j];
            membrane_set[i].insert(j+atom_count);
        }
        
        
        cout<<TMEM<<"Membrane "<<i<<TRESET<<endl;
        cout<<"Bond potential:";
        Bonds* bonds = convert_membrane_bond_info_to_openmm(membranes[i]);
        for (int j=0; j<membranes[i].get_num_of_node_pairs(); j++) {
            all_bonds[j+bond_count]=bonds[j];
            all_bonds[j+bond_count].atoms[0]=bonds[j].atoms[0]+atom_count;
            all_bonds[j+bond_count].atoms[1]=bonds[j].atoms[1]+atom_count;
            
        }
        
        cout<<"Dihedral potential:";
        Dihedrals* dihedrals = convert_membrane_dihedral_info_to_openmm(membranes[i]);
        for (int j=0; j<membranes[i].get_num_of_triangle_pairs(); j++) {
            all_dihedrals[j+dihe_count]=dihedrals[j];
            all_dihedrals[j+dihe_count].atoms[0]=dihedrals[j].atoms[0]+atom_count;
            all_dihedrals[j+dihe_count].atoms[1]=dihedrals[j].atoms[1]+atom_count;
            all_dihedrals[j+dihe_count].atoms[2]=dihedrals[j].atoms[2]+atom_count;
            all_dihedrals[j+dihe_count].atoms[3]=dihedrals[j].atoms[3]+atom_count;
        }
        
        //These parameters are used to shift the index of the atoms/bonds/dihedrals.
        atom_count += membranes[i].get_num_of_nodes();
        bond_count += membranes[i].get_num_of_node_pairs();
        dihe_count += membranes[i].get_num_of_triangle_pairs();
    }
}

void OpenMM_Actin_info_relay (vector<Actin>          acts,
                              vector<std::set<int> > &act_set,
                              MyAtomInfo*            all_atoms,
                              Bonds*                 all_bonds,
                              Dihedrals*             all_dihedrals,
                              int                    &atom_count,
                              int                    &bond_count,
                              int                    &dihe_count){
    for (int i=0; i<acts.size(); i++) {
        
        //Create a set of the atom index to use for OpenMM's custom non bond interaction set.
        cout<<TACT<<"Actin "<<i<<TRESET<<endl;
        MyAtomInfo* atoms = convert_Actin_position_to_openmm(acts[i]);
        for (int j=0;j<acts[i].get_num_of_nodes(); j++) {
            all_atoms[j+atom_count]=atoms[j];
            act_set[i].insert(j+atom_count);
        }
        
        
        cout<<"Bond potential:";
        Bonds* bonds = convert_Actin_bond_info_to_openmm(acts[i],atoms);
//        for (int j=0; j<(4*acts[i].get_num_of_node_pairs() + 4*acts[i].get_num_of_abp_pairs() +4*acts[i].get_num_of_MT_pairs() ); j++) {
        for (int j=0; j<acts[i].get_num_of_node_pairs(); j++) {
            all_bonds[j+bond_count]=bonds[j];
            all_bonds[j+bond_count].atoms[0]=bonds[j].atoms[0]+atom_count;
            all_bonds[j+bond_count].atoms[1]=bonds[j].atoms[1]+atom_count;
            
        }
        
        //These parameters are used to shift the index of the atoms/bonds/dihedrals.
        atom_count += acts[i].get_num_of_nodes();
//        bond_count += (4*acts[i].get_num_of_node_pairs() + 4*acts[i].get_num_of_abp_pairs() + 4*acts[i].get_num_of_MT_pairs() ) ;
        bond_count += acts[i].get_num_of_node_pairs();
        //        dihe_count += membranes[i].get_num_of_triangle_pairs();
    }
}




void OpenMM_ActMem_info_relay (vector<Actin>          acts,
                               vector<Membrane>       membranes,
                               Bonds*                 all_bonds,
                               int                    mem_atom_count,
                               int                    &bond_count){
    
    int mem_atom_counter = 0;
    int act_atom_counter = 0;
    
    for (int i=0; i<acts.size(); i++) {
        for (int k=0; k<membranes.size(); k++) {
            
            cout<<TACT<<"Actin "<<i<<TRESET<<" and "<<TMEM<<"Membrane "<<k<<TRESET<<" Bounded."<<endl;
            Bonds* bonds = convert_ActMem_bond_info_to_openmm(acts[i], k);
            for (int j=0; j<acts[i].return_num_of_actin_membrane_shared_nodes(k); j++) {
                all_bonds[j+bond_count]=bonds[j];
                //atom 0 == actin       atom 1 == membrane
                all_bonds[j+bond_count].atoms[0]=bonds[j].atoms[0]+mem_atom_count+act_atom_counter;
                all_bonds[j+bond_count].atoms[1]=bonds[j].atoms[1]+mem_atom_counter;
                
            }
            
            mem_atom_counter += membranes[k].get_num_of_nodes();
            bond_count += acts[i].return_num_of_actin_membrane_shared_nodes(k);
            
        }
        
        mem_atom_counter = 0;
        act_atom_counter += acts[i].get_num_of_nodes();
        
    }
}

void OpenMM_ECM_info_relay (vector<ECM>            ecms,
                            vector<std::set<int> > &ecm_set,
                            MyAtomInfo*            all_atoms,
                            Bonds*                 all_bonds,
                            Dihedrals*             all_dihedrals,
                            int                    &atom_count,
                            int                    &bond_count,
                            int                    &dihe_count){
    for (int i=0; i<ecms.size(); i++) {
        
        //Create a set of the atom index to use for OpenMM's custom non bond interaction set.
        cout<<TECM<<"ECM "<<i<<TRESET<<endl;
        MyAtomInfo* atoms = convert_ECM_position_to_openmm(ecms[i]);
        for (int j=0;j<ecms[i].get_num_of_nodes(); j++) {
            all_atoms[j+atom_count]=atoms[j];
            ecm_set[i].insert(j+atom_count);
        }
        
        
        cout<<"Bond potential:";
        Bonds* bonds = convert_ECM_bond_info_to_openmm(ecms[i] , atoms);
        
        for (int j=0; j<ecms[i].get_num_of_node_pairs(); j++) {
            all_bonds[j+bond_count]=bonds[j];
            all_bonds[j+bond_count].atoms[0]=bonds[j].atoms[0]+atom_count;
            all_bonds[j+bond_count].atoms[1]=bonds[j].atoms[1]+atom_count;
            
        }
        
        //These parameters are used to shift the index of the atoms/bonds/dihedrals.
        atom_count += ecms[i].get_num_of_nodes();
        bond_count += ecms[i].get_num_of_node_pairs();
        
        //        dihe_count += membranes[i].get_num_of_triangle_pairs();
    }
}

void OpenMM_Chromatin_info_relay (vector<Chromatin>                 chromos,
//                                  vector <std::set<int> > &chromatin_set,
                                  vector<vector <std::set<int> > > &chromatin_set,
                                  MyAtomInfo*            	        all_atoms,
                                  Bonds*                            all_bonds,
                                  Angles*                           all_angles,
                                  int                              &atom_count,
                                  int                              &bond_count,
                                  int                              &angle_count){
    for (int i=0; i<chromos.size(); i++) {
        
        //Create a set of the atom index to use for OpenMM's custom non bond interaction set.
        cout<<TCHR<<"Chromatin "<<i<<TRESET<<endl;
        MyAtomInfo* atoms = convert_Chromatin_position_to_openmm(chromos[i]);
        
        
        
        for (int j=0; j<chromos[i].get_num_of_nodes(); j++) {
            if (atoms[j].mass < 0.0000001) {
                atoms[j].vsite_atoms[0] += atom_count;
                atoms[j].vsite_atoms[1] += atom_count;
            }
            
            all_atoms[j+atom_count]=atoms[j];
            
            
            chromatin_set[i][chromos[i].get_node_type(j)].insert(j+atom_count);
        }
        
        //Check to if all node-types in the chromatin sets have a node associated to them. An empty node set will cause a problem in the interaction definition section since we use pointers to access the "set" and an empty set can may point to random numbers.
        for (int j=0; j<chromatin_set[i].size(); j++) {
            if (chromatin_set[i][j].size()==0) {
                cout<<"some node-types have no nodes associated with them. If a random distribution is used for node-type association this may be due to limited number of nodes.\n";
                exit(EXIT_FAILURE);
            }
        }
        cout<<"Bond potential:";
        Bonds* bonds = convert_Chromatin_bond_info_to_openmm(chromos[i]);
        
        for (int j=0; j<chromos[i].get_num_of_bonds(); j++) {
            all_bonds[j+bond_count]=bonds[j];
            all_bonds[j+bond_count].atoms[0]=bonds[j].atoms[0]+atom_count;
            all_bonds[j+bond_count].atoms[1]=bonds[j].atoms[1]+atom_count;
        }
        
        cout<<"Angle bond potential:";
        Angles* angles = convert_Chromatin_angle_bond_info_to_openmm(chromos[i]);
        
        for (int j=0; j<chromos[i].get_num_of_bonds(); j++) {
            all_angles[j+angle_count]=angles[j];
            all_angles[j+angle_count].atoms[0]=angles[j].atoms[0]+atom_count;
            all_angles[j+angle_count].atoms[1]=angles[j].atoms[1]+atom_count;
            all_angles[j+angle_count].atoms[2]=angles[j].atoms[2]+atom_count;
        }
        
        
        
        
        
        //These parameters are used to shift the index of the atoms/bonds/dihedrals.
        atom_count  += chromos[i].get_num_of_nodes();
        bond_count  += chromos[i].get_num_of_bonds();
        angle_count += chromos[i].get_num_of_angle_bonds();
        //        dihe_count += membranes[i].get_num_of_triangle_pairs();
    }
}

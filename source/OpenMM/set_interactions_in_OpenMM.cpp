
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include <stdlib.h>
#include <time.h>

using OpenMM::Vec3;
using std::vector;
using std::set;


void set_interactions(const MyAtomInfo                       atoms[],
                      Bonds*                                 bonds,
                      vector<std::set<int> >                &membrane_set,
                      vector<std::set<int> >                &actin_set,
                      vector<std::set<int> >                &ecm_set,
                      vector<vector<set<int> > >            &chromatin_set,
                      NonBondInteractionMap                 &interaction_map,
                      vector<OpenMM::CustomExternalForce*>  &ext_force,
                      vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                      vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                      OpenMM::System                        &system
                      ){
    
    
    vector<std::set<int>> sp_nodes;
    vector<std::set<int>> ns_nodes;
    sp_nodes.resize(1);
    ns_nodes.resize(1);
    
    
    vector<std::set<int>> ad_sites;
    vector<std::set<int>> na_sites;
    ad_sites.resize(1);
    na_sites.resize(1);
    
//    srand(time(NULL));
//
//    for(set<int>::iterator it=membrane_set[0].begin(); it != membrane_set[0].end(); ++it)
//    {
//        int r = rand() % 100 + 1;
//
//        if(atoms[*it].initPosInNm[0]>0)
//        {
//            if (r<75)
//            {
//                sp_nodes[0].insert(*it);
//            }
//            else
//            {
//                ns_nodes[0].insert(*it);
//            }
//        }
//        else
//        {
//            if (r<36)
//            {
//                sp_nodes[0].insert(*it);
//            }
//            else
//            {
//                ns_nodes[0].insert(*it);
//            }
//        }
//    }
//
//    cout<<"sp nodes"<<sp_nodes[0].size()<< '\n';
//    cout<<"ns nodes"<<ns_nodes[0].size()<< '\n';
    
    
    //srand(time(NULL));
    if (GenConst::Num_of_ECMs!=0) {
        for(set<int>::iterator it=ecm_set[0].begin(); it != ecm_set[0].end(); ++it)
        {
            
            if( atoms[*it].symbol == 'E')
            {
                ad_sites[0].insert(*it);
            }
            else
            {
                na_sites[0].insert(*it);
            }
            

        }
        
        cout<<"ad sites"<<ad_sites[0].size()<< '\n';
        cout<<"na sites"<<na_sites[0].size()<< '\n';
    }
    
    
    
    std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds);
   //Order: Membranes, Actins, ECMs, Chromatins, Point Particles
    for (int i=0; i< GenConst::Num_of_Membranes; i++) {
        //initialize external force for membrane i
        int ext_force_index;
        bool is_force = init_ext_force(ext_force, atoms, membrane_set, i , GenConst::Membrane_label);
        ext_force_index=int(ext_force.size()-1);
        if(is_force)
        {
            system.addForce(ext_force[ext_force_index]);
        }
        
        for (int j=0; j < i+1; j++) {
            
            std::string class_label_i=GenConst::Membrane_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i, j);
            interaction_map.set_force_label(i, j, GenConst::Membrane_label + std::to_string(i) +  GenConst::Membrane_label  + std::to_string(j));
            if (interaction_name == "Lennard-Jones") {
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, membrane_set, membrane_set, i, j , GenConst::Membrane_label , GenConst::Membrane_label);
                index = int(LJ_12_6_interactions.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                
                system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, membrane_set, membrane_set, i, j, GenConst::Membrane_label , GenConst::Membrane_label);
                
                index = int(ExcludedVolumes.size()-1);
                //
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(ExcludedVolumes[index], exclude_bonds);
                
                system.addForce(ExcludedVolumes[index]);
            }
            
        }
    }
    
    int class_count_i, class_count_j;
    
    for (int i=0; i < GenConst::Num_of_Actins; i++) {
        
        //initialize external force for actin i
        int ext_force_index;
        bool is_force = init_ext_force(ext_force, atoms, actin_set, i , GenConst::Actin_label);
        ext_force_index=int(ext_force.size()-1);
        if(is_force)
        {
            system.addForce(ext_force[ext_force_index]);
        }
        
        class_count_i = GenConst::Num_of_Membranes;
        class_count_j = 0;
        for (int j=0; j < GenConst::Num_of_Membranes; j++) {
            std::string class_label_i=GenConst::Actin_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i + class_count_i, j);
            interaction_map.set_force_label(i + class_count_i, j, GenConst::Actin_label + std::to_string(i) +  GenConst::Membrane_label  + std::to_string(j));
            
            if (interaction_name == "Lennard-Jones") {
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, actin_set, membrane_set, i, j , GenConst::Actin_label , GenConst::Membrane_label);
                index = int(LJ_12_6_interactions.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                
                system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, actin_set, membrane_set, i, j, GenConst::Actin_label , GenConst::Membrane_label);
                
                index = int(ExcludedVolumes.size()-1);
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(ExcludedVolumes[index], exclude_bonds);
                
                system.addForce(ExcludedVolumes[index]);
            }
        }
        
        class_count_j = GenConst::Num_of_Membranes;
        for (int j=class_count_j; j < class_count_j + i +1 ; j++) {
            
            std::string class_label_i=GenConst::Actin_label+std::to_string(i);
            std::string class_label_j=GenConst::Actin_label+std::to_string(j-class_count_j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i + class_count_i, j);
            interaction_map.set_force_label(i+ class_count_i, j, GenConst::Actin_label + std::to_string(i) +  GenConst::Actin_label  + std::to_string(j-class_count_j));
            if (interaction_name == "Lennard-Jones") {
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, actin_set, actin_set, i, j-class_count_j , GenConst::Actin_label , GenConst::Actin_label);
                index = int(LJ_12_6_interactions.size()-1);
                
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                
                system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, actin_set, actin_set, i, j-class_count_j, GenConst::Actin_label , GenConst::Actin_label);
                
                index = int(ExcludedVolumes.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(ExcludedVolumes[index], exclude_bonds);
                
                system.addForce(ExcludedVolumes[index]);
            }
            
        }
    }
    
    
    class_count_i = GenConst::Num_of_Membranes + GenConst::Num_of_Actins;
    class_count_j = 0;
    
    for (int i=0; i < GenConst::Num_of_ECMs; i++) {
        
        //initialize external force for ecm i
        int ext_force_index;
        bool is_force = init_ext_force(ext_force, atoms, ecm_set, i , GenConst::ECM_label);
        ext_force_index=int(ext_force.size()-1);
        if(is_force)
        {
            system.addForce(ext_force[ext_force_index]);
        }
        
        for (int j=0; j < GenConst::Num_of_Membranes; j++) {
            
            
            std::string class_label_i=GenConst::ECM_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i + class_count_i, j);
            interaction_map.set_force_label(i +class_count_i, j, GenConst::ECM_label + std::to_string(i) +  GenConst::Membrane_label  + std::to_string(j));
            if (interaction_name == "Lennard-Jones") {
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, membrane_set, i, j , GenConst::ECM_label , GenConst::Membrane_label);
                index = int(LJ_12_6_interactions.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                
                system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, ecm_set, membrane_set, i, j, GenConst::ECM_label , GenConst::Membrane_label);
                index = int(ExcludedVolumes.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(ExcludedVolumes[index], exclude_bonds);
                
                system.addForce(ExcludedVolumes[index]);
            } else if (interaction_name == "Lennard-Jones3"){
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, ns_nodes, i, j , GenConst::ECM_label , GenConst::Membrane_label);

                    index = int(LJ_12_6_interactions.size()-1);
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    system.addForce(LJ_12_6_interactions[index]);
                                   
                                   
                    initdouble_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, sp_nodes, i, j , GenConst::ECM_label , GenConst::Membrane_label);

                    index = int(LJ_12_6_interactions.size()-1);
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Lennard-Jones4"){
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, na_sites, membrane_set, i, j , GenConst::ECM_label , GenConst::Membrane_label);

                    index = int(LJ_12_6_interactions.size()-1);
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    system.addForce(LJ_12_6_interactions[index]);
                                   
                                   
                    initdouble_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ad_sites, membrane_set, i, j , GenConst::ECM_label , GenConst::Membrane_label);

                    index = int(LJ_12_6_interactions.size()-1);
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Lennard-Jones5"){
                init_LJ_4_2_interaction(LJ_12_6_interactions, atoms, ecm_set, membrane_set, i, j , GenConst::ECM_label , GenConst::Membrane_label);

                    index = int(LJ_12_6_interactions.size()-1);
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    system.addForce(LJ_12_6_interactions[index]);
                                   
                                   
                    initdouble_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ad_sites, membrane_set, i, j , GenConst::ECM_label , GenConst::Membrane_label);

                    index = int(LJ_12_6_interactions.size()-1);
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    system.addForce(LJ_12_6_interactions[index]);
            }
            
            
        }
        
        class_count_j = GenConst::Num_of_Membranes;
        for (int j=class_count_j; j < class_count_j+GenConst::Num_of_Actins; j++) {
            
            
            std::string class_label_i=GenConst::ECM_label+std::to_string(i);
            std::string class_label_j=GenConst::Actin_label+std::to_string(j-class_count_j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i + class_count_j, j);
            interaction_map.set_force_label(i +class_count_j, j, GenConst::ECM_label + std::to_string(i) +  GenConst::Actin_label  + std::to_string(j-class_count_j));
            if (interaction_name == "Lennard-Jones") {
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, actin_set, i, j-class_count_j, GenConst::ECM_label, GenConst::Actin_label);
                index = int(LJ_12_6_interactions.size()-1);
                
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                
                system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, ecm_set, actin_set, i, j-class_count_j, GenConst::ECM_label, GenConst::Actin_label);
                index = int(ExcludedVolumes.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(ExcludedVolumes[index], exclude_bonds);
                
                system.addForce(ExcludedVolumes[index]);
            }
            
        }
        
        
        class_count_j = GenConst::Num_of_Membranes + GenConst::Num_of_Actins;
        for (int j=class_count_j; j < class_count_j + i + 1; j++) {
            
            std::string class_label_i=GenConst::ECM_label+std::to_string(i);
            std::string class_label_j=GenConst::ECM_label+std::to_string(j-class_count_j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i + class_count_i, j);
            interaction_map.set_force_label(i + class_count_i, j, GenConst::ECM_label + std::to_string(i) +  GenConst::ECM_label  + std::to_string(j-class_count_j));
            if (interaction_name == "Lennard-Jones") {
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, ecm_set, i, j-class_count_j, GenConst::ECM_label , GenConst::ECM_label);
                index = int(LJ_12_6_interactions.size()-1);
                
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                
                system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, ecm_set, ecm_set, i, j-class_count_j, GenConst::ECM_label , GenConst::ECM_label);
                index = int(ExcludedVolumes.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                add_exclusion(ExcludedVolumes[index], exclude_bonds);
                
                system.addForce(ExcludedVolumes[index]);
            }
        }
    }
    
    class_count_i = GenConst::Num_of_Membranes + GenConst::Num_of_Actins + GenConst::Num_of_ECMs;
    class_count_j = 0;
    
    for (int i=0; i < GenConst::Num_of_Chromatins; i++) {
        for (int j=0; j < GenConst::Num_of_Membranes; j++) {
            
            std::string class_label_i=GenConst::Chromatin_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i + class_count_i, j);
            interaction_map.set_force_label(i + class_count_i, j, GenConst::Chromatin_label + std::to_string(i) +  GenConst::Membrane_label  + std::to_string(j));
            if (interaction_name == "Lennard-Jones") {
                for (int chromo_type=0; chromo_type< chromatin_set[i].size(); chromo_type++) {
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set[i][chromo_type], membrane_set, i, j, chromo_type, GenConst::Chromatin_label, GenConst::Membrane_label);
                    index = int(LJ_12_6_interactions.size()-1);
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                }
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, membrane_set, i, j, GenConst::Chromatin_label , GenConst::Membrane_label);
                index = int(ExcludedVolumes.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                
//                    ExcludedVolumes[index]->setForceGroup(7);
//                    cout<<TWWARN<<"Set forcegroup"<<TRESET<<endl;
                
                system.addForce(ExcludedVolumes[index]);
            }
            
            
            
        }
        
        class_count_j = GenConst::Num_of_Membranes;
        for (int j=class_count_j; j < class_count_j+GenConst::Num_of_Actins; j++) {
            
            
            std::string class_label_i=GenConst::Chromatin_label+std::to_string(i);
            std::string class_label_j=GenConst::Actin_label+std::to_string(j-class_count_j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i + class_count_j, j);
            interaction_map.set_force_label(i + class_count_i, j, GenConst::Chromatin_label + std::to_string(i) +  GenConst::Actin_label  + std::to_string(j-class_count_j));
            if (interaction_name == "Lennard-Jones") {
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, actin_set, i, j-class_count_j , GenConst::Chromatin_label , GenConst::Actin_label);
                index = int(LJ_12_6_interactions.size()-1);
                
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                
                system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, actin_set, i, j-class_count_j, GenConst::Chromatin_label , GenConst::Actin_label);
                
                index = int(ExcludedVolumes.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                
                
                system.addForce(ExcludedVolumes[index]);
            }
            
        }
        
        
        class_count_j = GenConst::Num_of_Membranes + GenConst::Num_of_Actins;
        for (int j=class_count_j; j < class_count_j + GenConst::Num_of_ECMs; j++) {
            
            std::string class_label_i=GenConst::Chromatin_label+std::to_string(i);
            std::string class_label_j=GenConst::ECM_label+std::to_string(j-class_count_j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i + class_count_i, j);
            interaction_map.set_force_label(i + class_count_i, j, GenConst::Chromatin_label + std::to_string(i) +  GenConst::ECM_label  + std::to_string(j-class_count_j));
            if (interaction_name == "Lennard-Jones") {
                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, ecm_set, i, j-class_count_j , GenConst::Chromatin_label , GenConst::ECM_label);
                index = int(LJ_12_6_interactions.size()-1);
                
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                
                system.addForce(LJ_12_6_interactions[index]);
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, ecm_set, i, j-class_count_j, GenConst::Chromatin_label , GenConst::ECM_label);
                index = int(ExcludedVolumes.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                
                
                system.addForce(ExcludedVolumes[index]);
            }
            
        }
        
        class_count_j = GenConst::Num_of_Membranes + GenConst::Num_of_Actins + GenConst::Num_of_ECMs;
        
        for (int j=class_count_j; j < class_count_j + i +1; j++) {
            
            std::string class_label_i=GenConst::Chromatin_label+std::to_string(i);
            std::string class_label_j=GenConst::Chromatin_label+std::to_string(j-class_count_j);
            
            //std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            string interaction_name = interaction_map.get_interacton_type(i + class_count_i, j);
            interaction_map.set_force_label(i + class_count_i, j, GenConst::Chromatin_label + std::to_string(i) +  GenConst::Chromatin_label  + std::to_string(j-class_count_j));
            if (interaction_name == "Lennard-Jones") {
                if (i == j-class_count_j) {
                    for (int chr_type_1=0; chr_type_1<chromatin_set[i].size(); chr_type_1++) {
                        for (int chr_type_2=chr_type_1; chr_type_2<chromatin_set[j-class_count_j].size(); chr_type_2++) {
                            
                            set<int> :: iterator it_1 = chromatin_set[i][chr_type_1].begin();
                            set<int> :: iterator it_2 = chromatin_set[j-class_count_j][chr_type_2].begin();
                            
                            if (atoms[*it_1].epsilon_LJ_12_6* atoms[*it_2].epsilon_LJ_12_6 != 0){
                                
                                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, chromatin_set, i, j-class_count_j, chr_type_1, chr_type_2, GenConst::Chromatin_label , GenConst::Membrane_label);
                                index = int(LJ_12_6_interactions.size()-1);
                                
                                // Add the list of atom pairs that are excluded from the excluded volume force.
                                LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                                
                                system.addForce(LJ_12_6_interactions[index]);
                            } else {
                                init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, chromatin_set, i, j-class_count_j, chr_type_1, chr_type_2, GenConst::Chromatin_label , GenConst::Chromatin_label);
                                index = int(ExcludedVolumes.size()-1);
                                // Add the list of atom pairs that are excluded from the excluded volume force.
                                    ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                                
                                system.addForce(ExcludedVolumes[index]);
                            }
                        }
                    }
                    
                    
                } else {
                    for (int chr_type_1=0; chr_type_1<chromatin_set[i].size(); chr_type_1++) {
                        for (int chr_type_2=0; chr_type_2<chromatin_set[j-class_count_j].size(); chr_type_2++) {
                            
                            set<int> :: iterator it_1 = chromatin_set[i][chr_type_1].begin();
                            set<int> :: iterator it_2 = chromatin_set[j-class_count_j][chr_type_2].begin();
                            
                            if (atoms[*it_1].epsilon_LJ_12_6* atoms[*it_2].epsilon_LJ_12_6 != 0){
                                
                                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, chromatin_set, i, j-class_count_j, chr_type_1, chr_type_2, GenConst::Chromatin_label , GenConst::Membrane_label);
                                index = int(LJ_12_6_interactions.size()-1);
                                
                                // Add the list of atom pairs that are excluded from the excluded volume force.
                                LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                                
                                system.addForce(LJ_12_6_interactions[index]);
                            } else {
                                init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, chromatin_set, i, j-class_count_j, chr_type_1, chr_type_2, GenConst::Chromatin_label , GenConst::Chromatin_label);
                                index = int(ExcludedVolumes.size()-1);
                                // Add the list of atom pairs that are excluded from the excluded volume force.
                                    ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                                
                                system.addForce(ExcludedVolumes[index]);
                            }
                        }
                    }
                }
            } else if (interaction_name == "Excluded-Volume"){
                init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, chromatin_set, i, j-class_count_j, GenConst::Chromatin_label , GenConst::Chromatin_label);
                index = int(ExcludedVolumes.size()-1);
                
                // Add the list of atom pairs that are excluded from the excluded volume force.
                ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                
                system.addForce(ExcludedVolumes[index]);
            } else if (interaction_name == "Lennard-JonesChromatinSpecial0"){
                if (i == j-class_count_j) {
                    for (int chr_type_1=0; chr_type_1<chromatin_set[i].size(); chr_type_1++) {
                        for (int chr_type_2=chr_type_1; chr_type_2<chromatin_set[j-class_count_j].size(); chr_type_2++) {
                            
                            set<int> :: iterator it_1 = chromatin_set[i][chr_type_1].begin();
                            set<int> :: iterator it_2 = chromatin_set[j-class_count_j][chr_type_2].begin();
                            
                            if (atoms[*it_1].epsilon_LJ_12_6* atoms[*it_2].epsilon_LJ_12_6 != 0){
                                
                                init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, chromatin_set, i, j-class_count_j, chr_type_1, chr_type_2, GenConst::Chromatin_label , GenConst::Membrane_label);
                                index = int(LJ_12_6_interactions.size()-1);
                                
                                // Add the list of atom pairs that are excluded from the excluded volume force.
                                LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                                
                                system.addForce(LJ_12_6_interactions[index]);
                            } else {
                                init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, chromatin_set, i, j-class_count_j, chr_type_1, chr_type_2, GenConst::Chromatin_label , GenConst::Chromatin_label);
                                index = int(ExcludedVolumes.size()-1);
                                // Add the list of atom pairs that are excluded from the excluded volume force.
                                    ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                                
                                system.addForce(ExcludedVolumes[index]);
                            }
                        }
                    }
                    
                    
                } else {
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, chromatin_set, i, j-class_count_j, GenConst::Chromatin_label , GenConst::Chromatin_label);
                    index = int(ExcludedVolumes.size()-1);
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(ExcludedVolumes[index]);
                }
            }
            
        }
    }
    
}


#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;


void set_interactions(const MyAtomInfo                       atoms[],
                      Bonds*                                 bonds,
                      vector<std::set<int> >                &membrane_set,
                      vector<std::set<int> >                &actin_set,
                      vector<std::set<int> >                &ecm_set,
                      vector<vector<set<int> > >            &chromatin_set,
                      vector<vector<int> >                   interaction_map,
                      vector<OpenMM::CustomExternalForce*>  &ext_force,
                      vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                      vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                      OpenMM::System                        &system
                      ){
   //Order: Membranes, Actins, ECMs, Chromatins, Point Particles
    for (int i=0; i< GenConst::Num_of_Membranes; i++) {
        //initialize external force for membrane i
        int ext_force_index;
        bool is_force = init_ext_force(ext_force, atoms, membrane_set, i , GenConst::Membrane_label);
        ext_force_index=ext_force.size()-1;
        if(is_force)
        {
            system.addForce(ext_force[ext_force_index]);
        }
        
        for (int j=0; j < i+1; j++) {
            
            std::string class_label_i=GenConst::Membrane_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            
            switch (interaction_map[i][j]) {
                
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, membrane_set, membrane_set, i, j , GenConst::Membrane_label , GenConst::Membrane_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        LJ_12_6_interactions[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, membrane_set, membrane_set, i, j, GenConst::Membrane_label , GenConst::Membrane_label);
                    
                    index = ExcludedVolumes.size()-1;
//                    cout<<"EV index = "<<index<<endl;
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        ExcludedVolumes[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    add_exclusion(ExcludedVolumes[index], exclude_bonds);
                    
                    system.addForce(ExcludedVolumes[index]);
                    
                    break;
                    
            }
        }
    }
    
    int class_count_i, class_count_j;
    
    for (int i=0; i < GenConst::Num_of_Actins; i++) {
        
        //initialize external force for actin i
        int ext_force_index;
        bool is_force = init_ext_force(ext_force, atoms, actin_set, i , GenConst::Actin_label);
        ext_force_index=ext_force.size()-1;
        if(is_force)
        {
            system.addForce(ext_force[ext_force_index]);
        }
        
        class_count_i = GenConst::Num_of_Membranes;
        class_count_j = 0;
        for (int j=0; j < GenConst::Num_of_Membranes; j++) {
            std::string class_label_i=GenConst::Actin_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            
            switch (interaction_map[i + class_count_i][j]) {
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, actin_set, membrane_set, i, j , GenConst::Actin_label , GenConst::Membrane_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        LJ_12_6_interactions[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, actin_set, membrane_set, i, j, GenConst::Actin_label , GenConst::Membrane_label);
                    
                    index = ExcludedVolumes.size()-1;
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        ExcludedVolumes[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(ExcludedVolumes[index], exclude_bonds);
                    
                    system.addForce(ExcludedVolumes[index]);
                    break;
            }
            
        }
        
        class_count_j = GenConst::Num_of_Membranes;
        for (int j=class_count_j; j < class_count_j + i +1 ; j++) {
            
            std::string class_label_i=GenConst::Actin_label+std::to_string(i);
            std::string class_label_j=GenConst::Actin_label+std::to_string(j-class_count_j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            
            
            switch (interaction_map[i+ class_count_i][j]) {
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, actin_set, actin_set, i, j-class_count_j , GenConst::Actin_label , GenConst::Actin_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        LJ_12_6_interactions[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, actin_set, actin_set, i, j-class_count_j, GenConst::Actin_label , GenConst::Actin_label);
                    
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        ExcludedVolumes[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(ExcludedVolumes[index], exclude_bonds);
                    
                    system.addForce(ExcludedVolumes[index]);
                    
                    break;
            }
        }
    }
    
    
    class_count_i = GenConst::Num_of_Membranes + GenConst::Num_of_Actins;
    class_count_j = 0;
    
    for (int i=0; i < GenConst::Num_of_ECMs; i++) {
        
        //initialize external force for ecm i
        int ext_force_index;
        bool is_force = init_ext_force(ext_force, atoms, ecm_set, i , GenConst::ECM_label);
        ext_force_index=ext_force.size()-1;
        if(is_force)
        {
            system.addForce(ext_force[ext_force_index]);
        }
        
        for (int j=0; j < GenConst::Num_of_Membranes; j++) {
            
            
            std::string class_label_i=GenConst::ECM_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            
            
            switch (interaction_map[i + class_count_i][j]) {
                
                case 1:
                    
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, membrane_set, i, j , GenConst::ECM_label , GenConst::Membrane_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        LJ_12_6_interactions[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, ecm_set, membrane_set, i, j, GenConst::ECM_label , GenConst::Membrane_label);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        ExcludedVolumes[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(ExcludedVolumes[index], exclude_bonds);
                    
                    system.addForce(ExcludedVolumes[index]);
                    break;
            }
            
        }
        
        class_count_j = GenConst::Num_of_Membranes;
        for (int j=class_count_j; j < class_count_j+GenConst::Num_of_Actins; j++) {
            
            
            std::string class_label_i=GenConst::ECM_label+std::to_string(i);
            std::string class_label_j=GenConst::Actin_label+std::to_string(j-class_count_j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            
            
            switch (interaction_map[i+ class_count_j][j]) {
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, actin_set, i, j-class_count_j, GenConst::ECM_label, GenConst::Actin_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        LJ_12_6_interactions[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, ecm_set, actin_set, i, j-class_count_j, GenConst::ECM_label, GenConst::Actin_label);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        ExcludedVolumes[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(ExcludedVolumes[index], exclude_bonds);
                    
                    system.addForce(ExcludedVolumes[index]);
                    
                    
                    break;
            }
        }
        
        
        class_count_j = GenConst::Num_of_Membranes + GenConst::Num_of_Actins;
        for (int j=class_count_j; j < class_count_j + i + 1; j++) {
            
            std::string class_label_i=GenConst::ECM_label+std::to_string(i);
            std::string class_label_j=GenConst::ECM_label+std::to_string(j-class_count_j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            
            
            switch (interaction_map[i + class_count_i][j]) {
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, ecm_set, i, j-class_count_j, GenConst::ECM_label , GenConst::ECM_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        LJ_12_6_interactions[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(LJ_12_6_interactions[index], exclude_bonds);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, ecm_set, ecm_set, i, j-class_count_j, GenConst::ECM_label , GenConst::ECM_label);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    //ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    /*for (int n=0; n<exclude_bonds.size(); ++n)
                    {
                        ExcludedVolumes[index]->addExclusion(exclude_bonds[n].first, exclude_bonds[n].second);
                    }*/
                    
                    add_exclusion(ExcludedVolumes[index], exclude_bonds);
                    
                    system.addForce(ExcludedVolumes[index]);
                    
                    
                    break;
            }
        }
    }
    
    class_count_i = GenConst::Num_of_Membranes + GenConst::Num_of_Actins + GenConst::Num_of_ECMs;
    class_count_j = 0;
    
    for (int i=0; i < GenConst::Num_of_Chromatins; i++) {
        
        for (int j=0; j < GenConst::Num_of_Membranes; j++) {
            
            std::string class_label_i=GenConst::Chromatin_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            
            
            switch (interaction_map[i + class_count_i][j]) {
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, membrane_set, i, j , GenConst::Chromatin_label , GenConst::Membrane_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, membrane_set, i, j, GenConst::Chromatin_label , GenConst::Membrane_label);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    
                    system.addForce(ExcludedVolumes[index]);
                    break;
            }
            
        }
        
        class_count_j = GenConst::Num_of_Membranes;
        for (int j=class_count_j; j < class_count_j+GenConst::Num_of_Actins; j++) {
            
            
            std::string class_label_i=GenConst::Chromatin_label+std::to_string(i);
            std::string class_label_j=GenConst::Actin_label+std::to_string(j-class_count_j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            
            
            switch (interaction_map[i+ class_count_j][j]) {
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, actin_set, i, j-class_count_j , GenConst::Chromatin_label , GenConst::Actin_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, actin_set, i, j-class_count_j, GenConst::Chromatin_label , GenConst::Actin_label);
                    
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    
                    system.addForce(ExcludedVolumes[index]);
                    
                    
                    break;
            }
        }
        
        
        class_count_j = GenConst::Num_of_Membranes + GenConst::Num_of_Actins;
        for (int j=class_count_j; j < class_count_j + GenConst::Num_of_ECMs; j++) {
            
            std::string class_label_i=GenConst::Chromatin_label+std::to_string(i);
            std::string class_label_j=GenConst::ECM_label+std::to_string(j-class_count_j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            
            
            switch (interaction_map[i + class_count_i][j]) {
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, ecm_set, i, j-class_count_j , GenConst::Chromatin_label , GenConst::ECM_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, ecm_set, i, j-class_count_j, GenConst::Chromatin_label , GenConst::ECM_label);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    
                    system.addForce(ExcludedVolumes[index]);
                    
                    
                    break;
            }
        }
        
        class_count_j = GenConst::Num_of_Membranes + GenConst::Num_of_Actins + GenConst::Num_of_ECMs;
        
        for (int j=class_count_j; j < class_count_j + i +1; j++) {
            
            std::string class_label_i=GenConst::Chromatin_label+std::to_string(i);
            std::string class_label_j=GenConst::Chromatin_label+std::to_string(j-class_count_j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            
            int index;
            
            switch (interaction_map[i + class_count_i][j]) {
                    
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, chromatin_set, i, j-class_count_j , GenConst::Chromatin_label , GenConst::Chromatin_label);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, chromatin_set, i, j-class_count_j, GenConst::Chromatin_label , GenConst::Chromatin_label);
                    index = ExcludedVolumes.size()-1;
//                    cout<<"EV index = "<<index<<endl;
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    
                    system.addForce(ExcludedVolumes[index]);
                    cout<<"exclusion complete\n";
                    
                    break;
            }
        }
    }
    
}
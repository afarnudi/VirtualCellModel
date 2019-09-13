#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"

const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;
/** -----------------------------------------------------------------------------
 *                      INITIALISE OpenMM DATA STRUCTURES
 * -----------------------------------------------------------------------------
 * We take these actions here:
 * (1) Load any available OpenMM plugins, e.g. Cuda and Brook.
 * (2) Allocate a MyOpenMMData structure to hang on to OpenMM data structures
 *     in a manner which is opaque to the caller.
 * (3) Fill the OpenMM::System with the force field parameters we want to
 *     use and the particular set of atoms to be simulated.
 * (4) Create an Integrator and a Context associating the Integrator with
 *     the System.
 * (5) Select the OpenMM platform to be used.
 * (6) Return the MyOpenMMData struct and the name of the Platform in use.
 *
 * Note that this function must understand the calling MD code's molecule and
 * force field data structures so will need to be customized for each MD code.
 */
MyOpenMMData* myInitializeOpenMM(const MyAtomInfo       atoms[],
                                 double                 stepSizeInFs,
                                 std::string&           platformName,
                                 TimeDependantData*     tdd,
                                 Bonds*                 bonds,
                                 Dihedrals*             dihedrals,
                                 vector<set<int> >      &membrane_set,
                                 vector<set<int> >      &actin_set,
                                 vector<set<int> >      &ecm_set,
                                 vector<set<int> >      &chromatin_set,
                                 vector<vector<int> >   interaction_map)
{
    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory
    (OpenMM::Platform::getDefaultPluginsDirectory());
    
    // Allocate space to hold OpenMM objects while we're using them.
    MyOpenMMData* omm = new MyOpenMMData();
    
    // Create a System and Force objects within the System.
    OpenMM::System&                 system      = *(omm->system = new OpenMM::System());
    
    // Retain a reference to each force object so we can fill in the forces.
    // Note: the System owns the force objects and will take care of deleting them;
    // don't do it yourself!
    
    
    
    cout<<"Defining interactions...\n";
    
    // Create an array of harmonic spring force objects to add to the system.
    bool HarmonicBondForce=false;
    OpenMM::HarmonicBondForce*      HarmonicBond = new OpenMM::HarmonicBondForce();
    
    bool Kelvin_VoigtBondForce=false;
    OpenMM::HarmonicBondForce*      Kelvin_VoigtBond = new OpenMM::HarmonicBondForce();
    
    // Create a vector of handles for the force objects. These handles will be added to the system. Each handle in the list will be associated with a class instance.
    
    vector<OpenMM::CustomBondForce*>X4harmonics;
    
    vector<OpenMM::CustomBondForce*> FENEs;

    vector<OpenMM::CustomNonbondedForce*> ExcludedVolumes;
    
    vector<OpenMM::CustomNonbondedForce*> LJ_12_6_interactions;

    std::vector< std::pair< int, int > > excluded_bonds;
    
    vector<OpenMM::CustomCompoundBondForce*> DihedralForces;
    
    
    //Order: Membranes, Actins, ECMs, Chromatins, Point Particles
    
    for (int i=0; i< GenConst::Num_of_Membranes; i++) {
        for (int j=0; j < i+1; j++) {
            
            std::string class_label_i=GenConst::Membrane_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            
            switch (interaction_map[i][j]) {
                
                case 1:
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, membrane_set, membrane_set, i, j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, membrane_set, membrane_set, i, j);
                    
                    index = ExcludedVolumes.size()-1;
//                    cout<<"EV index = "<<index<<endl;
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    system.addForce(ExcludedVolumes[index]);
                    
                    break;
                    
            }
        }
    }
    
    int class_count_i, class_count_j;
    
    for (int i=0; i < GenConst::Num_of_Actins; i++) {
        
        class_count_i = GenConst::Num_of_Membranes;
        class_count_j = 0;
        for (int j=0; j < GenConst::Num_of_Membranes; j++) {
            std::string class_label_i=GenConst::Actin_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            
            switch (interaction_map[i + class_count_i][j]) {
                case 1:
                    
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, actin_set, membrane_set, i, j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, actin_set, membrane_set, i, j);
                    
                    index = ExcludedVolumes.size()-1;
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    
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
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, actin_set, actin_set, i, j-class_count_j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, actin_set, actin_set, i, j-class_count_j);
                    
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    system.addForce(ExcludedVolumes[index]);
                    
                    break;
            }
        }
    }
    
    
    class_count_i = GenConst::Num_of_Membranes + GenConst::Num_of_Actins;
    class_count_j = 0;
    
    for (int i=0; i < GenConst::Num_of_ECMs; i++) {
        
        for (int j=0; j < GenConst::Num_of_Membranes; j++) {
            
            std::string class_label_i=GenConst::ECM_label+std::to_string(i);
            std::string class_label_j=GenConst::Membrane_label+std::to_string(j);
            
            std::vector< std::pair< int, int > > exclude_bonds=exclusion_list_generator(bonds, class_label_i, class_label_j);
            
            int index;
            
            
            switch (interaction_map[i + class_count_i][j]) {
                case 1:
                    
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, membrane_set, i, j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, ecm_set, membrane_set, i, j);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    
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
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, actin_set, i, j-class_count_j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, ecm_set, actin_set, i, j-class_count_j);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    
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
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, ecm_set, ecm_set, i, j-class_count_j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, ecm_set, ecm_set, i, j-class_count_j);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    
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
                    
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, membrane_set, i, j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, membrane_set, i, j);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    
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
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, actin_set, i, j-class_count_j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, actin_set, i, j-class_count_j);
                    
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    
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
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, ecm_set, i, j-class_count_j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, ecm_set, i, j-class_count_j);
                    index = ExcludedVolumes.size()-1;
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    
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
                    init_LJ_12_6_interaction(LJ_12_6_interactions, atoms, chromatin_set, chromatin_set, i, j-class_count_j);
                    index = LJ_12_6_interactions.size()-1;
                    
                    
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    LJ_12_6_interactions[index]->createExclusionsFromBonds(exclude_bonds, 0);
                    
                    system.addForce(LJ_12_6_interactions[index]);
                    break;
                case 2:
                    init_Excluded_volume_interaction(ExcludedVolumes, atoms, chromatin_set, chromatin_set, i, j-class_count_j);
                    index = ExcludedVolumes.size()-1;
//                    cout<<"EV index = "<<index<<endl;
                    // Add the list of atom pairs that are excluded from the excluded volume force.
                    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
                    ExcludedVolumes[index]->createExclusionsFromBonds(excluded_bonds, 0);
                    
                    
                    system.addForce(ExcludedVolumes[index]);
                    
                    
                    break;
            }
        }
    }
    
    
    
    
    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    std::vector<Vec3> initialPosInNm;
    std::vector<Vec3> initialVelInNmperPs;
    std::vector<double> sigma_ev;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        //        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atoms[n].mass);
        // Convert the initial position to nm and append to the array.
        const Vec3 posInNm(atoms[n].initPosInAng[0] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[1] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[2] * OpenMM::NmPerAngstrom);
        const Vec3 velocityInAngperPs(atoms[n].velocityInAngperPs[0],
                                      atoms[n].velocityInAngperPs[1],
                                      atoms[n].velocityInAngperPs[2]);
//        cout<<atoms[n].velocityInAngperPs[0]<<"\t"<<
//              atoms[n].velocityInAngperPs[1]<<"\t"<<
//              atoms[n].velocityInAngperPs[2]<<endl;
        initialPosInNm.push_back(posInNm);
        initialVelInNmperPs.push_back(velocityInAngperPs);
        
        //add particles to the excluded volume force. The number of particles should be equal to the number particles in the system. The exluded interaction lists should be defined afterwards.
        std::vector<double> sigma_ev;
        sigma_ev.push_back(atoms[n].radius
                           * OpenMM::NmPerAngstrom);
        for (int i=0; i<ExcludedVolumes.size(); i++) {
            ExcludedVolumes[i]->addParticle(sigma_ev);
        }
        for (int i=0; i<LJ_12_6_interactions.size(); i++) {
            LJ_12_6_interactions[i]->addParticle();
        }
        
        
        
    }
    
    
    set <std::string> FENE_classes;
    int FENE_index = -1;
    
    
    set <std::string> X4harmonic_classes;
    int X4harmonic_index = -1;
    
    for (int i=0; bonds[i].type != EndOfList; ++i) {
        const int*      atom = bonds[i].atoms;
        std::pair< int, int > temp;
        temp.first=bonds[i].atoms[0];
        temp.second=bonds[i].atoms[1];
        excluded_bonds.push_back(temp);
        
        
        
        
        
        
        switch (bonds[i].type) {
            case 1://FENE
            {
                auto FENE_item = FENE_classes.find(bonds[i].class_label);
                if (FENE_item == FENE_classes.end()) {
                    
                    FENE_classes.insert(bonds[i].class_label);
                    FENE_index++;
                    
                    FENEs.push_back(new OpenMM::CustomBondForce("(k_bond/(r-lmin))"));
                    FENEs[FENE_index]->addGlobalParameter("lmin",   bonds[i].FENE_lmin
                                                          * OpenMM::NmPerAngstrom);
                    FENEs[FENE_index]->addGlobalParameter("le0",   bonds[i].FENE_le0//);
                                                          * OpenMM::NmPerAngstrom);
//                    FENEs[FENE_index]->addGlobalParameter("le1",   bonds[i].FENE_le1
//                                                          * OpenMM::NmPerAngstrom);
                    FENEs[FENE_index]->addGlobalParameter("le1",   1.7
                                                          * OpenMM::NmPerAngstrom);

                    //    cout<<"bonds[0].FENE_lmax = "<<bonds[0].FENE_lmax<<endl;
                    FENEs[FENE_index]->addGlobalParameter("lmax",   bonds[i].FENE_lmax//);
                                                          * OpenMM::AngstromsPerNm);
                    FENEs[FENE_index]->addGlobalParameter("k_bond", bonds[i].stiffnessInKcalPerAngstrom2
                                                          * OpenMM::KJPerKcal//);
                                                          * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                    system.addForce(FENEs[FENE_index]);
                }
                FENEs[FENE_index]->addBond(atom[0], atom[1]);
            }
                break;
            case 2://Harmonic
            {
                HarmonicBondForce=true;
                // Note factor of 2 for stiffness below because Amber specifies the constant
                // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
                // it as used in the force term kx, with energy kx^2/2.
                HarmonicBond->addBond(atom[0], atom[1],
                                      bonds[i].nominalLengthInAngstroms
                                      * OpenMM::NmPerAngstrom,
                                      bonds[i].stiffnessInKcalPerAngstrom2
                                      * OpenMM::KJPerKcal
                                      * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                
            }
                break;
            case 3:// X4harmonic
            {
                auto X4harmonic_item = X4harmonic_classes.find(bonds[i].class_label);
                if (X4harmonic_item == X4harmonic_classes.end()) {
                    
                    X4harmonic_classes.insert(bonds[i].class_label);
                    X4harmonic_index++;
                    
                    X4harmonics.push_back( new OpenMM::CustomBondForce("0.25*k_bond*((r/r_rest)-1)^4"));
                    X4harmonics[X4harmonic_index]->addGlobalParameter("r_rest",   bonds[i].nominalLengthInAngstroms
                                                                      * OpenMM::NmPerAngstrom);
                    X4harmonics[X4harmonic_index]->addGlobalParameter("k_bond", bonds[i].stiffnessInKcalPerAngstrom4
                                                                      * OpenMM::KJPerKcal//);
                                                                      * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm*  OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                    system.addForce(X4harmonics[X4harmonic_index]);
                }
                X4harmonics[X4harmonic_index]->addBond(atom[0], atom[1]);
            }
                break;
                
            
            case 4://Kelvin-Voigt
            {
                //omm->Voigt = true;
                Kelvin_VoigtBondForce=true;
                // Note factor of 2 for stiffness below because Amber specifies the constant
                // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
                // it as used in the force term kx, with energy kx^2/2.
                Kelvin_VoigtBond->addBond(atom[0], atom[1],
                                   bonds[i].nominalLengthInAngstroms
                                   * OpenMM::NmPerAngstrom,
                                   bonds[i].stiffnessInKcalPerAngstrom2
                                   * OpenMM::KJPerKcal
                                   * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                
            }
                break;
                
            
                
        }
        
        
    }
    
    if (HarmonicBondForce) {
        system.addForce(HarmonicBond);
        omm->harmonic = HarmonicBond;
    }
    
    omm->EV = ExcludedVolumes;
    
    if (Kelvin_VoigtBondForce) {
        system.addForce(Kelvin_VoigtBond);
        tdd->Kelvin_VoigtBond = Kelvin_VoigtBond;
        tdd->Kelvin_Voigt = true;
        tdd->Kelvin_Nominal_length_calc();
    }
    
    // Add the list of atom pairs that are excluded from the excluded volume force.
    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
    
    
    
    //DFs = DihedralForces
    set <std::string> DFs_classes;
    int DFs_index = -1;
    
    //    cout<<"Here!!!!\n\n";
    
    for (int i=0; dihedrals[i].type != EndOfList; ++i) {
        
        auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
        if (DFs_item == DFs_classes.end()) {
            
            DFs_classes.insert(dihedrals[i].class_label);
            DFs_index++;
            
            DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "K_bend*(cos(dihedral(p1,p2,p3,p4)*0.5))"));
            DihedralForces[DFs_index]->addGlobalParameter("K_bend", dihedrals[i].bending_stiffness_value
                                                          * OpenMM::KJPerKcal
                                                          * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
            system.addForce(DihedralForces[DFs_index]);
        }
        DihedralForces[DFs_index]->addBond(dihedrals[i].atoms);
    }
    
    
    //Listing the names of all available platforms.
    cout<<"OpenMM available platforms:\nPlatform name  Estimated speed\n";
    for (int i = 0; i < OpenMM::Platform::getNumPlatforms(); i++) {
        OpenMM::Platform& platform = OpenMM::Platform::getPlatform(i);
        cout<<std::to_string(i)<<"- "<<platform.getName().c_str()<<"\t\t"<<platform.getSpeed()<<endl;
    }
    int platform_id=0;
    cout<<"Please choose a pltform (index): \n";
    std::cin>>platform_id;
    OpenMM::Platform& platform = OpenMM::Platform::getPlatform(platform_id);
    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero.
    
    
    switch (GenConst::Integrator_type) {
        case 0:
            omm->integrator = new OpenMM::VerletIntegrator(stepSizeInFs * OpenMM::PsPerFs);
            break;
            
        case 1:
            omm->integrator = new OpenMM::BrownianIntegrator(GenConst::temperature,
                                                             GenConst::frictionInPs,
                                                             stepSizeInFs * OpenMM::PsPerFs);
            break;
        case 2:
            omm->integrator = new OpenMM::LangevinIntegrator(GenConst::temperature,
                                                             GenConst::frictionInPs,
                                                             stepSizeInFs * OpenMM::PsPerFs);
            break;
    }
    
    

    omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform);
    omm->context->setPositions(initialPosInNm);
    omm->context->setVelocities(initialVelInNmperPs);
    platformName = omm->context->getPlatform().getName();
    
    return omm;
}

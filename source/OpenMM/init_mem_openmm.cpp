#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"

const int EndOfList=-1;
using OpenMM::Vec3;
using std::vector;
using std::set;
/** -----------------------------------------------------------------------------
 *                      INITIALIZE OpenMM DATA STRUCTURES
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
                                 Bonds*                 bonds,
                                 Dihedrals*             dihedrals,
                                 vector<set<int> >      &membrane_set,
                                 vector<set<int> >      &ecm_set,
                                 vector<vector<int> >   interaction_map)
{
    
    //======================================! Begin !======================================
    //                                  Force definitions
    //=====================================================================================
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
    
    //Define forces
    cout<<"Defining interactions...\n";
    // Create an array of harmonic spring force objects to add to the system.
    
    
    OpenMM::HarmonicBondForce*      HarmonicBond = new OpenMM::HarmonicBondForce();
   // system.addForce(HarmonicBond);
    vector<OpenMM::CustomBondForce*>X4harmonics;
    // Create a vector of handles for the FENE spring force object. These handles will be added to the system. Each handle in the list will be associated with a class instance.
    vector<OpenMM::CustomBondForce*> FENEs;
    
    
//    vector<OpenMM::CustomNonbondedForce*> ExcludedVolumes;
    OpenMM::CustomNonbondedForce* excluded_volume = new OpenMM::CustomNonbondedForce("10*(sigma/r)^6");
    
    // Creat a vector of node index pairs that form a bond. This list is later used to create a list of atoms that are excluded from the excluded volume force.
    std::vector< std::pair< int, int > > excluded_bonds;
    
    //excluded volume interaction strength.
    excluded_volume->addGlobalParameter("sigma",   atoms[0].radius
                                        * OpenMM::NmPerAngstrom);
    excluded_volume->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    excluded_volume->setCutoffDistance(2* atoms[0].radius
                                       * OpenMM::NmPerAngstrom);
//    cout<<"GenConst::Excluded_volume_interaction :: "<<GenConst::Excluded_volume_interaction<<endl;
    if (GenConst::Excluded_volume_interaction) {
//        cout<<"Ooops!\n\n";
        system.addForce(excluded_volume);
    }
    
    
    // Create a handle for 12 6 LJ inter class object interactions to add to the system.
    OpenMM::CustomNonbondedForce* LJ_12_6_interaction = new OpenMM::CustomNonbondedForce("epsilon*((sigma/r)^12-2*(sigma/r)^6)");
//    cout<<"GenConst::sigma_LJ_12_6  "<<GenConst::sigma_LJ_12_6<<endl;
    LJ_12_6_interaction->addGlobalParameter("sigma",   GenConst::sigma_LJ_12_6
                                                       * OpenMM::NmPerAngstrom);
    LJ_12_6_interaction->addGlobalParameter("epsilon",   GenConst::epsilon_LJ_12_6
                                                         * OpenMM::KJPerKcal
                                                         * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    LJ_12_6_interaction->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    LJ_12_6_interaction->setCutoffDistance(2 * GenConst::sigma_LJ_12_6
                                             * OpenMM::NmPerAngstrom);
   // system.addForce(LJ_12_6_interaction);
    
    
    vector<OpenMM::CustomCompoundBondForce*> DihedralForces;
    // We use the OpenMM's "CustomCompoundBondForce" to difine the bending force between two neighbouring membrane trinagles.
//    OpenMM::CustomCompoundBondForce* dihedral_force = new OpenMM::CustomCompoundBondForce(4, "K_bend*(cos(dihedral(p1,p2,p3,p4)))");
//    dihedral_force->addGlobalParameter("K_bend", dihedrals[0].bending_stiffness_value
//                                       * OpenMM::KJPerKcal
//                                       * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
//    system.addForce(dihedral_force);
    
    //======================================! End !======================================
    //                                  Force definitions
    //===================================================================================
    
    
    int interaction_count=0;
    //Order: Membranes, Actins, ECMs, Chromatins, Point Particles
    if (GenConst::Num_of_Membranes !=0) {
        
        for (int i=0; i< GenConst::Num_of_Membranes; i++) {
            
            for (int j=0; j < i+1; j++) {
                
                switch (interaction_map[i][j]) {
                    case 1:
//                        cout<<"mem-mem: LJ_12_6_interaction->(membrane_set[i], membrane_set[j]) = "<<i<<"\t"<<j<<endl;
                        LJ_12_6_interaction->addInteractionGroup(membrane_set[i], membrane_set[j]);
                        break;
                    case 2:
                        excluded_volume->addInteractionGroup(membrane_set[i], membrane_set[j]);
//                        cout<<"mem-mem: excluded_volume->(membrane_set[i], membrane_set[j]) = "<<i<<"\t"<<j<<endl;
                        break;
                    default:
                        break;
                }
            }
        }
    }
    interaction_count += GenConst::Num_of_Membranes;
    
    if (GenConst::Num_of_ECMs !=0) {
        
        for (int i=0; i < GenConst::Num_of_ECMs; i++) {
            int class_count=0;
            for (int j=0; j < GenConst::Num_of_Membranes; j++) {
//                cout<<"interaction_map[i + interaction_count][j]"<<interaction_map[i + interaction_count][j]<<endl;
                switch (interaction_map[i + interaction_count][j]) {
                    case 1:
//                        cout<<"ecm-mem: LJ_12_6_interaction->(ecm_set[i], membrane_set[j]) = "<<i<<"\t"<<j<<endl;
                        LJ_12_6_interaction->addInteractionGroup(ecm_set[i], membrane_set[j]);
                        break;
                    case 2:
//                        cout<<"ecm-mem: excluded_volume->(ecm_set[i], membrane_set[j]) = "<<i<<"\t"<<j<<endl;
                        excluded_volume->addInteractionGroup(ecm_set[i], membrane_set[j]);
                        break;
                    default:
                        break;
                }
            }
            
            class_count += GenConst::Num_of_Membranes;
            for (int j=class_count; j < class_count+GenConst::Num_of_ECMs; j++) {
//                cout<<"interaction_map[i + interaction_count][j]"<<interaction_map[i + interaction_count][j]<<endl;
                switch (interaction_map[i+ interaction_count][j]) {
                    case 1:
//                        cout<<"ecm-ecm: LJ_12_6_interaction->(ecm_set[i], ecm_set[j]) = "<<i<<"\t"<<j-class_count<<endl;
                        LJ_12_6_interaction->addInteractionGroup(ecm_set[i], ecm_set[j-class_count]);
                        break;
                    case 2:
//                        cout<<"ecm-ecm: excluded_volume->(ecm_set[i], ecm_set[j]) = "<<i<<"\t"<<j-class_count<<endl;
                        excluded_volume->addInteractionGroup(ecm_set[i], ecm_set[j-class_count]);
                        break;
                    default:
                        break;
                }
            }
        }
    }
    
    
    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    std::vector<Vec3> initialPosInNm;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        //        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atoms[n].mass);
        // Convert the initial position to nm and append to the array.
        const Vec3 posInNm(atoms[n].initPosInAng[0] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[1] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[2] * OpenMM::NmPerAngstrom);
        initialPosInNm.push_back(posInNm);
        
        //add particles to the excluded volume force. The number of particles should be equal to the number particles in the system. The exluded interaction lists should be defined afterwards.
        
        excluded_volume->addParticle();
        LJ_12_6_interaction->addParticle();
        
        
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
                    
                    FENEs.push_back(new OpenMM::CustomBondForce("k_bond*lmin*lmin*(((lmin/1.5)/(r-(lmin/1.5)))^6)*step(le1-r)+(-0.5*k_bond*lmax*lmax*log(1-(r*r/lmax*lmax)))*step(r-le0);"));
                    FENEs[FENE_index]->addGlobalParameter("lmin",   bonds[i].FENE_lmin
                                                 * OpenMM::NmPerAngstrom);
                    FENEs[FENE_index]->addGlobalParameter("le0",   bonds[i].FENE_le0//);
                                                 * OpenMM::NmPerAngstrom);
                    FENEs[FENE_index]->addGlobalParameter("le1",   bonds[i].FENE_le1
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
                                                 * OpenMM::KJPerKcal
                                                 * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm *  OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                    system.addForce(X4harmonics[X4harmonic_index]);
                }
                X4harmonics[X4harmonic_index]->addBond(atom[0], atom[1]);
            }
                break;
        }
        
        
    }
    
    // Add the list of atom pairs that are excluded from the excluded volume force.
    excluded_volume->createExclusionsFromBonds(excluded_bonds, 1);
    
    LJ_12_6_interaction->createExclusionsFromBonds(excluded_bonds, 1);
    
    
    //DFs = DihedralForces
    set <std::string> DFs_classes;
    int DFs_index = -1;
    
//    cout<<"Here!!!!\n\n";
    
    for (int i=0; dihedrals[i].type != EndOfList; ++i) {
        
        auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
        if (DFs_item == DFs_classes.end()) {
            
            DFs_classes.insert(dihedrals[i].class_label);
            DFs_index++;
            
            DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "K_bend*(cos(dihedral(p1,p2,p3,p4)))"));
            DihedralForces[DFs_index]->addGlobalParameter("K_bend", dihedrals[i].bending_stiffness_value
                                                          * OpenMM::KJPerKcal
                                                          * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
          //  system.addForce(DihedralForces[DFs_index]);
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
    
    
    
    omm->integrator = new OpenMM::VerletIntegrator(stepSizeInFs * OpenMM::PsPerFs);
    omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform);
    omm->context->setPositions(initialPosInNm);
    platformName = omm->context->getPlatform().getName();
    
    return omm;
}

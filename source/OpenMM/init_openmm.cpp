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
                                 TimeDependantData*     time_dependant_data,
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
    MyOpenMMData*       omm = new MyOpenMMData();

    // Create a System and Force objects within the System.
    OpenMM::System&     system = *(omm->system = new OpenMM::System());

    // Retain a reference to each force object so we can fill in the forces.
    // Note: the System owns the force objects and will take care of deleting them;
    // don't do it yourself!

    cout<<"Defining interactions...\n";

    
    // Create a vector of handles for the force objects. These handles will be added to the system. Each handle in the list will be associated with a class instance.
    vector<OpenMM::CustomNonbondedForce*> ExcludedVolumes;
    vector<OpenMM::CustomNonbondedForce*> LJ_12_6_interactions;
    vector<OpenMM::CustomExternalForce*> ext_force;

//    std::vector< std::pair< int, int > > excluded_bonds;
    
    
    set_interactions(atoms,
                     bonds,
                     membrane_set,
                     actin_set,
                     ecm_set,
                     chromatin_set,
                     interaction_map,
                     ext_force,
                     LJ_12_6_interactions,
                     ExcludedVolumes,
                     system);
    
    
    std::vector<Vec3> initialPosInNm;
    std::vector<Vec3> initialVelInNmperPs;
    add_particles_to_system_and_forces(atoms,
                                       initialPosInNm,
                                       initialVelInNmperPs,
                                       LJ_12_6_interactions,
                                       ExcludedVolumes,
                                       system);
    
    omm->EV = ExcludedVolumes;
    // Create an array of harmonic spring force objects to add to the system.
    
    //for time-dependant external force
    time_dependant_data->ext_force = ext_force;
    OpenMM::HarmonicBondForce*      HarmonicBond = new OpenMM::HarmonicBondForce();
    OpenMM::HarmonicBondForce*      Kelvin_VoigtBond = new OpenMM::HarmonicBondForce();
    vector<OpenMM::CustomBondForce*>X4harmonics;
    vector<OpenMM::CustomBondForce*> FENEs;
    
    set_bonded_forces(bonds,
                      HarmonicBond,
                      Kelvin_VoigtBond,
                      X4harmonics,
                      FENEs,
                      time_dependant_data,
                      system);
    
    omm->harmonic = HarmonicBond;
    omm->x4harmonic=X4harmonics[0];
    time_dependant_data->Kelvin_VoigtBond = Kelvin_VoigtBond;
    time_dependant_data->Kelvin_Nominal_length_calc();

    // Add the list of atom pairs that are excluded from the excluded volume force.
    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".


    
    
    vector<OpenMM::CustomCompoundBondForce*> DihedralForces;
    set_dihedral_forces(dihedrals,
                        DihedralForces,
                        system);

    if (dihedrals[0].type != EndOfList) {
        omm->Dihedral = DihedralForces;
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

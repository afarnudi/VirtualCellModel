#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"

const int EndOfList=-1;

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
MyOpenMMData* myInitializeOpenMM(const MyAtomInfo    atoms[],
                                 double              stepSizeInFs,
                                 std::string&        platformName,
                                 Bonds*              bonds,
                                 Dihedrals*          dihedrals)
{
    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory
    (OpenMM::Platform::getDefaultPluginsDirectory());
    
    // Allocate space to hold OpenMM objects while we're using them.
    MyOpenMMData* omm = new MyOpenMMData();
    
    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the System owns
    // the force objects and will take care of deleting them; don't do it yourself!
    OpenMM::System&                 system      = *(omm->system = new OpenMM::System());
    OpenMM::HarmonicBondForce&      bondStretch = *new OpenMM::HarmonicBondForce();
    system.addForce(&bondStretch);
    
    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    std::vector<Vec3> initialPosInNm;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        //        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atoms[0].mass);
        // Convert the initial position to nm and append to the array.
        const Vec3 posInNm(atoms[n].initPosInAng[0] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[1] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[2] * OpenMM::NmPerAngstrom);
        initialPosInNm.push_back(posInNm);
    }
    
    
    
    /**
     * Here we use the OpenMM's "CustomBondForce" to implament the FENE spring.
     */
    OpenMM::CustomBondForce* force2 = new OpenMM::CustomBondForce("0.5*k*(r-r0)^2");
    for (int i=0; bonds[i].type != EndOfList; ++i) {
        const int*      atom = bonds[i].atoms;
//        const BondType& bond = bondType[bonds[i].type];
        
        
        // Note factor of 2 for stiffness below because Amber specifies the constant
        // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
        // it as used in the force term kx, with energy kx^2/2.
        bondStretch.addBond(atom[0], atom[1],
                            bonds[i].nominalLengthInAngstroms
                            * OpenMM::NmPerAngstrom,
                            bonds[i].stiffnessInKcalPerAngstrom2
                            * 10 * OpenMM::KJPerKcal
                            * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    }
    
    
    for (int i=0; dihedrals[i].type != EndOfList; ++i) {
                /**
                 * Here we use the OpenMM's "CustomCompoundBondForce" to difine the bending force between two neighbouring membrane trinagles.
                 */
                OpenMM::CustomCompoundBondForce* force = new OpenMM::CustomCompoundBondForce(4, "K_bend*(cos(dihedral(p1,p2,p3,p4)))");
                force->addGlobalParameter("K_bend",dihedrals[i].bending_stiffness_value);
        
                force->addBond(dihedrals[i].atoms);
                system.addForce(force);
        
        }
    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero.
    omm->integrator = new OpenMM::VerletIntegrator(stepSizeInFs * OpenMM::PsPerFs);
    omm->context    = new OpenMM::Context(*omm->system, *omm->integrator);
    omm->context->setPositions(initialPosInNm);
    
    platformName = omm->context->getPlatform().getName();
    return omm;
}

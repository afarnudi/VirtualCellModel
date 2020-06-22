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
                                 vector<vector<set<int> >  >    &chromatin_set,
                                 //                                 vector<set<int> >      &chromatin_set,
                                 vector<vector<int> >   interaction_map)
{
    const string cbp_plugin_location="/scratch/alifarnudi/local/openmm/lib/plugins";
    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());
    //OpenMM::Platform::loadPluginsFromDirectory(cbp_plugin_location);
    
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
    vector<OpenMM::CustomExternalForce*>  ext_force;
    
    
    OpenMM::CMMotionRemover* comremover;
    comremover = new OpenMM::CMMotionRemover(GenConst::CMMotionRemoverStep);
    
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
    
    omm->LJ = LJ_12_6_interactions;
    omm->EV = ExcludedVolumes;
    // Create an array of harmonic spring force objects to add to the system.
    
    
    //for time-dependant external force
    time_dependant_data->ext_force = ext_force;
    OpenMM::HarmonicBondForce*      HarmonicBond = new OpenMM::HarmonicBondForce();
    OpenMM::HarmonicBondForce*      Kelvin_VoigtBond = new OpenMM::HarmonicBondForce();
    vector<OpenMM::CustomBondForce*> X4harmonics;
    vector<OpenMM::CustomBondForce*> FENEs;
    //for calculating the force between nano_particles which is induced by membrane.
    //OpenMM::HarmonicBondForce* calcforce=new OpenMM::HarmonicBondForce(); 
    
    //for creating a non-spherical nano-particle
    //OpenMM::HarmonicBondForce* nonspherical=new OpenMM::HarmonicBondForce();
    
     //calcforce *****************************
    //calcforce->addBond(2574,2587, 6* OpenMM::NmPerAngstrom, 0* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    //system.addForce(calcforce);
     
     
    /* //non-spherical
     nonspherical->addBond(2574,2587, 2* OpenMM::NmPerAngstrom, 4000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
     nonspherical->addBond(2573,2586, 2* OpenMM::NmPerAngstrom, 4000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
     nonspherical->addBond(2572,2585, 2* OpenMM::NmPerAngstrom, 4000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
     nonspherical->addBond(2569,2582, 2* OpenMM::NmPerAngstrom, 4000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
     nonspherical->addBond(2567,2580, 2* OpenMM::NmPerAngstrom, 4000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
     nonspherical->addBond(2564,2577, 2* OpenMM::NmPerAngstrom, 4000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
     nonspherical->addBond(2562,2575, 2* OpenMM::NmPerAngstrom, 4000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    system.addForce(nonspherical);
    /*
    nonspherical->addBond(12,25, 2* OpenMM::NmPerAngstrom, 10000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    nonspherical->addBond(11,24, 2* OpenMM::NmPerAngstrom, 10000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    nonspherical->addBond(10,23, 2* OpenMM::NmPerAngstrom, 10000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    nonspherical->addBond(7,20, 2* OpenMM::NmPerAngstrom, 10000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    nonspherical->addBond(5,18, 2* OpenMM::NmPerAngstrom, 10000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm); 
    nonspherical->addBond(2,15, 2* OpenMM::NmPerAngstrom, 10000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    nonspherical->addBond(0,13, 2* OpenMM::NmPerAngstrom, 10000* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    system.addForce(nonspherical);
    */
    set_bonded_forces(bonds,
                      HarmonicBond,
                      Kelvin_VoigtBond,
                      X4harmonics,
                      FENEs,
                      time_dependant_data,
                      system);
    
    
    omm->harmonic = HarmonicBond;
    //omm->calcforce=calcforce;
    omm->x4harmonic=X4harmonics;
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
    
    //cout<<"platform default directory path = "<<OpenMM::Platform::getDefaultPluginsDirectory()<<endl;
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
    
    
    std::vector<std::map<std::string, std::string> > device_properties;
    if (platform.getName() == "OpenCL") {
        cout<<"Available devices on the "<<platform.getName()<<" platform:\n";
        int counter=0;
        for (int i=0; i<10; i++) {
            for (int j=0; j<10; j++) {
                try {
                    std::map<std::string, std::string> temp_device_properties;
                    temp_device_properties["OpenCLPlatformIndex"]=std::to_string(i);
                    temp_device_properties["OpenCLDeviceIndex"]=std::to_string(j);
                    OpenMM::System temp_system;
                    temp_system.addParticle(1.0);
                    OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
                    OpenMM::Context temp_context(temp_system, temp_inegrator, platform, temp_device_properties);
                    std::vector<std::string> platform_devices = platform.getPropertyNames();
                    cout<<counter<<" : ";
                    for (auto & name : platform_devices){
                        if (name == "DeviceIndex" || name == "OpenCLPlatformIndex") {
                            continue;
                        } else {
                            cout<<"\t"<<name<<"\t"<<platform.getPropertyValue(temp_context, name)<<endl;
                        }
                    }
                    cout<<"------------------------"<<endl;
                    counter++;
                    device_properties.push_back(temp_device_properties);
                } catch (const std::exception& e) {
                    
                }
            }
        }
    } else if (platform.getName() == "CUDA") {
        cout<<"Available devices on the "<<platform.getName()<<" platform:\n";
        int counter=0;
        for (int i=0; i<10; i++) {
            for (int j=0; j<10; j++) {
                try {
                    std::map<std::string, std::string> temp_device_properties;
                    temp_device_properties["CudaPlatformIndex"]=std::to_string(i);
                    temp_device_properties["CudaDeviceIndex"]=std::to_string(j);
                    OpenMM::System temp_system;
                    temp_system.addParticle(1.0);
                    OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
                    OpenMM::Context temp_context(temp_system, temp_inegrator, platform, temp_device_properties);
                    std::vector<std::string> platform_devices = platform.getPropertyNames();
                    cout<<counter<<" : ";
                    for (auto & name : platform_devices){
                        if (name == "DeviceIndex" || name == "CUDAPlatformIndex") {
                            continue;
                        } else {
                            cout<<"\t"<<name<<"\t"<<platform.getPropertyValue(temp_context, name)<<endl;
                        }
                    }
                    cout<<"------------------------"<<endl;
                    counter++;
                    device_properties.push_back(temp_device_properties);
                } catch (const std::exception& e) {
                    
                }
            }
        }
    } else if (platform.getName() == "CPU") {
        OpenMM::System temp_system;
        temp_system.addParticle(1.0);
        OpenMM::VerletIntegrator temp_inegrator(stepSizeInFs * OpenMM::PsPerFs);
        OpenMM::Context temp_context(temp_system, temp_inegrator, platform);
        std::vector<std::string> platform_devices = platform.getPropertyNames();
        cout<<"CPU properties:\n";
        for (auto & name : platform_devices){
            cout<<"\t"<<name<<"\t"<<platform.getPropertyValue(temp_context, name)<<endl;
        }
        cout<<endl;
    }
    
    int device_id=0;
    if (device_properties.size()>1) {
        cout<<"Please choose a device (index): \n";
        std::cin>>device_id;
    }
    
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
            
//            omm->integrator = new OpenMM::LangevinIntegrator(GenConst::temperature,
//                                                             GenConst::frictionInPs,
//                                                             stepSizeInFs * OpenMM::PsPerFs);
            
            omm->Lintegrator = new OpenMM::LangevinIntegrator(GenConst::temperature,
                                                             GenConst::frictionInPs,
                                                             stepSizeInFs * OpenMM::PsPerFs);
            break;
    }
    if (GenConst::CMMotionRemover) {
        omm->system->addForce(comremover);
    }
    
    for (int i=0; atoms[i].type != EndOfList; i++) {
        if (atoms[i].mass < 0.0001) {
            if (atoms[i].class_label == "Chromatin") {
                omm->system->setParticleMass(i, 0);
                
                OpenMM::TwoParticleAverageSite* vsite_pars;
                vsite_pars =  new OpenMM::TwoParticleAverageSite(atoms[i].vsite_atoms[0], atoms[i].vsite_atoms[1], atoms[i].Vsite_weights[0], atoms[i].Vsite_weights[0]);
                
                omm->system->setVirtualSite(i, vsite_pars);
            }
        }
    }
    
    
    
    if (platform.getName() == "CPU" || platform_id==0) {
        switch (GenConst::Integrator_type) {
            case 0:
                omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform);
                break;
                
            case 1:
                omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform);
                break;
            case 2:
                
                omm->context    = new OpenMM::Context(*omm->system, *omm->Lintegrator, platform);
                break;
        }
        
    } else {
        switch (GenConst::Integrator_type) {
            case 0:
                omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform, device_properties[device_id]);
                break;
                
            case 1:
                omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform, device_properties[device_id]);
                break;
            case 2:
                
                omm->context    = new OpenMM::Context(*omm->system, *omm->Lintegrator, platform, device_properties[device_id]);
                break;
        }
        
    }
        
    omm->context->setPositions(initialPosInNm);
    omm->context->setVelocities(initialVelInNmperPs);
    
    platformName = omm->context->getPlatform().getName();
    
    return omm;
}

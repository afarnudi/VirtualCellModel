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
                                 vector<vector<set<int> > > &chromatin_set,
                                 ArgStruct_VCM           userinputs,
                                 NonBondInteractionMap  &interaction_map)
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
    
    cout<<TOMM<<"Defining interactions..."<<TRESET<<endl;;
    
    
    // Create a vector of handles for the force objects. These handles will be added to the system. Each handle in the list will be associated with a class instance.
    vector<OpenMM::CustomNonbondedForce*> ExcludedVolumes;
    vector<OpenMM::CustomNonbondedForce*> LJ_12_6_interactions;
    vector<OpenMM::CustomExternalForce*>  ext_force;
    
    
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
    vector<OpenMM::CustomBondForce*> Contractiles;
    vector<OpenMM::CustomBondForce*> HillBonds;
    vector<OpenMM::CustomBondForce*> Harmonic_minmax;
    vector<OpenMM::CustomBondForce*> KFs;
    //OpenMM::HarmonicAngleForce*     HarmonicAngle = new OpenMM::HarmonicAngleForce();
    
    
    
    
    set_bonded_forces(bonds,
                      HarmonicBond,
                      Kelvin_VoigtBond,
                      X4harmonics,
                      FENEs,
                      Contractiles,
                      KFs,
                      HillBonds,
                      Harmonic_minmax,
                      time_dependant_data,
                      system);
    
    
    omm->harmonic = HarmonicBond;
    //omm->calcforce=calcforce;
    omm->x4harmonic=X4harmonics;
    time_dependant_data->Kelvin_VoigtBond = Kelvin_VoigtBond;
    time_dependant_data->Hill_force = HillBonds;
    time_dependant_data->k_force = KFs;
    time_dependant_data->Kelvin_Nominal_length_calc();
    
    // Add the list of atom pairs that are excluded from the excluded volume force.
    // the second input is an integer, bondCutoff; OpenMM defines bondCutoff as "pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions".
    
    
    //pbcgel begin
//    cout<<"GenConst::Lbox "<<GenConst::Lbox<<endl;
    vector<vector<int> > pbcbonds;
    vector<int> pair;
    pair.resize(2,0);
    pair[0]=42; pair[1]=48;
    pbcbonds.push_back(pair);
    pair[0]=43; pair[1]=49;
    pbcbonds.push_back(pair);
    pair[0]=44; pair[1]=46;
    pbcbonds.push_back(pair);
    pair[0]=45; pair[1]=47;
    pbcbonds.push_back(pair);
    pair[0]=50; pair[1]=51;
    pbcbonds.push_back(pair);
    pair[0]=52; pair[1]=53;
    pbcbonds.push_back(pair);
    pair[0]=54; pair[1]=55;
    pbcbonds.push_back(pair);
    OpenMM::CustomBondForce* pbcharmonic = new OpenMM::CustomBondForce("k_pbc*r^2");
    pbcharmonic->addGlobalParameter("k_pbc",10000);
    system.addForce(pbcharmonic);
    
    for (int pbcind=0; pbcind<pbcbonds.size(); pbcind++) {
        pbcharmonic->addBond(pbcbonds[pbcind][0]+1002, pbcbonds[pbcind][1]+1002);
        if (GenConst::Periodic_box) {
            pbcharmonic->setUsesPeriodicBoundaryConditions(true);
        }
//        omm->harmonic->addBond(pbcbonds[pbcind][0]+1002, pbcbonds[pbcind][1]+1002,
//                              0,
//                              1000);
    }
//
//    if (GenConst::Periodic_box) {
//        HarmonicBond->setUsesPeriodicBoundaryConditions(true);
//    }
    //pbcgel end
    
    
    vector<OpenMM::CustomCompoundBondForce*> DihedralForces;
    set_dihedral_forces(dihedrals,
                        DihedralForces,
                        system);
    
    if (dihedrals[0].type != EndOfList) {
        omm->Dihedral = DihedralForces;
    }
    
//    GenConst::Periodic_box=true;
//    GenConst::Lbox=200;
    std::vector<Vec3> pbcxyz;
    
    if (GenConst::Periodic_box) {
        pbcxyz.resize(3);
        
        pbcxyz[0][0]=GenConst::Lbox;
        pbcxyz[0][1]=0;
        pbcxyz[0][2]=0;
        
        pbcxyz[1][0]=0;
        pbcxyz[1][1]=GenConst::Lbox;
        pbcxyz[1][2]=0;
        
        pbcxyz[2][0]=0;
        pbcxyz[2][1]=0;
        pbcxyz[2][2]=GenConst::Lbox;
        
        system.setDefaultPeriodicBoxVectors(pbcxyz[0], pbcxyz[1], pbcxyz[2]);
        
        pbcxyz[0][0]=0;
        pbcxyz[0][1]=0;
        pbcxyz[0][2]=0;
        
        pbcxyz[1][0]=0;
        pbcxyz[1][1]=0;
        pbcxyz[1][2]=0;
        
        pbcxyz[2][0]=0;
        pbcxyz[2][1]=0;
        pbcxyz[2][2]=0;
        
        system.getDefaultPeriodicBoxVectors(pbcxyz[0], pbcxyz[1], pbcxyz[2]);
        cout<<"Periodic Boundry Condition "<<TON<<"On"<<TRESET<<endl;
        cout<<"Periodic vectors (X,Y,Z):\n";
        cout<<TGRAY;
        for (int i=0; i<3; i++) {
            cout<<i<<": "<<pbcxyz[i][0]<<"\t"<<pbcxyz[i][1]<<"\t"<<pbcxyz[i][2]<<"\n";
        }
        cout<<TRESET;
    } else {
        cout<<"Periodic Boundry Condition "<<TOFF<<"Off"<<TRESET<<endl;
    }
    
    
    
    PlatformInfo platforminfo;
    if (userinputs.platforminput) {
        platforminfo = userinputs.platforminfo;
        get_platform_info(platforminfo);
    } else {
        platforminfo = get_platform_info();
    }
    
    OpenMM::Platform& platform = OpenMM::Platform::getPlatform(platforminfo.platform_id);
    generateHardwareReport(platforminfo);

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
        OpenMM::CMMotionRemover* comremover;
        comremover = new OpenMM::CMMotionRemover(GenConst::CMMotionRemoverStep);
        omm->system->addForce(comremover);
    }
    
    if (GenConst::MCBarostatFrequency != 0) {
        if (GenConst::MCBarostatTemperature<0) {
            GenConst::MCBarostatTemperature = GenConst::temperature;
        }
        OpenMM::MonteCarloBarostat* MCBarostat = new OpenMM::MonteCarloBarostat(GenConst::MCBarostatPressure, GenConst::MCBarostatTemperature, GenConst::MCBarostatFrequency);
        omm->system->addForce(MCBarostat);
    }
    
    for (int i=0; atoms[i].type != EndOfList; i++) {
        if (atoms[i].mass < 0.0001) {
            if (atoms[i].class_label == "Chromatin") {
                omm->system->setParticleMass(i, 0);
                
                OpenMM::TwoParticleAverageSite* vsite_pars;
                vsite_pars =  new OpenMM::TwoParticleAverageSite(atoms[i].vsite_atoms[0], atoms[i].vsite_atoms[1], atoms[i].Vsite_weights[0], atoms[i].Vsite_weights[0]);//I think the last one should be atoms[i].Vsite_weights[1]
                
                omm->system->setVirtualSite(i, vsite_pars);
            }
        }
    }
    
    
    
    if (platform.getName() == "CPU" || platforminfo.platform_id==0) {
        if ( omm->integrator != NULL ) {
            omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform);
        } else {
            omm->context    = new OpenMM::Context(*omm->system, *omm->Lintegrator, platform);
        }
        
        
        
    } else {
        if ( omm->integrator != NULL ) {
            omm->context    = new OpenMM::Context(*omm->system, *omm->integrator, platform, platforminfo.device_properties[platforminfo.platform_device_id]);
        } else {
            omm->context    = new OpenMM::Context(*omm->system, *omm->Lintegrator, platform, platforminfo.device_properties[platforminfo.platform_device_id]);
        }
        
        
    }
    
    omm->context->setPositions(initialPosInNm);
    omm->context->setVelocities(initialVelInNmperPs);
    
    
    platformName = omm->context->getPlatform().getName();
    
    const std::map <std::string, double> params = omm->context->getParameters();
    
    
    
    cout<<TGRAY<<params.size()<<endl;
    for(auto elem : params)
    {
        cout << elem.first << " " << elem.second << "\n";
    }
    cout<<"\n"<<TRESET;
    
    return omm;
}

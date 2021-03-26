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
                                 Angles*                angles,
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
    
    vector<string> loaderror = OpenMM::Platform::getPluginLoadFailures();
    for (auto &line: loaderror){
        cout<<line<<endl;
    }
    
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
    vector<OpenMM::CustomNonbondedForce*> WCAs;
    vector<OpenMM::CustomNonbondedForce*> WCAFCs;
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
                     WCAs,
                     WCAFCs,
                     generalParameters.MinimumForceDecleration,
                     system);

    
    
    std::vector<Vec3> initialPosInNm;
    std::vector<Vec3> initialVelInNmperPs;
    add_particles_to_system_and_forces(atoms,
                                       initialPosInNm,
                                       initialVelInNmperPs,
                                       LJ_12_6_interactions,
                                       ExcludedVolumes,
                                       WCAs,
                                       WCAFCs,
                                       interaction_map,
                                       system);
    
    omm->LJ  = LJ_12_6_interactions;
    omm->EV  = ExcludedVolumes;
    omm->WCA = WCAs;
    omm->WCAFC = WCAFCs;
    

    creatBondExclusion(bonds,
                       interaction_map,
                       LJ_12_6_interactions,
                       ExcludedVolumes,
                       WCAs,
                       WCAFCs);

    
    
    
    // Create an array of harmonic spring force objects to add to the system.
    
    
    //for time-dependant external force
    time_dependant_data->ext_force = ext_force;
    OpenMM::HarmonicBondForce*      HarmonicBond = new OpenMM::HarmonicBondForce();
    OpenMM::HarmonicBondForce*      Kelvin_VoigtBond = new OpenMM::HarmonicBondForce();
    vector<OpenMM::CustomBondForce*> X4harmonics;
    vector<OpenMM::CustomBondForce*> KremerGrests;
    vector<OpenMM::CustomBondForce*> Gompperbond;
    vector<OpenMM::CustomBondForce*> Gompperrep;
    vector<OpenMM::CustomBondForce*> Contractiles;
    vector<OpenMM::CustomBondForce*> HillBonds;
    vector<OpenMM::CustomBondForce*> Harmonic_minmax;
    vector<OpenMM::CustomBondForce*> KFs;
    //OpenMM::HarmonicAngleForce*     HarmonicAngle = new OpenMM::HarmonicAngleForce();
    
    if (!generalParameters.MinimumForceDecleration) {
        set_bonded_forces(bonds,
                          HarmonicBond,
                          Kelvin_VoigtBond,
                          X4harmonics,
                          KremerGrests,
                          Gompperbond,
                          Gompperrep,
                          Contractiles,
                          KFs,
                          HillBonds,
                          Harmonic_minmax,
                          time_dependant_data,
                          system);
    } else {
        set_bonded_forces(bonds,
//                          HarmonicBond,
//                          Kelvin_VoigtBond,
//                          X4harmonics,
                          KremerGrests,
//                          Gompperbond,
//                          Gompperrep,
//                          Contractiles,
//                          KFs,
//                          HillBonds,
//                          Harmonic_minmax,
//                          time_dependant_data,
                          system);
    }
    
    
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
    //    vector<vector<int> > pbcbonds;
    //    vector<int> pair;
    //    pair.resize(2,0);
    //    pair[0]=42; pair[1]=48;
    //    pbcbonds.push_back(pair);
    //    pair[0]=43; pair[1]=49;
    //    pbcbonds.push_back(pair);
    //    pair[0]=44; pair[1]=46;
    //    pbcbonds.push_back(pair);
    //    pair[0]=45; pair[1]=47;
    //    pbcbonds.push_back(pair);
    //    pair[0]=50; pair[1]=51;
    //    pbcbonds.push_back(pair);
    //    pair[0]=52; pair[1]=53;
    //    pbcbonds.push_back(pair);
    //    pair[0]=54; pair[1]=55;
    //    pbcbonds.push_back(pair);
    //    OpenMM::CustomBondForce* pbcharmonic = new OpenMM::CustomBondForce("k_pbc*r^2");
    //    pbcharmonic->addGlobalParameter("k_pbc",10000);
    //    system.addForce(pbcharmonic);
    //    int memnodes;
    //    if (membrane_set.size()!=0) {
    //        memnodes= membrane_set[0].size();
    //    } else {
    //        memnodes=0;
    //    }
    //    for (int pbcind=0; pbcind<pbcbonds.size(); pbcind++) {
    //        pbcharmonic->addBond(pbcbonds[pbcind][0]+memnodes, pbcbonds[pbcind][1]+memnodes);
    //        if (GenConst::Periodic_box) {
    //            pbcharmonic->setUsesPeriodicBoundaryConditions(true);
    //        }
    //    }
    //    MonteCarloAnisotropicBarostat(const Vec3 &defaultPressure, double defaultTemperature, bool scaleX = true, bool scaleY = true, bool scaleZ = true, int frequency = 25)
    //    bool anisotropicbarostat=true;
    //    if (anisotropicbarostat) {
    //        double basepressure = 0.0005;
    //        const Vec3 anisotropicpressure(basepressure,2*basepressure,2*basepressure);
    //        OpenMM::MonteCarloAnisotropicBarostat* AnisoMCBarostat = new OpenMM::MonteCarloAnisotropicBarostat(anisotropicpressure, 600,true,false,false,25);
    //        omm->system->addForce(AnisoMCBarostat);
    //    }
    
    //pbcgel end
    
    
    vector<OpenMM::CustomCompoundBondForce*> DihedralForces;
    set_dihedral_forces(dihedrals,
                        DihedralForces,
                        system);
    
    if (dihedrals[0].type != EndOfList) {
        omm->Dihedral = DihedralForces;
    }
    
    vector<OpenMM::CustomAngleForce*> AngleForces;
    if (!generalParameters.MinimumForceDecleration) {
        set_angle_forces(angles,
                         AngleForces,
                         system);
    } else {
        set_angle_forces_minimum(angles,
                                 AngleForces,
                                 system);
    }
    
    
    if (angles[0].type != EndOfList) {
        omm->Angle = AngleForces;
    }
    
    set_pbcvectors(system);
    
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
    
    
    
    if (generalParameters.Integrator_type=="Verlet") {
        omm->VerletIntegrator = new OpenMM::VerletIntegrator(stepSizeInFs * OpenMM::PsPerFs);
    } else if (generalParameters.Integrator_type=="Brownian"){
        omm->BrownianIntegrator = new OpenMM::BrownianIntegrator(generalParameters.temperature,
                                                                 generalParameters.frictionInPs,
                                                                 stepSizeInFs * OpenMM::PsPerFs);
    } else if (generalParameters.Integrator_type=="Langevin"){
        omm->LangevinIntegrator = new OpenMM::LangevinIntegrator(generalParameters.temperature,
                                                                 generalParameters.frictionInPs,
                                                                 stepSizeInFs * OpenMM::PsPerFs);
    } else if (generalParameters.Integrator_type=="Custom"){
//        set_customLangevin(omm, stepSizeInFs);
        set_multithermos(omm, interaction_map, stepSizeInFs, membrane_set, atoms);
    } else if (generalParameters.Integrator_type=="LangevinMinimise"){
        set_customLangevinforminimisation(omm, stepSizeInFs, generalParameters.MinimisationIntegraterRestriction);
    }
    
    
    
    if (generalParameters.CMMotionRemover) {
        OpenMM::CMMotionRemover* comremover;
        comremover = new OpenMM::CMMotionRemover(generalParameters.CMMotionRemoverStep);
        omm->system->addForce(comremover);
    }
    
    if (generalParameters.MCBarostatFrequency != 0) {
        if (generalParameters.MCBarostatTemperature<0) {
            generalParameters.MCBarostatTemperature = generalParameters.temperature;
        }
        OpenMM::MonteCarloBarostat* MCBarostat = new OpenMM::MonteCarloBarostat(generalParameters.MCBarostatPressure, generalParameters.MCBarostatTemperature, generalParameters.MCBarostatFrequency);
        omm->system->addForce(MCBarostat);
    }
    if (generalParameters.MCAnisoBarostatOn) {
        
        const Vec3 anisotropicpressure(generalParameters.MCAnisoBarostatPressure[0],generalParameters.MCAnisoBarostatPressure[1],generalParameters.MCAnisoBarostatPressure[2]);
        OpenMM::MonteCarloAnisotropicBarostat* AnisoMCBarostat = new OpenMM::MonteCarloAnisotropicBarostat(anisotropicpressure, generalParameters.MCAnisoBarostatTemperature, generalParameters.MCAnisoBarostatScaleXYZ[0], generalParameters.MCAnisoBarostatScaleXYZ[1], generalParameters.MCAnisoBarostatScaleXYZ[2],generalParameters.MCAnisoBarostatFrequency);
        omm->system->addForce(AnisoMCBarostat);
    }
    
//    ofstream output("system_fail.xml");
//    OpenMM::XmlSerializer::serialize<OpenMM::System>(omm->system, "System", output);
//    output.close();exit(0);
    
    if (platform.getName() == "CPU" || platforminfo.platform_id==0) {
        if (generalParameters.Integrator_type=="Verlet") {
            omm->context    = new OpenMM::Context(*omm->system, *omm->VerletIntegrator, platform);
        } else if (generalParameters.Integrator_type=="Brownian"){
            omm->context    = new OpenMM::Context(*omm->system, *omm->BrownianIntegrator, platform);
        } else if (generalParameters.Integrator_type=="Langevin"){
            omm->context    = new OpenMM::Context(*omm->system, *omm->LangevinIntegrator, platform);
        } else if (generalParameters.Integrator_type=="Custom"){
            omm->context    = new OpenMM::Context(*omm->system, *omm->CustomIntegrator, platform);
        } else if (generalParameters.Integrator_type=="LangevinMinimise"){
            omm->context    = new OpenMM::Context(*omm->system, *omm->LangevinMinimisation, platform);
        }
    } else {
        if ( generalParameters.Integrator_type=="Verlet" ) {
            omm->context    = new OpenMM::Context(*omm->system, *omm->VerletIntegrator, platform, platforminfo.device_properties[platforminfo.platform_device_id]);
        } else if (generalParameters.Integrator_type=="Brownian"){
            omm->context    = new OpenMM::Context(*omm->system, *omm->BrownianIntegrator, platform, platforminfo.device_properties[platforminfo.platform_device_id]);
        } else if (generalParameters.Integrator_type=="Langevin") {
            omm->context    = new OpenMM::Context(*omm->system, *omm->LangevinIntegrator, platform, platforminfo.device_properties[platforminfo.platform_device_id]);
        } else if (generalParameters.Integrator_type=="Custom"){
            omm->context    = new OpenMM::Context(*omm->system, *omm->CustomIntegrator, platform, platforminfo.device_properties[platforminfo.platform_device_id]);
        } else if (generalParameters.Integrator_type=="LangevinMinimise"){
            omm->context    = new OpenMM::Context(*omm->system, *omm->LangevinMinimisation, platform, platforminfo.device_properties[platforminfo.platform_device_id]);
        }
        
        
    }
    
    omm->context->setPositions(initialPosInNm);
    omm->context->setVelocities(initialVelInNmperPs);
    
    
    platformName = omm->context->getPlatform().getName();
    
    const std::map <std::string, double> params = omm->context->getParameters();
    
    cout<<flush;
    cout<<TGRAY<<params.size()<<endl;
    for(auto elem : params)
    {
        cout << elem.first << " " << elem.second << "\n";
    }
    cout<<"\n"<<TRESET;
    cout<<flush;
    
    
    
    
    
    
    
    
    
    return omm;
}

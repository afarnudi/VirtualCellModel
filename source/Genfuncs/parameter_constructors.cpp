#include <vector>
#include "maps.hpp"
#include "General_constants.h"

using namespace std;

struct GeneralParameters{
//    map<string, double> GP_double;
//    map<string, int>    GP_int;
//    map<string, bool>    GP_bool;
//    map<string, string> GP_string;
//    map<string, string> GP_help;
//    map<string, vector<double> > GP_vecdouble;
//
//    vector<double> Default_Periodic_box_vector0;
//    vector<double> Default_Periodic_box_vector1;
//    vector<double> Default_Periodic_box_vector2;
    map<string, vector<string> > GenParams;
    
    vector<string> values;
    GeneralParameters(){
        values.resize(2);
        
        values[0] ="10";
        values[1] ="#Simulation time length masured in pico seconds. Default value 10.";
        GenParams["SimulationTimeInPs"] = values;
        
        values[0] ="10";
        values[1] ="#Integration step size measured in femto seconds.  Default value 10.";
        GenParams["StepSizeInFs"] = values;
        
        values[0] ="1000";
        values[1] ="#Trajectory saving, data collection, etc wil be performed at the ReportIntervalInFs time intervales measured in femto seconds. Default value 1000.";
        GenParams["ReportIntervalInFs"] = values;
        
        values[0] ="1000";
        values[1] ="#Simulation enviroment size (cube). Default value 1000.";
        GenParams["SimulationBoxLength"] = values;
        
        values[0] ="5";
        values[1] ="#The friction coefficient which couples the system to the heat bath (in inverse pico seconds). Required by the Brownian and Langevin integrators. Default value 5.";
        GenParams["FrictionInPs"] = values;
        
        values[0] ="310";
        values[1] ="#The thermostat temperature (in Kelvin). Required by the Brownian and Langevin integrators. Default value 310";
        GenParams["Temperature"] = values;
        
        values[0] ="0";
        values[1] ="#The pressure acting on the system (in bar) through OpenMM's MonteCarloBarostat. Default 0";
        GenParams["MCBarostatPressure"] = values;
        
        values[0] ="0";
        values[1] ="#The temperature at which OpenMM's MonteCarloBarostat will think the system is being maintained (in Kelvin). Default thermostat's Temperature";
        GenParams["MCBarostatTemperature"] = values;
        
        values[0] ="0";
        values[1] ="#The frequency at which the MonteCarloBarostat pressure changes should be attempted (in time steps). If zero MCBarostat will be disabled. Default 0";
        GenParams["MCBarostatFrequency"] = values;
        
        values[0] ="100";
        values[1] ="#This option is not available. Note to the developer: This should be moved to the Membrane. Default value 100.";
        GenParams["MCStep"] = values;
        
        values[0] ="100";
        values[1] ="#This option is not available. Note to the developer: This should be moved to the Membrane. Default value 100.";
        GenParams["MemFluidity"] = values;
        
        values[0] ="10";
        values[1] ="#The number of steps that the centre of mass velocity is set to zero using OpenMM's CMMotionRemover. Default 10";
        GenParams["CMMotionRemoverStep"] = values;
        
        values[0] ="V";
        values[1] ="#Set the integrator type. 'Temperature' and 'FrictionCoeff' need to be set for Langevin and Brownian. V: Verlet. B: Brownian. L: Langevin. Default value 'V'.";
        GenParams["Integrator"] = values;
        
        values[0] ="true";
        values[1] ="#Collect energy parameters for the potentials (expensive) during each Report_Interval_In_Fs time point. Default true";
        GenParams["ReportEnergy"] = values;
        
        values[0] ="false";
        values[1] ="#Writes velocities and forces (cheap) of particles during each ReportIntervalInFs time point to the disk. Default false";
        GenParams["WriteVelocitiesForces"] = values;
        
        values[0] ="false";
        values[1] ="#Writes velocities and forces (cheap) of particles during each ReportIntervalInFs time point to the disk. Default false";
        GenParams["WriteVelocitiesForces"] = values;
        
        values[0] ="1000 0 0";
        values[1] ="#Periodic box vector (1 of 3). Default value (1000, 0, 0).";
        GenParams["PeriodicBoxVector0"] = values;
        
        values[0] ="0 1000 0";
        values[1] ="#Periodic box vector (2 of 3). Default value (0, 1000, 0).";
        GenParams["PeriodicBoxVector1"] = values;
        
        values[0] ="0 0 1000";
        values[1] ="#Periodic box vector (3 of 3). Default value (0, 0, 1000).";
        GenParams["PeriodicBoxVector2"] = values;
    }
    
//    GeneralParameters(){
//        GP_double["SimulationTimeInPs"] = 10;
//        GP_help  ["SimulationTimeInPs"] = "#Simulation time length masured in pico seconds. Default value 10.";
//
//        GP_double["StepSizeInFs"] = 10;
//        GP_help  ["StepSizeInFs"] = "#Integration step size measured in femto seconds.  Default value 10.";
//
//        GP_double["ReportIntervalInFs"] = 1000;
//        GP_help  ["ReportIntervalInFs"] = "#Trajectory saving, data collection, etc wil be performed at the ReportIntervalInFs time intervales measured in femto seconds. Default value 1000.";
//
//        GP_double["SimulationBoxLength"] = 1000;
//        GP_help  ["SimulationBoxLength"] = "#Simulation enviroment size (cube). Default value 1000.";
//
//        GP_double["FrictionInPs"] = 5;
//        GP_help  ["FrictionInPs"] = "#The friction coefficient which couples the system to the heat bath (in inverse pico seconds). Required by the Brownian and Langevin integrators. Default value 5.";
//
//        GP_double["Temperature"] = 310;
//        GP_help  ["Temperature"] = "#The thermostat temperature (in Kelvin). Required by the Brownian and Langevin integrators. Default value 310";
//
//        GP_double["MCBarostatPressure"] = 0;
//        GP_help  ["MCBarostatPressure"] = "#The pressure acting on the system (in bar) through OpenMM's MonteCarloBarostat. Default 0";
//
//        GP_double["MCBarostatTemperature"] = 0;
//        GP_help  ["MCBarostatTemperature"] = "#The temperature at which OpenMM's MonteCarloBarostat will think the system is being maintained (in Kelvin). Default thermostat's Temperature";
//
//        GP_double["MCBarostatFrequency"] = 0;
//        GP_help  ["MCBarostatFrequency"] = "#The frequency at which the MonteCarloBarostat pressure changes should be attempted (in time steps). If zero MCBarostat will be disabled. Default 0";
//
//
//
//
//        GP_int   ["MCStep"] = 100;
//        GP_help  ["MCStep"] = "#This option is not available. Note to the developer: This should be moved to the Membrane. Default value 100.";
//
//        GP_int   ["MemFluidity"] = 100;
//        GP_help  ["MemFluidity"] = "#This option is not available. Note to the developer: This should be moved to the Membrane. Default value 100.";
//
//        GP_int   ["CMMotionRemoverStep"] = 10;
//        GP_help  ["CMMotionRemoverStep"] = "#The number of steps that the centre of mass velocity is set to zero using OpenMM's CMMotionRemover. Default 10";
//
//
//
//
//
//        GP_string["Integrator"] = "V";
//        GP_help  ["Integrator"] = "#Set the integrator type. 'Temperature' and 'FrictionCoeff' need to be set for Langevin and Brownian. V: Verlet. B: Brownian. L: Langevin. Default value 'V'.";
//
//
//
//        GP_bool  ["ReportEnergy"] = true;
//        GP_help  ["ReportEnergy"] = "#Collect energy parameters for the potentials (expensive) during each Report_Interval_In_Fs time point. Default true";
//
//        GP_bool  ["WriteVelocitiesForces"] = false;
//        GP_help  ["WriteVelocitiesForces"] = "#Writes velocities and forces (cheap) of particles during each ReportIntervalInFs time point to the disk. Default false";
//
//        GP_bool  ["WriteBondsToPDB"] = false;
//        GP_help  ["WriteVelocitiesForces"] = "#Writes velocities and forces (cheap) of particles during each ReportIntervalInFs time point to the disk. Default false";
//
//
//
//
//
//        Default_Periodic_box_vector0.resize(3,0);
//        Default_Periodic_box_vector0[0]= 1000;
//        GP_vecdouble["PeriodicBoxVector0"] = Default_Periodic_box_vector0;
//        GP_help     ["PeriodicBoxVector0"] = "#Periodic box vector (1 of 3). Default value (1000, 0, 0).";
//
//        Default_Periodic_box_vector1.resize(3,0);
//        Default_Periodic_box_vector0[1]= 1000;
//        GP_vecdouble["PeriodicBoxVector1"] = Default_Periodic_box_vector1;
//        GP_help     ["PeriodicBoxVector1"] = "#Periodic box vector (2 of 3). Default value (0, 1000, 0).";
//
//        Default_Periodic_box_vector2.resize(3,0);
//        Default_Periodic_box_vector0[2]= 1000;
//        GP_vecdouble["PeriodicBoxVector2"] = Default_Periodic_box_vector2;
//        GP_help     ["PeriodicBoxVector2"] = "#Periodic box vector (3 of 3). Default value (0, 1000, 0).";
//    }
};

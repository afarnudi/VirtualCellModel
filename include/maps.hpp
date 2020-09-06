//
//  General_constants.h
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 26/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef maps_hpp
#define maps_hpp
#include <map>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <unordered_map>

using namespace std;

void read_general_parameters(string input_file_name, vector<string> &membrane_config_list, vector<string> &chromatin_config_list, vector<string> &actin_config_list, vector<string> &ecm_config_list, vector<string> &pointparticle_config_list);

void set_parameter(map<string, double> &general_param_map, string param_name, double param_value);

/**The interaction map specifies the class instances that are allowed to interact with one another and the nature of their interaction. The map is a text file that lists the class member instances of the enviroment followed by the interaction specifier that is represented with an integer. example:
 * For a simple cell that contains an outer layer of membrane, an actin network, a nucleus membrane, 1 chromatin, and an ECM substrait the interaction will be:
 * mem0 1
 * mem1 0   1
 * act0 1   1   1
 * ecm0 1   0   0   1
 * chr0 0   1   0   0   1
 *
 * Here we assume that mem0 is the outer membrane. The class instance indecies are set during runtime in the order in which the respective configuration file directory is written in the general configuration file. The list of labels on the left (mem0, mem1, etc) is ignored by the programme and its function is for the user to keep track of the columens. It should be noted that the programme expects to come across a word in each line (which is ignored) so the user must not delete the labels altogether. But the actual label written in the interaction map is up to the user as long as it is declared in a single word, for example 'abcdef1234ghi and not 'abc 12 def'.
 * The programme sets the interaction between the class instances in the following order: Membranes, Actins, ECMs, Chromatins, Point Particles. The order in which these interactions are specified in the map is important.
 */
void read_interaction_map(vector<vector<int> > &inter_map);

void configfile_generator(int status);

struct GeneralParameters{
    unordered_map<string, vector<string> > GenParams;
    
    vector<string> values;
    GeneralParameters(){
        values.resize(2);
        
        values[0] ="0 0 1000";
        values[1] ="#Periodic box vector (3 of 3). Default value (0, 0, 1000).";
        GenParams["PeriodicBoxVector2"] = values;
        
        values[0] ="0 1000 0";
        values[1] ="#Periodic box vector (2 of 3). Default value (0, 1000, 0).";
        GenParams["PeriodicBoxVector1"] = values;
        
        values[0] ="1000 0 0";
        values[1] ="#Periodic box vector (1 of 3). Default value (1000, 0, 0).";
        GenParams["PeriodicBoxVector0"] = values;
        
        values[0] ="false";
        values[1] ="#Writes velocities and forces (cheap) of particles during each ReportIntervalInFs time point to the disk. Default false";
        GenParams["WriteVelocitiesForces"] = values;
        
        values[0] ="true";
        values[1] ="#Collect energy parameters for the potentials (expensive) during each Report_Interval_In_Fs time point. Default true";
        GenParams["ReportEnergy"] = values;
        
        values[0] ="10";
        values[1] ="#The number of steps that the centre of mass velocity is set to zero using OpenMM's CMMotionRemover. Default 10";
        GenParams["CMMotionRemoverStep"] = values;
        
        values[0] ="100";
        values[1] ="#This option is not available. Note to the developer: This should be moved to the Membrane. Default value 100.";
        GenParams["MemFluidity"] = values;
        
        values[0] ="100";
        values[1] ="#This option is not available. Note to the developer: This should be moved to the Membrane. Default value 100.";
        GenParams["MCStep"] = values;
        
        values[0] ="0";
        values[1] ="#The frequency at which the MonteCarloBarostat pressure changes should be attempted (in time steps). If zero MCBarostat will be disabled. Default 0";
        GenParams["MCBarostatFrequency"] = values;
        
        values[0] ="0";
        values[1] ="#The temperature at which OpenMM's MonteCarloBarostat will think the system is being maintained (in Kelvin). Default thermostat's Temperature";
        GenParams["MCBarostatTemperature"] = values;
        
        values[0] ="0";
        values[1] ="#The pressure acting on the system (in bar) through OpenMM's MonteCarloBarostat. Default 0";
        GenParams["MCBarostatPressure"] = values;
        
        values[0] ="310";
        values[1] ="#The thermostat temperature (in Kelvin). Required by the Brownian and Langevin integrators. Default value 310";
        GenParams["Temperature"] = values;
        
        values[0] ="5";
        values[1] ="#The friction coefficient which couples the system to the heat bath (in inverse pico seconds). Required by the Brownian and Langevin integrators. Default value 5.";
        GenParams["FrictionInPs"] = values;
        
        values[0] ="V";
        values[1] ="#Set the integrator type. 'Temperature' and 'FrictionCoeff' need to be set for Langevin and Brownian. V: Verlet. B: Brownian. L: Langevin. Default value 'V'.";
        GenParams["Integrator"] = values;
        
        values[0] ="1000";
        values[1] ="#Simulation enviroment size (cube). Default value 1000.";
        GenParams["SimulationBoxLength"] = values;
        
        values[0] ="1000";
        values[1] ="#Trajectory saving, data collection, etc wil be performed at the ReportIntervalInFs time intervales measured in femto seconds. Default value 1000.";
        GenParams["ReportIntervalInFs"] = values;
        
        values[0] ="10";
        values[1] ="#Integration step size measured in femto seconds.  Default value 10.";
        GenParams["StepSizeInFs"] = values;
        
        values[0] ="10";
        values[1] ="#Simulation time length masured in pico seconds. Default value 10.";
        GenParams["SimulationTimeInPs"] = values;
    }
    
    unordered_map<string, vector<string> > get_unorderedmap(){
        return GenParams;
    }
};


#endif /* maps_hpp */

//
//  OpenMM_structs.h
//  Membrae
//
//  Created by Ali Farnudi on 10/06/2019.
//  Copyright © 2019 Ali Farnudi. All rights reserved.
//

#ifndef ConfigfileStructs_hpp
#define ConfigfileStructs_hpp

#include <string>
#include <vector>
#include <map>

using namespace std;

struct FLAGindex{
    map<string, int> FLAG;
    
    FLAGindex(){
        FLAG["-GeneralParameters"]=0;
        FLAG["-Membrane"]=1;
        FLAG["-Actin"]=2;
        FLAG["-ECM"]=3;
        FLAG["-Chromatin"]=4;
        FLAG["-InteractionTable"]=5;
        
    }
};



struct GeneralParameters{
    map<string, vector<string> > GenParams;
    vector<string> insertOrder;
    vector<string> values;
    
    
    
    GeneralParameters(){
        values.resize(2);
        
        values[0] ="VCProject";
        values[1] ="#Set the Project Directory. Outputfile format: Results/ProjectName/Date_time_instanceId/Date_time_instanceId.extension";
        GenParams["ProjectName"] = values;
        insertOrder.push_back("ProjectName");
        
        values[0] ="10";
        values[1] ="#Simulation time length masured in pico seconds. Default value 10.";
        GenParams["SimulationTimeInPs"] = values;
        insertOrder.push_back("SimulationTimeInPs");
        
        values[0] ="10";
        values[1] ="#Integration step size measured in femto seconds.  Default value 10.";
        GenParams["StepSizeInFs"] = values;
        insertOrder.push_back("StepSizeInFs");
        
        values[0] ="1000";
        values[1] ="#Trajectory saving, data collection, etc wil be performed at the ReportIntervalInFs time intervales measured in femto seconds. Default value 1000.";
        GenParams["ReportIntervalInFs"] = values;
        insertOrder.push_back("ReportIntervalInFs");
        
        values[0] ="0";
        values[1] ="#Simulation enviroment size (cube). Default value 0.";
        GenParams["SimulationBoxLength"] = values;
        insertOrder.push_back("SimulationBoxLength");
        
        values[0] ="V";
        values[1] ="#Set the integrator type. 'Temperature' and 'FrictionCoeff' need to be set for Langevin and Brownian. V: Verlet. B: Brownian. L: Langevin. Default value 'V'.";
        GenParams["Integrator"] = values;
        insertOrder.push_back("Integrator");
        
        values[0] ="5";
        values[1] ="#The friction coefficient which couples the system to the heat bath (in inverse pico seconds). Required by the Brownian and Langevin integrators. Default value 5.";
        GenParams["FrictionInPs"] = values;
        insertOrder.push_back("FrictionInPs");
        
        values[0] ="310";
        values[1] ="#The thermostat temperature (in Kelvin). Required by the Brownian and Langevin integrators. Default value 310";
        GenParams["Temperature"] = values;
        insertOrder.push_back("Temperature");
        
        values[0] ="0";
        values[1] ="#The pressure acting on the system (in bar) through OpenMM's MonteCarloBarostat. Default 0";
        GenParams["MCBarostatPressure"] = values;
        insertOrder.push_back("MCBarostatPressure");
        
        values[0] ="0";
        values[1] ="#The temperature at which OpenMM's MonteCarloBarostat will think the system is being maintained (in Kelvin). Default thermostat's Temperature";
        GenParams["MCBarostatTemperature"] = values;
        insertOrder.push_back("MCBarostatTemperature");
        
        values[0] ="0";
        values[1] ="#The frequency at which the MonteCarloBarostat pressure changes should be attempted (in time steps). If zero MCBarostat will be disabled. Default 0";
        GenParams["MCBarostatFrequency"] = values;
        insertOrder.push_back("MCBarostatFrequency");
        
        values[0] ="0 0 0";
        values[1] ="#The pressure acting on each axis (in bar) through OpenMM's MonteCarloAnisotropicBarostat. Default 0 0 0";
        GenParams["MCAnisoBarostatPressure"] = values;
        insertOrder.push_back("MCAnisoBarostatPressure");
        
        values[0] ="0";
        values[1] ="#The temperature at which OpenMM's MonteCarloAnisotropicBarostat will think the system is being maintained (in Kelvin). Default thermostat's Temperature";
        GenParams["MCAnisoBarostatTemperature"] = values;
        insertOrder.push_back("MCAnisoBarostatTemperature");
        
        values[0] ="false false false";
        values[1] ="#Allow (or not allow) the X,Y, or Z dimension of the periodic box to change size through OpenMM's MonteCarloAnisotropicBarostat. Default fasle false false";
        GenParams["MCAnisoBarostatScaleXYZ"] = values;
        insertOrder.push_back("MCAnisoBarostatScaleXYZ");
        
        values[0] ="0";
        values[1] ="#The frequency at which the MonteCarloAnisotropicBarostat pressure changes should be attempted (in time steps). If zero MCAnisoBarostat will be disabled. Default 0";
        GenParams["MCAnisoBarostatFrequency"] = values;
        insertOrder.push_back("MCAnisoBarostatFrequency");
        
        values[0] ="0";
        values[1] ="#The number of steps that the centre of mass velocity is set to zero using OpenMM's CMMotionRemover. Default 0";
        GenParams["CMMotionRemoverStep"] = values;
        insertOrder.push_back("CMMotionRemoverStep");
        
        values[0] ="true";
        values[1] ="#Collect energy parameters for the potentials (expensive) during each Report_Interval_In_Fs time point. Default true";
        GenParams["ReportEnergy"] = values;
        insertOrder.push_back("ReportEnergy");
        
        values[0] ="false";
        values[1] ="#Write the particle velocities (cheap) fora each  ReportIntervalInFs time point. Default false";
        GenParams["WantVelocity"] = values;
        insertOrder.push_back("WantVelocity");
        
        values[0] ="1000 0 0";
        values[1] ="#Periodic box vector (1 of 3). Default value (1000, 0, 0).";
        GenParams["PeriodicBoxVector0"] = values;
        insertOrder.push_back("PeriodicBoxVector0");
        
        values[0] ="0 1000 0";
        values[1] ="#Periodic box vector (2 of 3). Default value (0, 1000, 0).";
        GenParams["PeriodicBoxVector1"] = values;
        insertOrder.push_back("PeriodicBoxVector1");
        
        values[0] ="0 0 1000";
        values[1] ="#Periodic box vector (3 of 3). Default value (0, 0, 1000).";
        GenParams["PeriodicBoxVector2"] = values;
        insertOrder.push_back("PeriodicBoxVector2");
        
        values[0] ="0";
        values[1] ="#This option is not available. Note to the developer: This should be moved to the Membrane. Default value 100.";
        GenParams["MCStep"] = values;
        insertOrder.push_back("MCStep");
        
        values[0] ="0";
        values[1] ="#This option is not available. Note to the developer: This should be moved to the Membrane. Default value 100.";
        GenParams["MemFluidity"] = values;
        insertOrder.push_back("MemFluidity");
        
    }
    
    map<string, vector<string> > get_map(){
        return GenParams;
    }
    
    vector<string> get_insertOrder(){
        return insertOrder;
    }
};

#endif /* OpenMM_structs_h */


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
    
    string ProjectName;
    /**A string containing the selected platform information during runtime.*/
    string hardwareReport;
    /**Name of the output file. Please note that the date and time the file is generated will be attached to this name.*/
    string trajectory_file_name;
    string Membrane_label="mem";
    string Actin_label="act";
    string ECM_label="ecm";
    string Chromatin_label="chr";
    
    bool MinimumForceDecleration;
    
    /**Boltzmann's constant set to 0.008314459920816468 KJ/mol.kelvin*/
    double BoltzmannKJpermolkelvin = 0.008314459920816468;
    /**Simulation time (in picoseconds). If this parameter is not set in the general config file by the user, or the value is set to zero, it will be calculate during runtime by multiplying the 'step size' by the 'total number of steps'.*/
    double Simulation_Time_In_Ps;
    /**Trajectory saving rate in femtoseconds. If this parameter is not set by the user in the general config file, or the value is set to zero, the interval will be calculated by multiplying the MD_traj_save_step by the Step_Size_In_Fs.*/
    double Report_Interval_In_Fs;
    /** Integration step size (fs).  Default 0.001*/
    double Step_Size_In_Fs;
    
    bool Minimise = false;
    double  MinimiseTolerance;
    int  MinimiseMaxIterations;
    
    int Num_of_Membranes;
    int Num_of_Actins;
    int Num_of_Chromatins;
    int Num_of_ECMs;
    
    bool   WantEnergy;
    /**Collect forces (cheap) during each Report_Interval_In_Fs time point. Default true*/
    bool   WantForce;
    /**Writes velocities and forces (cheap) of particles during each Report_Interval_In_Fs time point to the disk. Default false*/
    bool   WantVelocity;
    /**Make the velocity of the centre of mass (COM) zero by subtracting the COM velocity from all the particles' velocity after every CMMotionRemoverStep step. Default false*/
    bool CMMotionRemover;
    /**The number of steps where the centre of mass velocity is set to zero using the CMMotionRemover. Default 10*/
    int CMMotionRemoverStep;
    //1/0 true/false if you want/don't want to calculate and write the voronoi area associated with Membrane nodes in the properties output file. Default 0
    bool Wantvoronoi;
    //When in test mode, most of the console prints (std::cout) of the programme will be turned off for a better viewing of the test report. Deafault 0 (off)
    bool Testmode;
    
    /**The pressure acting on the system (in bar) through OpenMM's MonteCarloBarostat. Default 0*/
    double MCBarostatPressure;
    /**The temperature at which OpenMM's MonteCarloBarostat will think the system is being maintained (in Kelvin). Default Thermostat temperature*/
    double MCBarostatTemperature;
    /**The frequency at which the MonteCarloBarostat pressure changes should be attempted (in time steps). If zero MCBarostat will be disabled. Default 0.*/
    int MCBarostatFrequency;
    
    /**The status of OpenMM's MonteCarloAnisotropicBarostat. Default off*/
    bool MCAnisoBarostatOn;
    /**The pressure acting on each axis (in bar) through OpenMM's MonteCarloAnisotropicBarostat. Default 0 0 0*/
    vector<double> MCAnisoBarostatPressure;
    /**The temperature at which OpenMM's MonteCarloAnisotropicBarostat will think the system is being maintained (in Kelvin). Default Thermostat temperature*/
    double MCAnisoBarostatTemperature;
    /**Allow (or not allow) the X,Y, or Z dimension of the periodic box to change size through OpenMM's MonteCarloAnisotropicBarostat. Default fasle false false*/
    vector<bool> MCAnisoBarostatScaleXYZ;
    /**The frequency at which the MonteCarloAnisotropicBarostat pressure changes should be attempted (in time steps). If zero MCAnisoBarostat will be disabled. Default 0*/
    int MCAnisoBarostatFrequency;
    
    /**Set the integrator type. This flag is for the OpenMM integrators and will not function if the OpenMM engine is not selected.
     * V: Verlet
     * B: Brownian, temperature and frictionCoeff need to be set as well.
     * L: Langevin, temperature and frictionCoeff need to be set as well.
     * C: Langevin, temperature and frictionCoeff need to be set as well.
     *Default V*/
    string Integrator_type;
    /**Set the friction coefficient which couples the system to the heat bath (in inverse picoseconds). Default 5*/
    double frictionInPs;
    /**Set the temperature of the heat bath (in Kelvin). Default 300*/
    double temperature;
    
    /**Set the custom thermostat temperature of the heat bath (in Kelvin) by writing it after 'C' in the intergrator type category. If left empty the tempreture of the integrator will be used. */
    double customtemperature=-1;
    
    
    double Simulation_box_length;
    /**The simulation uses periodic boundary condition. Default False*/
    bool Periodic_condtion_status;
    vector<double> PeriodicBoxVector0;
    vector<double> PeriodicBoxVector1;
    vector<double> PeriodicBoxVector2;
    
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
        values[1] ="#Set the integrator type. 'Temperature' and 'FrictionCoeff' needs to be set for Langevin and Brownian. V: Verlet. B: Brownian. L: Langevin. Default value 'V'.";
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
        
        values[0] ="false";
        values[1] ="#Use 'true' or 'false' to specify th use of OpenMM's Minimize function. It searches for a new set of particle positions that represent a local minimum of the potential energy. The search is performed using the L-BFGS algorithm. Distance constraints are enforced during minimization by adding a harmonic restraining force to the potential function. Default false";
        GenParams["Minimise"] = values;
        insertOrder.push_back("Minimise");
        
        values[0] ="10";
        values[1] ="#This specifies how precisely the energy minimum must be located. Minimization will be halted once the root-mean-square value of all force components reaches this tolerance. The default value is 10.";
        GenParams["MinimiseTolerance"] = values;
        insertOrder.push_back("MinimiseTolerance");
        
        values[0] ="0";
        values[1] ="The maximum number of iterations to perform. If this is 0, minisation is continued until the results converge without regard to how many iterations it takes. The default value is 0.";
        GenParams["MinimiseMaxIterations"] = values;
        insertOrder.push_back("MinimiseMaxIterations");
        
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
        
        values[0] ="true";
        values[1] ="#When true, sll interaction potentails will be declared using the minimum number of 'Forces'. This will result in better performance for large systems. If false, multiple force groups will be defined for each class. This will come in handy when wanting to look at the evolution of selected forces on a class.";
        GenParams["MinimumForceDecleration"] = values;
        insertOrder.push_back("MinimumForceDecleration");
        
    }
    
    map<string, vector<string> > get_map(){
        return GenParams;
    }
    
    vector<string> get_insertOrder(){
        return insertOrder;
    }
};

#endif /* OpenMM_structs_h */


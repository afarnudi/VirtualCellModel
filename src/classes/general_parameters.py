from classes.general.settings import Settings


class GeneralParameters:
    default_settings = {
        "ProjectName": Settings(
            "VCMProject",
            "#Set the Project Directory. The ProjectName can be a single "
            "directory or a path. Output file format: "
            "Results/ProjectName/Date_time_instanceId/Date_time_instanceId.extension"
            "Default value: VCMProject.",
        ),
        "SimulationTimeInPs": Settings(
            "10",
            "#Simulation total runtime duration measured in picoseconds. "
            "Floating point values are allowed. Default value 10.",
        ),
        "StepSizeInFs": Settings(
            "10",
            "#Integration step size measured in femtoseconds. Floating point"
            " values are allowed. Default value: 10.",
        ),
        "ReportIntervalInFs": Settings(
            "1000",
            "#Trajectory saving, data collection, etc wil be performed at the"
            " ReportIntervalInFs time intervals measured in femtoseconds.\n"
            "It can also be set to record data at exponential time intervals "
            "by specifying 'EXP Step' followed by the number of samples you "
            "want to save. Example: ReportIntervalInFs EXP Step 200, the "
            "simulation outputs will be recorded every t = int(exp(m*delta)),"
            " where t is an integer and the maximum value is "
            "T=SimulationTimeInPs*1000/StepSizeInFs, m is an integer "
            "multiplier, and delta=log(T)/N. Note log(T) is the natural "
            "logarithm and N is the number of report intervals requested "
            "(here 200). \nAlso, you can just set an exponent for the "
            "exponential data sampling. Example ReportIntervalInFs EXP 0.98, "
            "the simulation outputs will be recorded every t = exp(m*delta), "
            "where t is an integer and the maximum value is "
            "T=SimulationTimeInPs*1000/StepSizeInFs, m is an integer "
            "multiplier. Note log(T) is the natural logarithm. "
            "Default value: 1000",
        ),
        "SimulationBoxLength": Settings(
            "0",
            "#Simulation box size (cube). When value is zero, the periodic "
            "boundary condition will be switched off. If value is not zero, "
            "the simulation will be performed in a periodic cube with length "
            'SimulationBoxLength. Set the "PeriodicBoxVector"s when '
            "simulating a non-cubic box with periodic boundary condition. "
            "Default value: 0.",
        ),
        "Integrator": Settings(
            "V",
            "#Set the integrator type. Default value 'V'.\n"
            "#\tV: Verlet (default integrator). \n"
            "#\tB: Brownian. *needs: FrictionInInvertPs and Temperature\n"
            "#\tL: Langevin. *needs: FrictionInInvertPs and Temperature\n"
            "#\tLM: This is an Integrator which simulates a System using "
            "Langevin dynamics, with the LFMiddle discretization. *needs: "
            "FrictionInInvertPs and Temperature\n"
            '#\tGJF: The GJF thermostat based on "A simple and effective '
            'Verlet-type algorithm for simulating Langevin dynamics" by Niels'
            " Grønbech-Jensen  & Oded Farago Published online: 14 Feb 2013 "
            "DOI:10.1080/00268976.2012.760055  *needs: FrictionInInvertPs "
            "and Temperature\n"
            "#\tGJF20: Use: GJF20 A (or B, C) An update to the GJF thermostat"
            ' based on "Defining velocities for accurate kinetic statistics '
            'in the GJF thermostat" by Niels Grønbech-Jensen and Oded Farago '
            "DOI: 10.1103/PhysRevE.101.022123 *needs: FrictionInInvertPs and "
            "Temperature\n"
            "#\tGlobal: Bussi et al global thermostat DOI: "
            "http://dx.doi.org/10.1016/j.cpc.2008.01.006 *needs: "
            "FrictionInInvertPs and Temperature\n",
        ),
        "FrictionInInvertPs": Settings(
            "0.01",
            "#The friction coefficient which couples the system to the heat "
            "bath (in inverse picoseconds). Required by the Brownian and "
            "Langevin integrators. Default value: 0.01",
        ),
        "TemperatureInKelvin": Settings(
            "310",
            "#The thermostat temperature (in Kelvin). Most integrators "
            "require a temperature for particle velocity adjustment. Default"
            " value: 310",
        ),
        "setVelocitiesToTemperature": Settings(
            "false",
            "#Set all particle initial velocities to random values taken from"
            " a Boltzmann distribution at a given temperature. Default value:"
            " false",
        ),
        "Seed": Settings(
            "0",
            "Most integrators require to generate sudo random numbers (SRN) "
            "to add noise to the system. Each seed can be used to generate a "
            "unique set of SRN. Here is more information on how OpenMM uses "
            'the seed for integration. #OpenMM Documentation: "Set the '
            "random number seed. The precise meaning of this parameter is "
            "undefined, and is left up to each Platform to interpret in an "
            "appropriate way. It is guaranteed that if two simulations are "
            "run with different random number seeds, the sequence of random "
            "forces will be different. On the other hand, no guarantees are "
            "made about the behavior of simulations that use the same seed. "
            "In particular, Platforms are permitted to use non-deterministic "
            "algorithms which produce different results on successive runs, "
            'even if those runs were initialized identically.". The quote '
            "was taken from the OpenMM's BrownianIntegrator documentation on "
            "http://docs.openmm.org/latest/api-python/generated/openmm.openmm.BrownianIntegrator.html#openmm.openmm.BrownianIntegrator.setRandomNumberSeed"
            " Default value: 0",
        ),
        "Minimise": Settings(
            "false",
            "#If 'true', use OpenMM's Minimize function. Search for a new set"
            " of particle positions that corresponding to a local minimum of "
            "the potential energy. The search uses the L-BFGS algorithm. "
            "Distance constraints are enforced during minimization by adding "
            "a harmonic restraining force to the potential function. Default "
            "value: false",
        ),
        "MinimiseTolerance": Settings(
            "10",
            "#Specify how precisely the energy minimum must be located. "
            "Minimisation will be halted once the root-mean-square value of "
            "all force components reaches this tolerance. Default value: 10.",
        ),
        "MinimiseMaxIterations": Settings(
            "0",
            "The maximum number of iterations to perform during minimisation."
            " If zero, minisation is continued until the results converge "
            "without regard to how many iterations it takes. Default value: 0",
        ),
        "MCBarostatPressure": Settings(
            "0",
            "#The pressure acting on the system (in bar) through OpenMM's "
            "MonteCarloBarostat. Default value: 0",
        ),
        "MCBarostatTemperature": Settings(
            "TemperatureInKelvin",
            "#The temperature at which OpenMM's MonteCarloBarostat will think"
            " the system is being maintained (in Kelvin). "
            "Default value: TemperatureInKelvin",
        ),
        "MCBarostatFrequency": Settings(
            "0",
            "#The attempt frequency of the MonteCarloBarostat pressure "
            "changes (in time steps). If zero, MCBarostat will be disabled. "
            "Default value: 0",
        ),
        "MCAnisoBarostatPressure": Settings(
            "0 0 0",
            "#The pressure acting on each axis (in bar) through OpenMM's "
            "MonteCarloAnisotropicBarostat. Default value: 0 0 0",
        ),
        "MCAnisoBarostatTemperature": Settings(
            "TemperatureInKelvin",
            "#The temperature at which OpenMM's MonteCarloAnisotropicBarostat"
            " will think the system is being maintained (in Kelvin). "
            "Default value: TemperatureInKelvin",
        ),
        "MCAnisoBarostatScaleXYZ": Settings(
            "false false false",
            "#Allow (or not allow) the X,Y, or Z dimension of the periodic "
            "box to change size through OpenMM's "
            "MonteCarloAnisotropicBarostat. Default value: fasle false false",
        ),
        "MCAnisoBarostatFrequency": Settings(
            "0",
            "#The attempt frequency which MonteCarloAnisotropicBarostat "
            "pressure changes should be attempted (in time steps). If zero,"
            " MCAnisoBarostat will be disabled. "
            "Default value: 0",
        ),
        "CMMotionRemoverStep": Settings(
            "0",
            "#The number of steps that the centre of mass velocity is set to "
            "zero using OpenMM's CMMotionRemover. "
            "Default value: 0",
        ),
        "ReportEnergy": Settings(
            "true",
            "#Collect energy parameters for the potentials (expensive) "
            "during each Report_Interval_In_Fs time point. "
            "Default value: true",
        ),
        "PeriodicBoxVector0": Settings(
            "1000 0 0",
            "#Periodic box vector (1 of 3). Default value: 1000, 0, 0",
        ),
        "PeriodicBoxVector1": Settings(
            "0 1000 0",
            "#Periodic box vector (2 of 3). Default value: 0, 1000, 0",
        ),
        "PeriodicBoxVector2": Settings(
            "0 0 1000",
            "#Periodic box vector (3 of 3). Default value: 0, 0, 1000",
        ),
        "TextOutputs": Settings(
            "PSF XYZ",
            "#Simulation outputs (text format): You can add as many output"
            " formats as you wish. PSF: psf connectivity file, XYZ: xyz file"
            ", PDB: pdb file, VEL: xyz file style for velocities, FORCE: xyz"
            " file style for forces, N: no text output. "
            "Default value: PSF XYZ",
        ),
        "BinOutputs": Settings(
            "N",
            "#Simulation outputs (binary format): You can add as many output"
            " formats as you wish. XYZ: xyz coordinates, VEL: xyz components "
            "of velocities, TPK: time, potential energy, and kinetic energy, "
            "N: No binary output. "
            "Default value: N ",
        ),
        "Precision": Settings(
            "single",
            "#Set the precision of the calculations. If the precision is not"
            " supported by the platform the default precision of the "
            "platform will be used. "
            "Default value: single ",
        ),
        "bondCutoff": Settings(
            "0",
            "#For non bonded forces, pairs of particles that are separated "
            "by this many bonds or fewer are added to the list of exclusions."
            " Default value: 0",
        ),
    }

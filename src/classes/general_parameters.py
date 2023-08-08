import numpy as np
from classes.general.configuration import Settings
from classes.general.print_colors import TerminalColors as tc
from general.datatype_conversion import string_to_float
from general.datatype_conversion import string_to_int
from general.datatype_conversion import string_to_bool


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
        "SetVelocitiesToTemperature": Settings(
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
        "BondCutoff": Settings(
            "0",
            "#For non bonded forces, pairs of particles that are separated "
            "by this many bonds or fewer are added to the list of exclusions."
            " Default value: 0",
        ),
    }
    caller_title = "GeneralParameters parser"

    def __init__(self):
        self.settings = GeneralParameters.default_settings.copy()
        self.lower_to_cap_key_dict = None
        self.build_lower_to_cap_key_dict()

        self.simulation_time_in_ps = None
        self.exp_sampling = None
        self.exp_sampling_exponent = None
        self.exp_sampling_num_of_steps = None
        self.report_interval_in_fs = None
        self.step_size_in_fs = None
        self.simulation_box_length = None
        self.periodic_boundary_condition = None
        self.available_integrators = [
            "V",
            "B",
            "M",
            "L",
            "LM",
            "LFLangevinMultiT",
            "LFLangevinMultiTDropN3",
            "GJF",
            "GJFMT",
            "GJFMTDropN3",
            "GJF20",
            "Global",
            "GlobalMTDropN3",
        ]
        self.integrator = None
        self.minimisation_integrator_restriction = None
        self.custom_temperature = None
        self.GJF_case = None
        self.friction_in_invert_ps = None
        self.seed = None
        self.bond_cut_off = None
        self.temperature = None
        self.want_energy = None
        self.set_velocities_to_temperature = None
        self.minimise = None
        self.minimise_tolerance = None
        self.minimise_max_iterations = None
        self.cmm_motion_remover_step = None
        self.mc_barostat_pressure = None
        self.mc_barostat_temperature = None
        self.mc_barostat_frequency = None
        self.mc_aniso_barostat_pressure = None
        self.mc_aniso_barostat_temperature = None
        self.mc_aniso_barostat_scale_xyz = None
        self.mc_aniso_barostat_frequency = None
        self.mc_aniso_barostat_is_on = False
        self.mc_aniso_barostat_scale_is_on = False
        self.periodic_box_vector0 = None
        self.periodic_box_vector1 = None
        self.periodic_box_vector2 = None
        self.project_name = None
        self.want_psf = None
        self.want_pdb = None
        self.want_xyz = None
        self.want_vel = None
        self.want_force = None
        self.want_curvature = None
        self.want_xy_bin = None
        self.want_vel_bin = None
        self.want_tpk_bin = None
        self.precision = None

    def build_lower_to_cap_key_dict(self):
        cap_keys = self.settings.keys()
        self.lower_to_cap_key_dict = {key.lower(): key for key in cap_keys}

    def set_value(self, key, val):
        keys = self.lower_to_cap_key_dict.keys()
        if key.lower() in keys:
            cap_key = self.lower_to_cap_key_dict[key.lower()]
            self.settings[cap_key].value = val
        else:
            raise RuntimeError(
                f'{tc.TFAILED} "{tc.TFILE}{key}{tc.TFAILED}" is not a GeneralParameters setting. Use the template generator to get a list of available settings.{tc.TRESET}'
            )

    def parse_settings(self):
        self.parse_values()
        self.run_consistency_check()

    def msg_title(self, key):
        return f"{GeneralParameters.caller_title}: {key}:"

    def parse_values(self):
        keys = self.default_settings.keys()
        for key in keys:
            vals = self.settings[key].value.split()
            help = self.settings[key].help
            msg = self.msg_title(key)

            if key == "SimulationTimeInPs":
                self.simulation_time_in_ps = string_to_float(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "ReportIntervalInFs":
                if vals[0] == "EXP":
                    self.exp_sampling = True
                    if len(vals) < 2:
                        raise TypeError(
                            f'{tc.TWARN}{msg} You need to specify an exponent. Example: ReportIntervalInFs EXP 12.2 or specify the number of steps to be saved: EXP Step 200"{tc.TRESET}'
                        )
                    else:
                        if vals[1] == "Step":
                            if len(vals) < 3:
                                raise TypeError(
                                    f'{tc.TWARN}{msg} You need to specify the number of steps. Example: ReportIntervalInFs EXP Step 200"{tc.TRESET}'
                                )
                            self.exp_sampling_num_of_steps = string_to_int(
                                vals[2],
                                msg,
                                help,
                            )
                        else:
                            self.exp_sampling_exponent = string_to_float(
                                vals[1],
                                msg,
                                help,
                            )
                else:
                    self.exp_sampling = False
                    self.report_interval_in_fs = string_to_float(
                        vals[0],
                        msg,
                        help,
                    )
            elif key == "StepSizeInFs":
                self.step_size_in_fs = string_to_float(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "SimulationBoxLength":
                self.simulation_box_length = string_to_float(
                    vals[0],
                    msg,
                    help,
                )
                if vals[0] != "0":
                    self.periodic_boundary_condition = True
                else:
                    self.periodic_boundary_condition = False
            elif key == "Integrator":
                if vals[0] not in self.available_integrators:
                    raise TypeError(
                        f'{tc.TFAILED}GeneralParameters parser: Integrator: I don\'t understand integrator option. "{tc.TFILE}{vals[0]}{tc.TFAILED}". Please consult the setting description:\n{tc.TRESET}{help}'
                    )
                if vals[0] == "V":
                    self.integrator = "Verlet"
                elif vals[0] == "B":
                    self.integrator = "Brownian"
                elif vals[0] == "M":
                    self.integrator = "LangevinMinimise"
                    if len(vals) != 2:
                        raise TypeError(
                            f'{tc.TFAILED}GeneralParameters parser: Integrator: You need to specify a restriction for the position evolution for the "Langevin minimise" integrator. Example: Integrator M 12.2857142857. Please consult the setting description:\n{tc.TRESET}{help}'
                        )
                    self.minimisation_integrator_restriction = string_to_float(
                        vals[1],
                        msg,
                        help,
                    )
                elif vals[0] == "L":
                    self.integrator = "Langevin"
                elif vals[0] == "LM":
                    self.integrator = "LangevinMiddle"
                elif vals[0] == "LFLangevinMultiT":
                    self.integrator = "LFLangevinMulti-thermos"
                    if len(vals) > 1:
                        self.custom_temperature = string_to_float(
                            vals[1],
                            msg,
                            help,
                        )
                elif vals[0] == "LFLangevinMultiTDropN3":
                    self.integrator = "LFLangevinMulti-thermosDropNewton3"
                    if len(vals) > 1:
                        self.custom_temperature = string_to_float(
                            vals[1],
                            msg,
                            help,
                        )
                elif vals[0] == "GJF":
                    self.integrator = "GJF"
                elif vals[0] == "GJFMT":
                    self.integrator = "GJF2013Multi-thermos"
                    if len(vals) > 1:
                        self.custom_temperature = string_to_float(
                            vals[1],
                            msg,
                            help,
                        )
                elif vals[0] == "GJFMTDropN3":
                    self.integrator = "GJF2013Multi-thermosDropNewton3"
                    if len(vals) > 1:
                        self.custom_temperature = string_to_float(
                            vals[1],
                            msg,
                            help,
                        )
                elif vals[0] == "GJF20":
                    self.integrator = "GJF2020"
                    if len(vals) != 2:
                        raise TypeError(
                            f"{tc.TFAILED}GeneralParameters parser: Integrator: You need to specify the 'Case'. Choices are between A and B. Example: Integrator GJF20 A. Please consult the setting description:\n{tc.TRESET}{help}"
                        )
                    self.GJF_case = vals[1]
                elif vals[0] == "Global":
                    self.integrator = "Bussi2008"
                elif vals[0] == "GlobalMTDropN3":
                    self.integrator = "Bussi2008Multi-thermosDropNewton3"
                    if len(vals) > 1:
                        self.custom_temperature = string_to_float(
                            vals[1],
                            msg,
                            help,
                        )
            elif key == "FrictionInInvertPs":
                self.friction_in_invert_ps = string_to_float(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "Seed":
                self.seed = string_to_int(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "BondCutoff":
                self.bond_cut_off = string_to_int(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "TemperatureInKelvin":
                self.temperature = string_to_float(
                    vals[0],
                    msg,
                    help,
                )
                if self.custom_temperature is None:
                    self.custom_temperature = self.temperature
            elif key == "ReportEnergy":
                self.want_energy = string_to_bool(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "SetVelocitiesToTemperature":
                self.set_velocities_to_temperature = string_to_bool(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "Minimise":
                self.minimise = string_to_bool(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "MinimiseTolerance":
                self.minimise_tolerance = string_to_float(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "MinimiseMaxIterations":
                self.minimise_max_iterations = string_to_int(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "CMMotionRemoverStep":
                self.cmm_motion_remover_step = string_to_int(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "MCBarostatPressure":
                self.mc_barostat_pressure = string_to_float(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "MCBarostatTemperature":
                if vals[0] == "TemperatureInKelvin":
                    self.mc_barostat_temperature = self.temperature
                else:
                    self.mc_barostat_temperature = string_to_float(
                        vals[0],
                        msg,
                        help,
                    )
            elif key == "MCBarostatFrequency":
                self.mc_barostat_frequency = string_to_float(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "MCAnisoBarostatPressure":
                if len(vals) < 3:
                    raise TypeError(
                        f"{tc.TFAILED}GeneralParameters parser: MCAnisoBarostatPressure: Three arguments required. Please consult the setting description:\n{tc.TRESET}{help}"
                    )
                self.mc_aniso_barostat_pressure = [
                    string_to_float(
                        vals[0],
                        msg,
                        help,
                    ),
                    string_to_float(
                        vals[1],
                        msg,
                        help,
                    ),
                    string_to_float(
                        vals[2],
                        msg,
                        help,
                    ),
                ]
            elif key == "MCAnisoBarostatTemperature":
                if vals[0] == "TemperatureInKelvin":
                    self.mc_aniso_barostat_temperature = self.temperature
                else:
                    self.mc_aniso_barostat_temperature = string_to_float(
                        vals[0],
                        msg,
                        help,
                    )
            elif key == "MCAnisoBarostatScaleXYZ":
                if len(vals) < 3:
                    raise TypeError(
                        f"{tc.TFAILED}GeneralParameters parser: MCAnisoBarostatPressure: Three arguments required. Please consult the setting description:\n{tc.TRESET}{help}"
                    )
                self.mc_aniso_barostat_scale_xyz = [
                    string_to_bool(
                        vals[0],
                        msg,
                        help,
                    ),
                    string_to_bool(
                        vals[1],
                        msg,
                        help,
                    ),
                    string_to_bool(
                        vals[2],
                        msg,
                        help,
                    ),
                ]
            elif key == "MCAnisoBarostatFrequency":
                self.mc_aniso_barostat_frequency = string_to_float(
                    vals[0],
                    msg,
                    help,
                )
            elif key == "PeriodicBoxVector0":
                if len(vals) < 3:
                    raise TypeError(
                        f"{tc.TFAILED}GeneralParameters parser: PeriodicBoxVector0: Three arguments required. Please consult the setting description:\n{tc.TRESET}{help}"
                    )
                self.periodic_box_vector0 = [
                    string_to_float(
                        vals[0],
                        msg,
                        help,
                    ),
                    string_to_float(
                        vals[1],
                        msg,
                        help,
                    ),
                    string_to_float(
                        vals[2],
                        msg,
                        help,
                    ),
                ]
            elif key == "PeriodicBoxVector1":
                if len(vals) < 3:
                    raise TypeError(
                        f"{tc.TFAILED}GeneralParameters parser: PeriodicBoxVector1: Three arguments required. Please consult the setting description:\n{tc.TRESET}{help}"
                    )
                self.periodic_box_vector1 = [
                    string_to_float(
                        vals[0],
                        msg,
                        help,
                    ),
                    string_to_float(
                        vals[1],
                        msg,
                        help,
                    ),
                    string_to_float(
                        vals[2],
                        msg,
                        help,
                    ),
                ]
            elif key == "PeriodicBoxVector2":
                if len(vals) < 3:
                    raise TypeError(
                        f"{tc.TFAILED}GeneralParameters parser: PeriodicBoxVector2: Three arguments required. Please consult the setting description:\n{tc.TRESET}{help}"
                    )
                self.periodic_box_vector2 = [
                    string_to_float(
                        vals[0],
                        msg,
                        help,
                    ),
                    string_to_float(
                        vals[1],
                        msg,
                        help,
                    ),
                    string_to_float(
                        vals[2],
                        msg,
                        help,
                    ),
                ]
            elif key == "ProjectName":
                self.project_name = vals[0]
            elif key == "TextOutputs":
                for val in vals:
                    if val in ["PSF", "PDB", "XYZ", "VEL", "FORCE", "CURVE"]:
                        if val == "PSF":
                            self.want_psf = True
                        elif val == "PDB":
                            self.want_pdb = True
                        elif val == "XYZ":
                            self.want_xyz = True
                        elif val == "VEL":
                            self.want_vel = True
                        elif val == "FORCE":
                            self.want_force = True
                        elif val == "CURVE":
                            self.want_curvature = True
                    else:
                        raise TypeError(
                            f'{tc.TFAILED}GeneralParameters parser: TextOutputs: I don\'t understand TextOutputs option. "{tc.TFILE}{val}{tc.TFAILED}". Please consult the setting description:\n{tc.TRESET}{help}'
                        )
            elif key == "BinOutputs":
                for val in vals:
                    if val in ["XYZ", "VEL", "TPK"]:
                        if val == "XYZ":
                            self.want_xyz_bin = True
                        elif val == "VEL":
                            self.want_vel_bin = True
                        elif val == "TPK":
                            self.want_tpk_bin = True
                    else:
                        raise TypeError(
                            f'{tc.TFAILED}GeneralParameters parser: BinOutputs: I don\'t understand BinOutputs option. "{tc.TFILE}{val}{tc.TFAILED}". Please consult the setting description:\n{tc.TRESET}{help}'
                        )
            elif key == "Precision":
                if vals[0].lower() in ["single", "double"]:
                    self.precision = vals[0]
                else:
                    raise TypeError(
                        f'{tc.TFAILED}GeneralParameters parser: Precision: I don\'t understand Precision option. "{tc.TFILE}{vals[0]}{tc.TFAILED}". Choices are "single" or "double". Please consult the setting description:\n{tc.TRESET}{help}'
                    )

    def run_consistency_check(self):
        for i in range(3):
            if np.isclose(self.mc_aniso_barostat_pressure[i], 0):
                self.mc_aniso_barostat_is_on = True
            if self.mc_aniso_barostat_scale_xyz[i]:
                self.mc_aniso_barostat_scale_is_on = True
        if self.mc_aniso_barostat_is_on and self.mc_aniso_barostat_scale_is_on is False:
            raise RuntimeError(
                f"{tc.TWARN}GeneralParameters parser: consistency check: MCAnisoBarostat and MCAnisoBarostatScaleXYZ: Non of the axes are allowed to scale. At least one axis has to be allowed to scale. Please edit the configuration file and try again.\n"
            )
        if (
            self.mc_aniso_barostat_is_on
            and np.isclose(self.mc_barostat_pressure, 0) is False
        ):
            raise RuntimeError(
                f"{tc.TWARN}GeneralParameters parser: consistency check: MCAnisoBarostat and MCBarostatPressure: VCM does not have a feature that allows the MCAnisoBarostat and MCBarostat to work simultaneously. Please choose one of the methods and edit the configurations and try again.\n"
            )

import numpy as np
import pytest
from classes.general_parameters import GeneralParameters


def test_gen_param_parsing_defaults():
    gen_param = GeneralParameters()
    gen_param.parse_settings()
    assert gen_param.project_name == "VCMProject"
    assert np.isclose(gen_param.simulation_time_in_ps, 10.0)
    assert np.isclose(gen_param.step_size_in_fs, 10.0)
    assert np.isclose(gen_param.report_interval_in_fs, 1000.0)
    assert gen_param.exp_sampling == False
    assert gen_param.exp_sampling_exponent == None
    assert gen_param.exp_sampling_num_of_steps == None
    assert np.isclose(gen_param.simulation_box_length_in_nm, 0.0)
    assert gen_param.periodic_boundary_condition == False
    assert gen_param.integrator == "Verlet"
    assert gen_param.minimisation_integrator_restriction == None
    assert gen_param.custom_temperature == gen_param.temperature
    assert gen_param.GJF_case == None
    assert np.isclose(gen_param.friction_in_invert_ps, 0.01)
    assert np.isclose(gen_param.temperature, 310.0)
    assert gen_param.set_velocities_to_temperature == False
    assert gen_param.seed == 0
    assert type(gen_param.seed) == int
    assert gen_param.minimise == False
    assert np.isclose(gen_param.minimise_tolerance, 10.0)
    assert gen_param.minimise_max_iterations == 0
    assert type(gen_param.minimise_max_iterations) == int
    assert np.isclose(gen_param.mc_barostat_pressure, 0.0)
    assert np.isclose(gen_param.mc_barostat_temperature, gen_param.temperature)
    assert np.isclose(gen_param.mc_barostat_frequency, 0.0)
    assert gen_param.mc_aniso_barostat_pressure == [0.0, 0.0, 0.0]
    assert np.isclose(gen_param.mc_aniso_barostat_temperature, gen_param.temperature)
    assert gen_param.mc_aniso_barostat_scale_xyz == [False, False, False]
    assert gen_param.mc_aniso_barostat_is_on == False
    assert np.isclose(gen_param.mc_aniso_barostat_frequency, 0.0)
    assert gen_param.cmm_motion_remover_step == 0
    assert type(gen_param.cmm_motion_remover_step) == int
    assert gen_param.want_energy == True
    assert gen_param.periodic_box_vector0 == [1000.0, 0.0, 0.0]
    assert gen_param.periodic_box_vector1 == [0.0, 1000.0, 0.0]
    assert gen_param.periodic_box_vector2 == [0.0, 0.0, 1000.0]
    assert gen_param.want_curvature == False
    assert gen_param.want_vel == False
    assert gen_param.want_psf == True
    assert gen_param.want_xyz == True
    assert gen_param.want_force == False
    assert gen_param.want_pdb == False
    assert gen_param.want_tpk_bin == False
    assert gen_param.want_vel_bin == False
    assert gen_param.want_xyz_bin == False
    assert gen_param.precision == "single"
    assert gen_param.bond_cut_off == 0
    assert type(gen_param.bond_cut_off) == int


def test_gen_param_parsing_non_defaults():
    gen_param = GeneralParameters()
    gen_param.set_value("ProjectName", "My/path/to/project")
    gen_param.set_value("SimulationTimeInPs", "100")
    gen_param.set_value("StepSizeInFs", "100")
    gen_param.set_value("ReportIntervalInFs", "10")
    gen_param.set_value("SimulationBoxLengthInNm", "10")
    gen_param.set_value("Integrator", "LM")
    gen_param.set_value("FrictionInInvertPs", "100")
    gen_param.set_value("TemperatureInKelvin", "100")
    gen_param.set_value("SetVelocitiesToTemperature", "true")
    gen_param.set_value("Seed", "123")
    gen_param.set_value("Minimise", "True")
    gen_param.set_value("MinimiseTolerance", "0.1")
    gen_param.set_value("MinimiseMaxIterations", "1000")
    gen_param.set_value("MCBarostatPressure", "1000")
    gen_param.set_value("MCBarostatTemperature", "2")
    gen_param.set_value("MCBarostatFrequency", "20")
    gen_param.set_value("MCAnisoBarostatPressure", "0 2 -200")
    gen_param.set_value("MCAnisoBarostatTemperature", "-5")
    gen_param.set_value("MCAnisoBarostatScaleXYZ", "false True true")
    gen_param.set_value("MCAnisoBarostatFrequency", "3")
    gen_param.set_value("CMMotionRemoverStep", "30")
    gen_param.set_value("ReportEnergy", "False")
    gen_param.set_value("PeriodicBoxVector0", "12 34 56")
    gen_param.set_value("PeriodicBoxVector1", "-78 910.11 -12.13")
    gen_param.set_value("PeriodicBoxVector2", "14 -151.6 -.1718")
    gen_param.set_value("TextOutputs", "N")
    gen_param.set_value("BinOutputs", "XYZ VEL TPK")
    gen_param.set_value("Precision", "double")
    gen_param.set_value("BondCutoff", "3")

    gen_param.parse_settings()

    assert gen_param.project_name == "My/path/to/project"
    assert np.isclose(gen_param.simulation_time_in_ps, 100.0)
    assert np.isclose(gen_param.step_size_in_fs, 100.0)
    assert np.isclose(gen_param.report_interval_in_fs, 10.0)
    assert gen_param.exp_sampling == False
    assert gen_param.exp_sampling_exponent == None
    assert gen_param.exp_sampling_num_of_steps == None
    assert np.isclose(gen_param.simulation_box_length_in_nm, 10.0)
    assert gen_param.periodic_boundary_condition == True
    assert gen_param.integrator == "LangevinMiddle"
    assert gen_param.minimisation_integrator_restriction == None
    assert gen_param.custom_temperature == gen_param.temperature
    assert gen_param.GJF_case == None
    assert np.isclose(gen_param.friction_in_invert_ps, 100)
    assert np.isclose(gen_param.temperature, 100.0)
    assert gen_param.set_velocities_to_temperature == True
    assert gen_param.seed == 123
    assert type(gen_param.seed) == int
    assert gen_param.minimise == True
    assert np.isclose(gen_param.minimise_tolerance, 0.1)
    assert gen_param.minimise_max_iterations == 1000
    assert type(gen_param.minimise_max_iterations) == int
    assert np.isclose(gen_param.mc_barostat_pressure, 1000.0)
    assert np.isclose(gen_param.mc_barostat_temperature, 2)
    assert np.isclose(gen_param.mc_barostat_frequency, 20.0)
    assert gen_param.mc_aniso_barostat_pressure == [0.0, 2.0, -200.0]
    assert np.isclose(gen_param.mc_aniso_barostat_temperature, -5)
    assert gen_param.mc_aniso_barostat_scale_xyz == [False, True, True]
    assert gen_param.mc_aniso_barostat_is_on == False
    assert np.isclose(gen_param.mc_aniso_barostat_frequency, 3.0)
    assert gen_param.cmm_motion_remover_step == 30
    assert type(gen_param.cmm_motion_remover_step) == int
    assert gen_param.want_energy == False
    assert gen_param.periodic_box_vector0 == [12.0, 34.0, 56.0]
    assert gen_param.periodic_box_vector1 == [-78.0, 910.11, -12.13]
    assert gen_param.periodic_box_vector2 == [14.0, -151.6, -0.1718]
    assert gen_param.want_curvature == False
    assert gen_param.want_vel == False
    assert gen_param.want_psf == False
    assert gen_param.want_xyz == False
    assert gen_param.want_force == False
    assert gen_param.want_pdb == False
    assert gen_param.want_tpk_bin == True
    assert gen_param.want_vel_bin == True
    assert gen_param.want_xyz_bin == True
    assert gen_param.precision == "double"
    assert gen_param.bond_cut_off == 3
    assert type(gen_param.bond_cut_off) == int


def test_gen_param_parsing_ReportIntervalInFs_EXP():
    gen_param = GeneralParameters()
    gen_param.set_value("ReportIntervalInFs", "EXP")
    with pytest.raises(TypeError) as exc_info:
        gen_param.parse_settings()
        assert (
            "GeneralParameters parser: ReportIntervalInFs: You need to specify an exponent. Example: ReportIntervalInFs EXP 12.2 or specify the number of steps to be saved: EXP Step 200"
            in str(exc_info)
        )


def test_gen_param_parsing_ReportIntervalInFs_EXP_12():
    gen_param = GeneralParameters()
    gen_param.set_value("ReportIntervalInFs", "EXP 12.2")
    gen_param.parse_settings()
    assert gen_param.exp_sampling == True
    assert gen_param.exp_sampling_exponent == 12.2
    assert gen_param.exp_sampling_num_of_steps is None


def test_gen_param_parsing_ReportIntervalInFs_EXP_step():
    gen_param = GeneralParameters()
    gen_param.set_value("ReportIntervalInFs", "EXP Step")
    with pytest.raises(TypeError) as exp_info:
        gen_param.parse_settings()
        assert (
            "GeneralParameters parser: ReportIntervalInFs: You need to specify the number of steps. Example: ReportIntervalInFs EXP Step 200"
            in str(exp_info)
        )


def test_gen_param_parsing_ReportIntervalInFs_EXP_step_12():
    gen_param = GeneralParameters()
    gen_param.set_value("ReportIntervalInFs", "EXP step 12")
    with pytest.raises(TypeError) as exc_info:
        gen_param.parse_settings()
        assert (
            '" when expected "Step". You need to specify the number of steps. Example: ReportIntervalInFs EXP Step 200'
            in str(exc_info)
        )


def test_gen_param_parsing_ReportIntervalInFs_EXP_Step_12_point_3():
    gen_param = GeneralParameters()
    gen_param.set_value("ReportIntervalInFs", "EXP Step 12.3")
    with pytest.raises(TypeError) as exc_info:
        gen_param.parse_settings()
        assert (
            "GeneralParameters parser: ReportIntervalInFs: Was not able to convert "
            in str(exc_info)
        )


def test_gen_param_parsing_TemperatureInKelvin_negative_values():
    gen_param = GeneralParameters()
    gen_param.set_value("TemperatureInKelvin", "-1")
    with pytest.raises(ValueError) as exc_info:
        gen_param.parse_settings()
        assert (
            "GeneralParameters parser: TemperatureInKelvin: Negative values not allowed for the system temperature."
            in str(exc_info)
        )


def test_gen_param_parsing_TemperatureInKelvin_two_values():
    gen_param = GeneralParameters()
    gen_param.set_value("TemperatureInKelvin", "10 100")
    gen_param.parse_settings()
    assert gen_param.temperature == 10
    assert type(gen_param.temperature) == float

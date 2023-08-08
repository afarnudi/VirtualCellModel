import numpy as np
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

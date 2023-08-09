from classes.general.print_colors import TerminalColors as tc
from classes.general_parameters import GeneralParameters


def get_section_configs(configs, sec_name):
    conf_dicts = []
    for conf in configs:
        if conf.get_name() == sec_name:
            conf_dicts.append(conf)
    return conf_dicts


def set_general_parameter_configs(gp_input):
    gp = GeneralParameters()
    input_keys = gp_input.get_key_names()
    for key in input_keys:
        val = gp_input.get_value(key)
        gp.set_value(key, val)
    return gp


def get_general_parameter_section(gp_inputs):
    if len(gp_inputs) != 1:
        raise RuntimeError(
            f"{tc.TFAILED}One (and only one) GeneralParameters section must be declared in the configuration file. I found {tc.TFILE}{len(gp_input)}{tc.TRESET}"
        )
    return gp_inputs[0]


def get_general_parameters_from_configs(configs):
    gp_inputs = get_section_configs(configs, "GeneralParameters")
    gp_input = get_general_parameter_section(gp_inputs)
    gen_params = set_general_parameter_configs(gp_input)
    gen_params.parse_settings()
    return gen_params

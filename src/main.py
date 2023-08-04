import time

from general.argument_parser import analyse_parser_arguments
from general.argument_parser import create_parser
from general.report_time import print_runtime
from general.report_time import get_simulation_time_stamp
from general.configfile_tools import import_configurations
from classes.general_parameters import GeneralParameters
from general_parameters.general_parameters_funcs import get_general_parameters_from_configs


def run():
    simulation_time_stamp = get_simulation_time_stamp()

    parser = create_parser()
    user_args = parser.parse_args()
    user_inputs = analyse_parser_arguments(user_args, parser)
    configs = import_configurations(user_inputs.config_file_path)
    general_parameters = get_general_parameters_from_configs(configs)

    wall_clock_time_start = time.monotonic()
    cpu_time_start = time.process_time()

    time.sleep(2)

    wall_clock_time_end = time.monotonic()
    cpu_time_end = time.process_time()

    print_runtime(
        wall_clock_time_end - wall_clock_time_start,
        "Wall clock time of the simulation",
    )
    print_runtime(
        cpu_time_end - cpu_time_start,
        "CPU time of the simulation",
    )


if __name__ == "__main__":
    run()

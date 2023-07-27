import time

from source.general.argument_parser import analyse_parser_arguments
from source.general.argument_parser import create_parser





def main():
    wall_clock_time_start = time.monotonic()
    cpu_time_start = time.process_time()

    parser = create_parser()
    user_args = parser.parse_args()
    user_inputs = analyse_parser_arguments(user_args, parser)
    time.sleep(2)

    wall_clock_time_end = time.monotonic()
    cpu_time_end = time.process_time()

    print_runtime(
        wall_clock_time_end - wall_clock_time_start,
        "Wall clock time of the simulation:",
    )
    print_runtime(
        cpu_time_end - cpu_time_start,
        "CPU time of the simulation:",
    )


if __name__ == "__main__":
    main()

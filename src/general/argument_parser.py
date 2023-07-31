import argparse
import sys

from general.classes.UserInputs import UserInputs
from general.classes.PlatformInfo import print_available_platforms
from general.template_generator import config_file_template_generator
from general.resume_file_path_parser import find_resume_config_file


def create_parser():
    """Create a parser for VCM.

    Set parser with VCM's required user input arguments.

    Returns:
        UserInputs: Class containing parsed user inputs
    """
    parser = argparse.ArgumentParser(
        prog="VCM",
        description="Welcome to the Virtual Cell Model (VCM). VCM can simulate the mechanical behaviour of living cells and tissues using configuration files. If you already have a configuration file, please provide the path. To generate a configuration file, use the configuration generator option.",
    )
    parser.add_argument(
        "-g",
        "--generate-template",
        action="store_true",
        help="Generate a configuration file that contains all of VCM classes and their configurable options with default values and description.",
    )
    parser.add_argument(
        "-c",
        "--configfile-path",
        help="Path to the configuration file. example: Path/to/the/configuration/file",
        type=str,
        metavar="",
    )
    parser.add_argument(
        "-r",
        "--resume",
        type=str,
        metavar="",
        help="The directory path of the interrupted simulation that contains the simulation outputs. The simulation can only resume on the same machine and hardware it was originally running on to insure a consist run.",
    )
    parser.add_argument(
        "-a",
        "--available-platforms",
        action="store_true",
        help="View available platforms and print the flags for selected platform.",
    )
    parser.add_argument(
        "--platform",
        type=str,
        help='Platform (name or index) used for the simulation. Use "VCM -a, --available-platforms" to view the list of the available platforms.',
        metavar="",
    )
    parser.add_argument(
        "--platform-device-ID",
        type=int,
        help='ID of platform device to use for the simulation. Use "VCM -a, --available-platforms" to view the list of the available platform devices.',
        # default=0,
        metavar="",
    )
    # parser.add_argument(
    #     "--openmm-plugin-path",
    #     type=str,
    #     help="Path to OpenMM's plugins. Usually located at '/lib/plugins' in OpenMM's installation path",
    #     default="/usr/local/openmm/lib/plugins",
    #     metavar="",
    # )
    parser.add_argument(
        "--write-at-end",
        action="store_true",
        help="Write all outputs at the end of simulation. Warning: Resume will not be supported.",
    )
    return parser


def analyse_parser_arguments(user_args, parser):
    """Interpret arguments.

    Analyse and interpret user command line arguments and invoke appropriate action.

    Args:
        user_args (argparse.Namespace): User argument values.
        parser (argparse.ArgumentParser): Customised parser for VCM's required arguments.

    Returns:
        UserInputs: Inputs interpreted from user command line arguments.
    """
    user_inputs = UserInputs()
    # user_inputs.user_selected_platform = user_args.platform
    # user_inputs.user_selected_device = user_args.platform_device_ID
    # user_inputs.platform_info.openmm_plugin_path = user_args.openmm_plugin_path
    user_inputs.write_at_end = user_args.write_at_end
    if user_args.generate_template:
        config_file_template_generator()
        sys.exit(0)
    elif user_args.available_platforms:
        print_available_platforms(user_args.platform, user_args.platform_device_ID)
        sys.exit(0)
    elif user_args.resume is not None:
        user_inputs.resume_path = user_args.resume
        user_inputs.config_file_path = find_resume_config_file()
    elif user_args.config_file_path is not None:
        user_inputs.config_file_path = user_args.config_file_path
    else:
        parser.print_help()
        sys.exit(1)

    return user_inputs

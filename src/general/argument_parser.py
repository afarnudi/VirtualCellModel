import argparse
import sys

from general.classes.UserInputs import UserInputs
from general.platform_selection_tools import print_available_platforms_and_devices
from general.template_generator import config_file_template_generator
from general.configfile_tools import find_resume_config_file
from general.platform_selection_tools import get_platform_index_and_name
from general.platform_selection_tools import get_platform_device
from general.configfile_tools import check_file_path
from general.configfile_tools import parse_dir_path


def create_parser():
    """Create a parser for VCM.

    Set parser with VCM's required user input arguments.

    Returns:
        UserInputs: Class containing parsed user inputs
    """
    parser = argparse.ArgumentParser(
        prog="VCM",
        description="Welcome to the Virtual Cell Model (VCM). VCM can simulate the mechanical behaviour of living cells and tissues using configuration files. If you already have a configuration file, please provide the path. To generate a configuration file, use the configuration template generator option.",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-g",
        "--generate-template",
        action="store_true",
        help="Generate a configuration file that contains all of VCM classes and their configurable options with default values and description.",
    )
    group.add_argument(
        "-a",
        "--available-platforms",
        action="store_true",
        help="View available platforms and print the flags for selected platform.",
    )
    group.add_argument(
        "-c",
        "--configfile-path",
        help="Path to the configuration file. example: Path/to/the/configuration/file",
        type=str,
        metavar="",
        default=None,
    )
    group.add_argument(
        "-r",
        "--resume",
        type=str,
        metavar="",
        help="The directory path of the interrupted simulation that contains the simulation outputs. The simulation can only resume on the same machine and hardware it was originally running on to insure a consist run.",
    )
    parser.add_argument(
        "--platform",
        type=str,
        help='Platform (name or index) used for the simulation. Use "VCM -a, --available-platforms" to view the list of the available platforms.',
        metavar="",
        default=None,
    )
    parser.add_argument(
        "--platform-device-ID",
        type=int,
        help='ID of platform device to use for the simulation. Use "VCM -a, --available-platforms" to view the list of the available platform devices.',
        # default=0,
        metavar="",
        default=None,
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
    user_inputs.write_at_end = user_args.write_at_end
    if user_args.generate_template:
        config_file_template_generator()
        sys.exit(0)
    if user_args.available_platforms:
        print_available_platforms_and_devices(
            user_args.platform, user_args.platform_device_ID
        )
        sys.exit(0)
    if user_args.resume is not None:
        user_inputs.resume_path = parse_dir_path(user_args.resume)
        user_inputs.config_file_path = find_resume_config_file()
    if user_args.configfile_path is not None:
        user_inputs.config_file_path = check_file_path(user_args.configfile_path)
    
    user_inputs.platform_info.index, user_inputs.platform_info.name = get_platform_index_and_name(user_args.platform)
    if user_inputs.platform_info.name != "Reference":
        user_inputs.platform_info.device_ID = get_platform_device(user_inputs.platform_info.name, user_args.platform_device_ID)

    return user_inputs

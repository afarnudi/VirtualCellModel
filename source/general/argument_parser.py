import argparse
import sys

from source.general.classes.UserInputs import UserInputs
from source.general.classes.PlatformInfo import print_platform_info
from source.general.template_generator import config_file_template_generator
from source.general.resume_file_path_parser import find_resume_config_file


def create_parser():
    """
    Create a parser

    Set parser with VCM's required user input arguments.

    Returns:
        UserInputs: Class containing parsed user inputs
    """
    parser = argparse.ArgumentParser(
        prog="VCM",
        description="The Virtual Cell Model (VCM)\nVCM can simulate the mechanical behaviour of living cells and tissues. To generate a configuration file, use the configuration generator. If you already have a configuration file, please provide teh path.",
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
        help='Platform used for the simulation. Use "VCM -a, --available-platforms" to view the list of the available platforms.',
        metavar="",
    )
    parser.add_argument(
        "--platform-device-ID",
        type=int,
        help='ID of platform device to use for the simulation. Use "VCM -a, --available-platforms" to view the list of the available platform devices.',
        metavar="",
    )
    parser.add_argument(
        "--openmm-plugin-path",
        type=str,
        help="Path to OpenMM's plugins. Usually located at '/lib/plugins' in OpenMM's installation path",
        default="/usr/local/openmm/lib/plugins",
        metavar="",
    )
    parser.add_argument(
        "--write-at-end",
        action="store_true",
        help="Write all outputs at the end of simulation. Warning: Resume will not be supported.",
    )
    return parser


def analyse_parser_argumetns(user_args, parser):
    """
    Interpret arguments

    Analyse and interpret user command line arguments and invoke appropriate action.

    Args:
        user_args : argparse.Namespace
            user argument values
        parser : argparse.ArgumentParser
            customised parser for VCM's required arguments

    Returns:
        _type_: _description_
    """
    user_inputs = UserInputs()
    if user_args.generate_template:
        config_file_template_generator()
        sys.exit(0)
    elif user_args.available_platforms:
        print_platform_info(user_inputs)
        sys.exit(0)
    elif user_args.resume is not None:
        user_inputs.resume_path = user_args.resume
        user_inputs.config_file_path = find_resume_config_file()
    elif user_args.configfile_path is not None:
        user_inputs.config_file_path = user_args.configfile_path
    else:
        parser.print_help()
        sys.exit(1)

    user_inputs.platform_info.platform_name = user_args.platform
    user_inputs.platform_info.platform_device_ID = user_args.platform_device_ID

    return user_inputs

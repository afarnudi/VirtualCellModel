import argparse
import sys

from classes.general.UserInputs import UserInputs
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
        description="Generate a configuration file that contains all of VCM classes and their configurable options with default values and description.",
    )
    subparsers = parser.add_subparsers(required=True, metavar="", dest="command")

    # Template parser def
    parser_generate_template = subparsers.add_parser(
        "generate-template",
        help="Preprocess the labelled images for shape deformation analysis.",
    )

    # Available platforms parser def
    parser_available_platforms = subparsers.add_parser(
        "available-platforms",
        help="View available platforms and print the flags for selected platform.",
    )

    # Run parser def
    parser_run = subparsers.add_parser(
        "run",
        help="Run simulations.",
    )

    # Resume parser def
    parser_resume = subparsers.add_parser(
        "resume",
        help="Resume simulations.",
    )

    for p in [parser_available_platforms, parser_run, parser_resume]:
        p.add_argument(
            "-p",
            "--platform",
            type=str,
            help='Platform (name or index) used for the simulation. Use "VCM -a, --available-platforms" to view the list of the available platforms.',
            metavar="",
            default=None,
        )
        p.add_argument(
            "-d",
            "--platform-device-ID",
            type=int,
            help='ID of platform device to use for the simulation. Use "VCM -a, --available-platforms" to view the list of the available platform devices.',
            # default=0,
            metavar="",
            default=None,
        )
        p.add_argument(
            "--openmm-plugin-path",
            type=str,
            help="Path to OpenMM's plugins. Usually located at '/lib/plugins' in OpenMM's installation path",
            default="/usr/local/openmm/lib/plugins",
            metavar="",
        )

    for p in [parser_run, parser_resume]:
        p.add_argument(
            "configfile_path",
            help="Path to the configuration file. example: Path/to/the/configuration/file",
            type=str,
            default=None,
        )
        p.add_argument(
            "--write-at-end",
            action="store_true",
            help="Write all outputs at the end of simulation. Warning: Resume will not be supported.",
        )
    return parser


def analyse_parser_arguments(user_args):
    """Interpret arguments.

    Analyse and interpret user command line arguments and invoke appropriate action.

    Args:
        user_args (argparse.Namespace): User argument values.

    Returns:
        UserInputs: Inputs interpreted from user command line arguments.
    """
    user_inputs = UserInputs()

    if user_args.command == "generate-template":
        config_file_template_generator()
        sys.exit(0)
    if user_args.command == "available-platforms":
        print_available_platforms_and_devices(
            user_args.platform, user_args.platform_device_ID
        )
        sys.exit(0)

    user_inputs.write_at_end = user_args.write_at_end
    if user_args.command == "resume":
        user_inputs.resume_path = parse_dir_path(user_args.configfile_path)
        user_inputs.config_file_path = find_resume_config_file()
    if user_args.command == "run":
        user_inputs.config_file_path = check_file_path(user_args.configfile_path)

    (
        user_inputs.platform_info.index,
        user_inputs.platform_info.name,
    ) = get_platform_index_and_name(user_args.platform)
    if user_inputs.platform_info.name != "Reference":
        user_inputs.platform_info.device_ID = get_platform_device(
            user_inputs.platform_info.name, user_args.platform_device_ID
        )

    return user_inputs

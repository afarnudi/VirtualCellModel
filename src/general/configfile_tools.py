import os
import sys
from classes.general.configuration import Configuration
from classes.general.configuration import is_section_declaration
from classes.general.configuration import SECTION_LIST
from classes.general.print_colors import TerminalColors as tc


def find_resume_config_file():
    print("find_resume_config_file is under development")
    sys.exit()


def check_file_path(file_path):
    """Check if file path exists.

    Args:
        file_path (str): The path to the file.

    Raises:
        FileNotFoundError: If the file was not found.

    Returns:
        str: Path to the file.
    """
    if os.path.isfile(file_path) is False:
        raise FileNotFoundError(
            f'{tc.TRESET}No such file or directory: {tc.TFILE}"{file_path}"{tc.TRESET}.'
        )
    return file_path


def parse_dir_path(dir_path):
    """Check if directory exists.

    Args:
        dir_path (str): Directory path.

    Raises:
        FileNotFoundError: If directory does not exist._

    Returns:
        str: Directory path guarantied to end with "/". f dir_path is a file path, returns directory containing the file.
    """
    if os.path.isdir(dir_path):
        path = dir_path
    else:
        if os.path.isfile(dir_path):
            path = os.path.dirname(dir_path)
        else:
            raise FileNotFoundError(
                f'{tc.TRESET}No such file or directory: {tc.TFILE}"{dir_path}"{tc.TRESET}.'
            )
    if path.endswith("/") is False:
        path += "/"
    return path


def read_config_file_lines(config_file_path):
    """Open and read the lines of a text based configuration file.

    Args:
        config_file_path (str): Configuration file path.

    Returns:
        list: list of strings containing the configuration file lines.
    """
    with open(config_file_path, "r") as r:
        lines = r.readlines()
    return lines


def get_configurations(conf_lines):
    """Create VCM configuration objects from the configuration file content.

    Check configuration lines for section decelerations, create VCM configuration objects by adding the user input settings.

    Args:
        conf_lines (list[str]): List of strings containing the configuration file cleaned up lines.

    Raises:
        TypeError: If the first line of the configuration does not contain a section declaration.

    Returns:
        list: List of Configuration objects with the section configurations.
    """

    configurations = []
    for line in conf_lines:
        if is_section_declaration(line):
            configurations.append(Configuration(line))
        elif configurations:
            configurations[-1].add(line)
        else:
            continue
    return configurations


def configurations_consistency_check(configurations):
    """Consistency check for user configuration inputs.

    Args:
        configurations (list[Configuration]): List of user selected VCM sections.

    Raises:
        TypeError: When there is no GeneralParameters sections or it is not uniquely defined.
        TypeError: When there are no sections aside from a GeneralParameters and InteractionTable section.
        TypeError: When there is more than one InteractionTable section.


    Returns:
        bool: True if all checks pass.
    """

    section_name_classes = SECTION_LIST.copy()
    section_name_classes.remove("InteractionTable")
    section_name_classes.remove("GeneralParameters")

    section_names = [c.get_name() for c in configurations]

    if section_names.count("GeneralParameters") != 1:
        raise TypeError(
            f"{tc.TFAILED} The configuration file must contain {tc.TBOLD}one{tc.TRESET}{tc.TFAILED} GeneralParameters section.{tc.TRESET}"
        )
    section_names.remove("GeneralParameters")

    interaction_table_count = section_names.count("InteractionTable")
    if interaction_table_count == 0:
        print(
            f"{tc.TBLINK}{tc.TWARN}No interaction section in the configuration file. Non bonded interactions will be switched off.{tc.TRESET}"
        )
    elif interaction_table_count > 1:
        raise TypeError(
            f"{tc.TFAILED} The configuration file must contain {tc.TBOLD}at most one{tc.TRESET}{tc.TFAILED} InteractionTable section.{tc.TRESET}"
        )
    else:
        section_names.remove("InteractionTable")

    for sec in section_names:
        if sec not in section_name_classes:
            raise TypeError(
                f'{tc.TRESET}Apart from the "GeneralParameters" section, the configuration file must also contain at least one of the {section_name_classes} sections.'
            )

    return True


def import_configurations(config_file_path):
    """Import user configurations from file.

    Args:
        config_file_path (str): Path to the configuration file.

    Returns:
        list[Configuration]: List of Configuration objects containing user sections and settings.
    """
    conf_lines = read_config_file_lines(config_file_path)
    configurations = get_configurations(conf_lines)
    configurations_consistency_check(configurations)
    return configurations

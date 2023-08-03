import os
import sys
from classes.general.configuration import Configuration
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


def strip_of_comments(line, comment="#"):
    """Remove characters to the right of the comment symbol (including the comment symbol).

    Args:
        line (str): String containing statements and/or comments.
        comment (str, optional): String identified as the comment symbol. Defaults to "#".

    Returns:
        str: String content that lies on the left of the comment symbol.
    """
    return line.split(comment)[0].rstrip().lstrip()


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


def clean_lines(conf_lines):
    """Preprocess and clean up the configuration lines lines.

    During preprocessing remove comments, empty lines, and trailing white space on the right and left of the line.

    Args:
        conf_lines (list): list of configuration lines (strings).

    Returns:
        list: list of lines in the form of strings.
    """
    clean_lines = []
    for line in conf_lines:
        line = strip_of_comments(line)
        line = line.rstrip()
        line = line.lstrip()
        if len(line) == 0:
            continue
        clean_lines.append(line)
    return clean_lines


def get_configurations(conf_lines):
    """Extract sections and associated configuration lines from the list of configurations.

    Args:
        conf_lines (list): List of strings containing the configuration file cleaned up lines.

    Raises:
        TypeError: If the first line of the configuration does not contain a section declaration.

    Returns:
        list: List of Configuration objects with the section configurations.
    """
    first_line = conf_lines[0]
    configurations = []
    if Configuration.check_string_for_section_declaration(first_line):
        configurations.append(Configuration(first_line))
        active_conf = configurations[-1]
    else:
        raise TypeError(
            f'{tc.TRESET}{tc.TFAILED}Configuration file parsing: "{tc.TFILE}{first_line}{tc.TFAILED}" is not a section name. First non-comment statement of a configuration file must begin with a section declaration.\n{Configuration.section_declaration_message}{tc.TRESET}'
        )
    for line in conf_lines[1:]:
        if Configuration.check_string_for_section_declaration(line):
            configurations.append(Configuration(line))
            active_conf = configurations[-1]
            continue
        active_conf.add(line)
    return configurations

def check_configurations(configurations):
    section_name_classes = Configuration.SECTION_LIST.copy()
    section_name_classes.remove("InteractionTable")
    section_name_classes.remove("GeneralParameters")

    section_names = [c.get_name() for c in configurations]
    if "GeneralParameters" in section_names:
        section_names.remove("GeneralParameters")
        try:
            section_names.remove("InteractionTable")
        except ValueError:
            print(f"{tc.TBLINK}{tc.TWARN}No interaction section in the configuration file. Non bonded interactions will be switched off.{tc.TRESET}")
        for sec in section_names:
            if sec not in section_name_classes:
                raise TypeError(f"{tc.TRESET}Apart from the \"GeneralParameters\" section, the configuration file must also contain at least one of the {section_name_classes} sections.")
    else:
        raise TypeError(f"{tc.TFAILED} The configuration file must contain a GeneralParameters section.")
    return True


def import_configurations(config_file_path):
    conf_lines = read_config_file_lines(config_file_path)
    conf_lines = clean_lines(conf_lines)
    configurations = get_configurations(conf_lines)
    check_configurations(configurations)
    return configurations

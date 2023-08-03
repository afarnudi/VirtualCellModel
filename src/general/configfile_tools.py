import os
import sys
from general.classes.configuration import Configuration
from general.classes.print_colors import TerminalColors as tc


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
    return line.split(comment)[0].rstrip().lstrip()


def read_config_file_lines(config_file_path):
    with open(config_file_path, "r") as r:
        lines = r.readlines()
    return lines


def clean_lines(conf_lines):
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
    first_line = conf_lines[0]
    configurations = []
    if Configuration.check_string_for_class_declaration(first_line):
        configurations.append(Configuration(first_line))
        active_conf = configurations[-1]
    else:
        raise TypeError(
            f'{tc.TRESET}{tc.TFAILED}Configuration file parsing: "{tc.TFILE}{first_line}{tc.TFAILED}" is not a class name. First non-comment statement of a configuration file must begin with a class declaration.\nClass declarations must be at the beginning of the line and in the following format "-ClassName". ClassName choices are: {Configuration.CLASS_LIST}{tc.TRESET}'
        )
    for line in conf_lines[1:]:
        if Configuration.check_string_for_class_declaration(line):
            configurations.append(Configuration(line))
            active_conf = configurations[-1]
            continue
        active_conf.add(line)
    return configurations


def import_configurations(config_file_path):
    conf_lines = read_config_file_lines(config_file_path)
    conf_lines = clean_lines(conf_lines)
    configurations = get_configurations(conf_lines)
    return configurations

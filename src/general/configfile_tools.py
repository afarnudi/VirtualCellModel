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


def get_configurations(config_file_path):
    conf_lines = read_config_file_lines(config_file_path)
    configurations = []
    for line in conf_lines:
        line = strip_of_comments(line)
        if len(line) == 0:
            continue
        if Configuration.check_string_for_class_declaration(line):
            configurations.append(Configuration(line))
            active_conf = configurations[-1]
            continue
        active_conf.add(line)
    return configurations
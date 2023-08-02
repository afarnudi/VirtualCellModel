import os
import sys
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

import pytest
from src.general.configfile_tools import check_file_path
from src.general.configfile_tools import parse_dir_path
from src.general.configfile_tools import strip_of_comments


def test_check_file_path_file():
    file_path = "tests/test_general/test_configfile_tools.py"
    assert check_file_path(file_path) == file_path


def test_check_file_path_file():
    file_path = "tests/test_general/NO_SUCH_FILE"
    with pytest.raises(FileNotFoundError) as exc_info:
        check_file_path(file_path)
        assert "No such file or directory: " in str(exc_info)


def test_check_file_path_dir():
    file_path = "tests/test_general/"
    with pytest.raises(FileNotFoundError) as exc_info:
        check_file_path(file_path)
        assert "No such file or directory: " in str(exc_info)


def test_parse_dir_path_file():
    dir_path = "tests/test_general/test_configfile_tools.py"
    answer = "tests/test_general/"
    assert parse_dir_path(dir_path) == answer


def test_parse_dir_path_non_existing_dir():
    dir_path = "tests/test_general/NO_SUCH_DIR"
    with pytest.raises(FileNotFoundError) as exc_info:
        parse_dir_path(dir_path)
        assert "No such file or directory: " in str(exc_info)


def test_parse_dir_path_empty_dir():
    dir_path = ""
    with pytest.raises(FileNotFoundError) as exc_info:
        parse_dir_path(dir_path)
        assert "No such file or directory: " in str(exc_info)


def test_parse_dir_current_dir():
    dir_path = "."
    answer = "./"
    assert parse_dir_path(dir_path) == answer


def test_parse_dir_upper_dir():
    dir_path = ".."
    answer = "../"
    assert parse_dir_path(dir_path) == answer


def test_parse_dir_from_upper_dir():
    dir_path = "src/../tests/"
    answer = "src/../tests/"
    print(parse_dir_path(dir_path))
    assert parse_dir_path(dir_path) == answer


def test_strip_of_comments_empty():
    assert strip_of_comments("") == ""


def test_strip_of_comments_comment():
    assert strip_of_comments("#") == ""


def test_strip_of_comments_long_comment():
    assert strip_of_comments("# some comment") == ""


def test_strip_of_comments_no_comment():
    line = "statement"
    assert strip_of_comments(line) == line


def test_strip_of_comments_no_comment_long():
    line = "statement about something."
    assert strip_of_comments(line) == line


def test_strip_of_comments_no_comment_long():
    line = "statement about something."
    assert strip_of_comments(line) == line


def test_strip_of_comments_statement_and_comment():
    line = "statement about something."
    comment = "# some comment"
    assert strip_of_comments(line + comment) == line

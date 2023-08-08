import pytest
from src.classes.general.configuration import Configuration
from src.classes.general.configuration import clean_line
from src.classes.general.configuration import strip_of_comments
from src.classes.general.configuration import is_section_declaration


def test_configuration_init_correct_name():
    for name in Configuration.SECTION_LIST:
        temp = Configuration(f"-{name}")
        assert temp.get_name() == name


def test_configuration_init_incorrect_name():
    for name in Configuration.SECTION_LIST:
        with pytest.raises(TypeError) as exc_info:
            temp = Configuration(name)
            assert (
                'Section declarations must be at the beginning of the line and in the following format "-SectionName". SectionName choices are: '
                in str(exc_info)
            )


def test_is_section_declaration_correct_name():
    for name in Configuration.SECTION_LIST:
        assert is_section_declaration(f"-{name}") == True


def test_is_section_declaration_normal_line():
    for name in Configuration.SECTION_LIST:
        assert is_section_declaration(name) == False


def test_is_section_declaration_incorrect_name():
    for name in Configuration.SECTION_LIST:
        with pytest.raises(TypeError) as exc_info:
            is_section_declaration(f"-{name}s")
            assert (
                'Section declarations must be at the beginning of the line and in the following format "-SectionName". SectionName choices are: '
                in str(exc_info)
            )

def test_all_conf_values_are_strings():
    name = Configuration.SECTION_LIST[0]
    temp = Configuration(f"-{name}")
    temp.add("conf1 0 0 0")
    temp.add("conf2 set1 set2")
    temp.add("conf3 Au")
    temp.add("conf4 true")
    temp.add("conf5 Au true")
    for key in temp.get_key_names():
        assert type(temp.get_value(key)) == str



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


def test_clean_line():
    lines = [
        "#some comments",
        "   #comments",
        "",
        "-Section Name",
        "    -Section Name",
        "-Section Name     ",
        "    -Section Name    ",
        "    -Section Name   #asdfasdfadsf ",
        "-Section Name   #asdfasdfadsf ",
        "-Section Name   stuff here #asdfasdfadsf ",
        "   -Section Name   stuff here #asdfasdfadsf ",
        "-Section Name   stuff here ",
        "-Section Name   stuff here",
        "   -Section Name   stuff here",
        "configuration",
        "   configuration   ",
        "   configuration",
        "configuration     ",
        "configuration    2134 #comments here ",
        "   configuration    2134 234 #comments here ",
        "   configuration    Au #comments here ",
        "   configuration    1,2,3,4 #comments here ",
        "   configuration     #comments here ",
        "   configuration     #comments here ",
    ]
    answers = [
        "",
        "",
        "",
        "-Section Name",
        "-Section Name",
        "-Section Name",
        "-Section Name",
        "-Section Name",
        "-Section Name",
        "-Section Name   stuff here",
        "-Section Name   stuff here",
        "-Section Name   stuff here",
        "-Section Name   stuff here",
        "-Section Name   stuff here",
        "configuration",
        "configuration",
        "configuration",
        "configuration",
        "configuration    2134",
        "configuration    2134 234",
        "configuration    Au",
        "configuration    1,2,3,4",
        "configuration",
        "configuration",
    ]
    for line, answer in zip(lines, answers):
        cleaned_line = clean_line(line)
        assert answer == cleaned_line
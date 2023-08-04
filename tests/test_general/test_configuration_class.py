import pytest
from src.classes.general.configuration import Configuration


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


def test_check_string_for_section_declaration_correct_name():
    for name in Configuration.SECTION_LIST:
        assert Configuration.check_string_for_section_declaration(f"-{name}") == True


def test_check_string_for_section_declaration_normal_line():
    for name in Configuration.SECTION_LIST:
        assert Configuration.check_string_for_section_declaration(name) == False


def test_check_string_for_section_declaration_incorrect_name():
    for name in Configuration.SECTION_LIST:
        with pytest.raises(TypeError) as exc_info:
            Configuration.check_string_for_section_declaration(f"-{name}s")
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
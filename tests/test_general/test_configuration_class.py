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


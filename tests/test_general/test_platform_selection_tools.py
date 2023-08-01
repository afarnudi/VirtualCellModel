from unittest import mock
import pytest
from src.general.platform_selection_tools import get_platform_device_from_user
from src.general.platform_selection_tools import get_platform_device_properties
from src.general.platform_selection_tools import parse_platform_device_selection
from src.general.platform_selection_tools import check_index_type
from openmm import OpenMMException


def test_get_platform_device_from_user():
    for val in ["0", "1", "2", "3", "-1", "4", "5", "6548", "A", "#", "ali", "Akbar"]:
        with mock.patch("builtins.input", return_value=f"{val}"):
            assert get_platform_device_from_user() == val


def test_check_index_type_string_int():
    for val in ["-1", "4", "5", "6548"]:
        with mock.patch("builtins.input", return_value=f"{val}"):
            assert check_index_type(val) == int(val)


def test_check_index_type_string_non_int():
    for val in ["3.4", "Ali", "#", "bad device index", "ALLCAP", "0.4", "-1."]:
        with pytest.raises(ValueError) as exc_info, mock.patch(
            "builtins.input", return_value=f"{val}"
        ):
            check_index_type(val)
            assert " not an integer index." in str(exc_info.value)


def test_parse_platform_device_selection_correct_int():
    for val in ['0', '1', '2', '3', '4']:
        result = parse_platform_device_selection([0, 1, 2, 3, 4], val)
        assert result == int(val)


def test_parse_platform_device_selection_incorrect_int():
    for val in ['-1', '5', '658', '7', '13']:
        with pytest.raises(ValueError) as exc_info:
            parse_platform_device_selection([0, 1, 2, 3, 4], val)
            assert "not in list of available device indices:" in str(exc_info)

def test_parse_platform_device_selection_incorrect_int():
    for val in ["3.4", "Ali", "#", "bad device index", "ALLCAP", "0.4", "-1."]:
        with pytest.raises(ValueError) as exc_info:
            parse_platform_device_selection([0, 1, 2, 3, 4], val)
            assert " not an integer index." in str(exc_info.value)


def test_parse_platform_device_selection_charachter():
    for val in ["a", "A", "!", "$", "# *", "*", "-"]:
        with pytest.raises(ValueError) as exc_info, mock.patch(
            "builtins.input", return_value="Ali"
        ):
            parse_platform_device_selection([0, 1, 2, 3], val)
            assert " not an integer index." in str(exc_info.value)


def test_get_platform_device_properties_Reference():
    result = get_platform_device_properties("Reference")
    assert len(result) == 0


def test_get_platform_device_properties_CPU():
    result = get_platform_device_properties("CPU")
    assert len(result) != 0


# def test_get_platform_device_properties_wrong_platform_name():
#     with pytest.raises(OpenMMException) as exc_info, mock.patch(
#             "builtins.input", return_value=f"{val}"
#         ):
#         result = get_platform_device_properties("Bad platform name")
#     assert "There is no registered Platform called" in str(exc_info)




from unittest import mock
import pytest
from src.general.platform_selection_tools import get_platform_device_from_user
from src.general.platform_selection_tools import get_platform_device_properties
from src.general.platform_selection_tools import parse_platform_device_selection
from openmm import OpenMMException


def test_get_platform_device_from_user_correct_int():
    for val in [0, 1, 2, 3]:
        with mock.patch("builtins.input", return_value=f"{val}"):
            assert get_platform_device_from_user([0, 1, 2, 3]) == val


def test_get_platform_device_from_user_incorrect_int():
    for val in [-1, 4, 5, 6548]:
        with mock.patch("builtins.input", return_value=f"{val}"):
            assert get_platform_device_from_user([0, 1, 2, 3]) == val


def test_get_platform_device_from_user_float():
    for val in [-1.0, 0.4, 5.2, 65.48]:
        with pytest.raises(ValueError) as exc_info, mock.patch(
            "builtins.input", return_value=f"{val}"
        ):
            get_platform_device_from_user([0, 1, 2, 3])
            assert " not an integer index." in str(exc_info.value)


def test_get_platform_device_from_user_string():
    for val in ["ali", "ALI", "Elham", "DAvid", "John Doe", "Jane_Doe", "Akbar"]:
        with pytest.raises(ValueError) as exc_info, mock.patch(
            "builtins.input", return_value="Ali"
        ):
            get_platform_device_from_user([0, 1, 2, 3])
            assert " not an integer index." in str(exc_info.value)


def test_get_platform_device_from_user_charachter():
    for val in ["a", "A", "!", "$", "# *", "*", "-"]:
        with pytest.raises(ValueError) as exc_info, mock.patch(
            "builtins.input", return_value="Ali"
        ):
            get_platform_device_from_user([0, 1, 2, 3])
            assert " not an integer index." in str(exc_info.value)


def test_get_platform_device_properties_Reference():
    result = get_platform_device_properties("Reference")
    assert len(result) == 0


def test_get_platform_device_properties_CPU():
    result = get_platform_device_properties("CPU")
    assert len(result) != 0


def test_get_platform_device_properties_wrong_platform_name():
    with pytest.raises(OpenMMException) as exc_info:
        result = get_platform_device_properties("Bad platform name")
    assert "There is no registered Platform called" in str(exc_info)


def test_parse_platform_device_selection_correct_int():
    for val in [0, 1, 2, 3, 4]:
        result = parse_platform_device_selection([0, 1, 2, 3, 4], val)
        assert result == val


def test_parse_platform_device_selection_incorrect_int():
    for val in [-1, 5, 658, 7, 13]:
        with pytest.raises(ValueError) as exc_info:
            parse_platform_device_selection([0, 1, 2, 3, 4], val)
            assert "not in list of available device indices:" in str(exc_info)
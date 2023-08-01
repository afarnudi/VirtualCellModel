from unittest import mock
import pytest
from src.general.platform_selection_tools import get_platform_device_from_user
from src.general.platform_selection_tools import get_platform_device_properties
from src.general.platform_selection_tools import parse_platform_device_selection
from src.general.platform_selection_tools import check_index_type
from src.general.platform_selection_tools import get_list_of_platform_names
from src.general.platform_selection_tools import parse_platform_index
from src.general.platform_selection_tools import parse_platform_name
from src.general.platform_selection_tools import parse_selected_platform
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
    for val in ["0", "1", "2", "3", "4"]:
        result = parse_platform_device_selection([0, 1, 2, 3, 4], val)
        assert result == int(val)


def test_parse_platform_device_selection_incorrect_int():
    for val in ["-1", "5", "658", "7", "13"]:
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


def test_get_list_of_platform_names():
    name_list = get_list_of_platform_names()
    for name in ["Reference", "CPU"]:
        assert name in name_list


def test_parse_platform_index_non_int():
    for val in ["0.", "ali", "platform", "%", "#", "0.5", "-10.2", "38.15"]:
        assert (None, None) == parse_platform_index(val)


def test_parse_platform_index_Reference():
    assert (0, "Reference") == parse_platform_index("0")


def test_parse_platform_index_CPU():
    assert (1, "CPU") == parse_platform_index("1")


def test_parse_platform_index_non_platform():
    with pytest.raises(ValueError) as exc_info:
        parse_platform_index("10")
        assert " not in list of available platforms indices. " in str(exc_info.value)


def test_parse_platform_name_non_platform_name():
    for val in [
        "cpu",
        "reference",
        "opencl",
        "cuda",
        "Akbar",
        "ali",
        "platform",
        "%",
        "#",
        "0.5",
        "-10.2",
        "38.15",
        "0",
        "1",
    ]:
        with pytest.raises(ValueError) as exc_info:
            parse_platform_name(val)
            assert " not in list of available platforms: " in str(exc_info.value)


def test_parse_platform_name_Reference():
    assert (0, "Reference") == parse_platform_name("Reference")


def test_parse_platform_name_CPU():
    assert (1, "CPU") == parse_platform_name("CPU")


def test_parse_selected_platform_Reference():
    for val in ["0", "Reference"]:
        assert (0, "Reference") == parse_selected_platform(val)


def test_parse_selected_platform_CPU():
    for val in ["1", "CPU"]:
        assert (1, "CPU") == parse_selected_platform(val)

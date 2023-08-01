from unittest import mock
import pytest
from src.general.platform_selection_tools import get_platform_device_from_user


def test_get_platform_device_from_user_correct_int():
    for val in [0,1,2,3]:
        with mock.patch('builtins.input', return_value=f"{val}"):
            assert get_platform_device_from_user([0,1,2,3]) == val


def test_get_platform_device_from_user_incorrect_int():
    for val in [-1,4,5,6548]:
        with mock.patch('builtins.input', return_value=f"{val}"):
            assert get_platform_device_from_user([0,1,2,3]) == val


def test_get_platform_device_from_user_float():
    for val in [-1.,.4,5.2,65.48]:
        with pytest.raises(ValueError)  as exc_info, mock.patch('builtins.input', return_value=f"{val}"):
            get_platform_device_from_user([0,1,2,3])
            assert " not an integer index." in str(exc_info.value)

def test_get_platform_device_from_user_string():
    for val in ['ali', 'ALI', 'Elham', 'DAvid', 'John Doe', 'Jane_Doe', "Akbar"]:
        with pytest.raises(ValueError)  as exc_info, mock.patch('builtins.input', return_value="Ali"):
            get_platform_device_from_user([0,1,2,3])
            assert " not an integer index." in str(exc_info.value)

def test_get_platform_device_from_user_charachter():  
    for val in ['a', 'A', '!', '$', '# *', '*', "-"]:  
        with pytest.raises(ValueError)  as exc_info, mock.patch('builtins.input', return_value="Ali"):
            get_platform_device_from_user([0,1,2,3])
            assert " not an integer index." in str(exc_info.value)

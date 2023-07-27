import pytest
from src.general.argument_parser import create_parser
from src.general.argument_parser import analyse_parser_arguments


def test_analyse_parser_template_generate_short_flag_register():
    parser = create_parser()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        user_args = parser.parse_args(["-g"])
        user_inputs = analyse_parser_arguments(user_args, parser)
    # assert user_inputs.generate_template == True
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_analyse_parser_template_generate_long_flag_register():
    parser = create_parser()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        user_args = parser.parse_args(["--generate-template"])
        user_inputs = analyse_parser_arguments(user_args, parser)
    # assert user_inputs.generate_template == True
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0

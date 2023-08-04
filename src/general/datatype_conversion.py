from classes.general.print_colors import TerminalColors as tc


def string_to_float(string, caller_title, help=None):
    try:
        result = float(string)
    except:
        raise TypeError(
            f'{tc.TFAILED}{caller_title} Was not able to convert "{tc.TFILE}{string}{tc.TFAILED}" to a float. Please consult the setting description:\n{tc.TRESET}{help}'
        )
    return result


def string_to_int(string, caller_title, help=None):
    try:
        result = int(string)
    except:
        raise TypeError(
            f'{tc.TFAILED}{caller_title} Was not able to convert "{tc.TFILE}{string}{tc.TFAILED}" to an integer. Please consult the setting description:\n{tc.TRESET}{help}'
        )
    return result


def string_to_bool(string, caller_title, help=None):
    if string.lower() in ["true", "1"]:
        return True
    elif string.lower() in ["false", "0"]:
        return False
    else:
        raise TypeError(
            f'{tc.TFAILED}{caller_title} Could not parse "{tc.TFILE}{string}{tc.TFAILED}". Please use "true" or "false". Please consult the setting description:\n{tc.TRESET}{help}'
        )


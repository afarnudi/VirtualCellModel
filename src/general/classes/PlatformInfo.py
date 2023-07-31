import sys
from openmm.openmm import Platform
from general.classes.print_colors import TerminalColors as tc


class PlatformInfo:
    def __init__(self):
        self.name = None
        self.index = None
        self.device_ID = 0
        self.openmm_plugin_path = None


def get_list_of_platforms():
    num_platforms = Platform.getNumPlatforms()
    platform_names = []
    for i in range(num_platforms):
        platform = Platform.getPlatform(i)
        platform_names.append(platform.getName().lower())
    return platform_names


def parse_selected_platform(selected_platform):
    selected_platform = selected_platform.lower()
    platform_names = get_list_of_platforms()
    num_platforms = len(platform_names)

    if selected_platform in platform_names:
        platform_name = selected_platform
        platform_index = platform_names.index(selected_platform)
    else:
        try:
            user_platform_id = int(selected_platform)
        except:
            print(
                f'{tc.TRESET}"{tc.TFILE}{selected_platform}{tc.TRESET}" not in list of available platforms: {platform_names}'
            )
            raise SystemExit
        if user_platform_id in range(num_platforms):
            platform = Platform.getPlatform(user_platform_id)
            platform_name = platform.getName()
            platform_index = user_platform_id
        else:
            print(
                f'{tc.TRESET}"{tc.TFILE}{user_platform_id}{tc.TRESET}" not in list of available platforms indices. {list(range(num_platforms))}'
            )
            raise SystemExit

    print(tc.TRESET, end="")
    return platform_index, platform_name


def get_platform_from_user():
    platform_names = get_list_of_platforms()

    print(
        f"{tc.TOMM}\nOpenMM available platforms:\n{tc.TGRAY}Index Name \t  Speed (Estimated){tc.TRESET}"
    )
    for index, name in enumerate(platform_names):
        platform = Platform.getPlatform(index)
        print(f"({tc.TBOLD}{index}{tc.TRESET})  {name}\t   {platform.getSpeed():0.0f}")
    user_platform_input = input(f"Please choose a platform: \n{tc.TFILE}")
    print(tc.TRESET, end="")
    platform_index, platform_name = parse_selected_platform(user_platform_input)
    return platform_index, platform_name


def get_platform_index_and_name(selected_platform):
    if selected_platform is None:
        platform_index, platform_name = get_platform_from_user()
    else:
        platform_index, platform_name = parse_selected_platform(selected_platform)
    return platform_index, platform_name


def print_available_platforms(user_inputs):
    platform_index, platform_name = get_platform_index_and_name(user_inputs.user_selected_platform)
    
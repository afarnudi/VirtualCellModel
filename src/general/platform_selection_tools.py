from openmm.openmm import Platform
from general.classes.print_colors import TerminalColors as tc

from general.openmm_tools import create_dummy_context


def get_platform_device_from_user(device_property_list, selected_device):
    num_of_devices = len(device_property_list)
    if num_of_devices > 1:
        if selected_device is None:
            selected_device = input(f"Please choose a device (index): \n{tc.TFILE}")
            print(tc.TRESET, end="")
            try:
                selected_device = int(selected_device)
            except:
                print(
                    f'"{tc.TFILE}{selected_device}{tc.TRESET}" not an integer index. Your choices: {list(range(num_of_devices))}'
                )
                raise SystemExit
        if selected_device in range(num_of_devices):
            return selected_device
        else:
            print(
                f'{tc.TRESET}"{tc.TFILE}{selected_device}{tc.TRESET}" not in list of available device indices: {list(range(num_of_devices))}'
            )
            raise SystemExit
    else:
        selected_device = 0
    return selected_device


def print_device_property_list(platform_name, device_property_list):
    print(f"{tc.TRESET}-----------------------------------------")
    tc_color = {}
    tc_color["OpenCL"] = tc.TOCL
    tc_color["CUDA"] = tc.TCUD
    tc_color["CPU"] = tc.TCPU
    color = tc_color[platform_name]

    if platform_name in ["OpenCL", "CUDA"]:
        print(f"Available devices on the {color}{platform_name}{tc.TRESET} platform:")
    elif platform_name == "CPU":
        print(f"{color}{platform_name}{tc.TRESET} platform properties:")
    for i, dict in enumerate(device_property_list):
        print(f"{tc.TBOLD}{color}{i}{tc.TRESET} : ", end="")
        for key, val in dict.items():
            if key not in [
                "DeviceIndex",
                "OpenCLPlatformIndex",
                "TempDirectory",
                "CudaHostCompiler",
            ]:
                print(f"\t{key}\t{color}{val}{tc.TRESET}")
        print(f"{tc.TRESET}-----------------------------------------")


def get_CPU_property_list():
    platform, context = create_dummy_context(platform)
    property_names = platform.getPropertyNames()
    props = {}
    for prop_name in property_names:
        props[prop_name] = platform.getPropertyValue(context, prop_name)
    return [props]


def get_CUDA_property_list():
    device_property_list = []
    list_depth = 20
    for plat_ind in range(list_depth):
        for precision in ["single", "double", "mixed"]:
            properties = {
                "DeviceIndex": str(plat_ind),
                "Precision": precision,
            }
            try:
                platform, context = create_dummy_context(platform, properties)
                property_names = platform.getPropertyNames()
                props = {}
                for prop_name in property_names:
                    props[prop_name] = platform.getPropertyValue(context, prop_name)
                device_property_list.append(props)
            except:
                pass
    return device_property_list


def get_OpenCL_property_list():
    device_property_list = []
    list_depth = 20
    for plat_ind in range(list_depth):
        for dev_ind in range(list_depth):
            for precision in ["single", "double", "mixed"]:
                properties = {
                    "OpenCLPlatformIndex": str(plat_ind),
                    "OpenCLDeviceIndex": str(dev_ind),
                    "Precision": precision,
                }
                try:
                    platform, context = create_dummy_context(platform, properties)
                    property_names = platform.getPropertyNames()
                    props = {}
                    for prop_name in property_names:
                        props[prop_name] = platform.getPropertyValue(context, prop_name)
                    device_property_list.append(props)
                except:
                    pass
    return device_property_list


def get_platform_device_properties(platform_name):
    platform = Platform.getPlatformByName(platform_name)
    platform_name = platform.getName()
    if platform_name == "OpenCL":
        device_property_list = get_OpenCL_property_list()
    elif platform_name == "CUDA":
        device_property_list = get_CUDA_property_list()
    elif platform_name == "CPU":
        device_property_list = get_CPU_property_list()
    if len(device_property_list) == 0 and platform_name != "Reference":
        print(
            f'{tc.TRESET}Could not find available devices on the "{tc.TBOLD}{platform_name}{tc.TRESET}" platform. Check list of available platforms with VCM or check OpenMM\'s installation and plugins.'
        )
        raise
    return device_property_list


def get_platform_device(platform_name, selected_device):
    device_property_list = get_platform_device_properties(platform_name)
    print_device_property_list(platform_name, device_property_list)
    platform_device = get_platform_device_from_user(
        device_property_list, selected_device
    )
    return platform_device


def parse_platform_index(selected_platform):
    platform_names = get_list_of_platforms()
    num_platforms = len(platform_names)
    platform_index, platform_name = None, None
    try:
        user_platform_id = int(selected_platform)
    except:
        return platform_index, platform_name
    if user_platform_id in range(num_platforms):
        platform = Platform.getPlatform(user_platform_id)
        platform_name = platform.getName()
        platform_index = user_platform_id
    else:
        print(
            f'{tc.TRESET}"{tc.TFILE}{user_platform_id}{tc.TRESET}" not in list of available platforms indices. {list(range(num_platforms))}'
        )
        raise SystemExit
    return platform_index, platform_name


def parse_platform_name(selected_platform):
    platform_names = get_list_of_platforms()
    if selected_platform in platform_names:
        platform_name = selected_platform
        platform_index = platform_names.index(selected_platform)
    else:
        print(
            f'{tc.TRESET}"{tc.TFILE}{selected_platform}{tc.TRESET}" not in list of available platforms: {platform_names}'
        )
        raise SystemExit
    return platform_index, platform_name


def parse_selected_platform(selected_platform):
    platform_index, platform_name = parse_platform_index(selected_platform)
    if platform_name is None:
        platform_index, platform_name = parse_platform_name(selected_platform)
    return platform_index, platform_name


def get_list_of_platforms():
    num_platforms = Platform.getNumPlatforms()
    platform_names = []
    for i in range(num_platforms):
        platform = Platform.getPlatform(i)
        platform_names.append(platform.getName())
    return platform_names


def print_available_platforms():
    platform_names = get_list_of_platforms()
    print(
        f"{tc.TOMM}\nOpenMM available platforms:\n{tc.TGRAY}Index Name \t  Speed (Estimated){tc.TRESET}"
    )
    for index, name in enumerate(platform_names):
        platform = Platform.getPlatform(index)
        print(f"({tc.TBOLD}{index}{tc.TRESET})  {name}\t   {platform.getSpeed():0.0f}")


def get_platform_from_user():
    user_platform_input = input(f"Please choose a platform: \n{tc.TFILE}")
    print(tc.TRESET, end="")
    platform_index, platform_name = parse_selected_platform(user_platform_input)
    return platform_index, platform_name


def get_platform_index_and_name(selected_platform):
    if selected_platform is None:
        print_available_platforms()
        platform_index, platform_name = get_platform_from_user()
    else:
        platform_index, platform_name = parse_selected_platform(selected_platform)
    return platform_index, platform_name


def print_available_platforms(user_selected_platform, user_selected_device):
    platform_index, platform_name = get_platform_index_and_name(user_selected_platform)
    user_flags = f"\nSelected device flags:\n--platform {platform_name}"
    if platform_name != "Reference":
        platform_device_ID = get_platform_device(platform_name, user_selected_device)
        user_flags += f" --platform-device-ID {platform_device_ID}"
    print(user_flags)

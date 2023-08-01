from openmm.openmm import Platform
from general.classes.print_colors import TerminalColors as tc

from general.openmm_tools import create_dummy_context


def check_index_type(index):
    """Check if input is integer.

    Args:
        index (str, int): String or integer representing an index.

    Raises:
        ValueError: When input cannot be converted into an int.

    Returns:
        int: index
    """
    try:
        index = int(index)
    except:
        raise ValueError(f'"{tc.TFILE}{index}{tc.TRESET}" not an integer index.')
    return index


def get_user_input(message):
    """Get platform device index from the user.

    Raises:
        ValueError: When user input is not an int.

    Returns:
        str: user selected platform device index.
    """
    user_input = input(f"{message}\n{tc.TFILE}")
    print(tc.TRESET, end="")
    return user_input


def print_device_property_list(platform_name, device_property_list):
    """Print the properties of devices belonging to a platform.

    Go through property names and values of devices on the platform and print them on the screen together with an index for user selection/identification.

    Args:
        platform_name (str): platform name (Reference, CPU, CUDA, or OpenCL).
        device_property_list (list): list of dictionaries containing property name and values of devices that belong to the platform.
    """
    print(f"{tc.TRESET}-----------------------------------------")
    color = tc.tc_color[platform_name]
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


def get_CPU_property_list(platform):
    """Get the properties of the CPU platform.

    Args:
        platform (openmm.Platform): OpenMM platform referencing the CPU platform.

    Returns:
        list: list of dictionaries containing property name and values of devices that belong to the platform.
    """
    platform, context = create_dummy_context(platform)
    property_names = platform.getPropertyNames()
    props = {}
    for prop_name in property_names:
        props[prop_name] = platform.getPropertyValue(context, prop_name)
    return [props]


def get_CUDA_property_list(platform):
    """Get the properties of the CUDA platform.

    Args:
        platform (openmm.Platform): OpenMM platform referencing the CUDA platform.

    Returns:
        list: list of dictionaries containing property name and values of devices that belong to the platform.
    """
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


def get_OpenCL_property_list(platform):
    """Get the properties of the OpenCL platform.

    Args:
        platform (openmm.Platform): OpenMM platform referencing the OpenCL platform.

    Returns:
        list: list of dictionaries containing property name and values of devices that belong to the platform.
    """
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
    """Get the list of properties for each available device on the platform.

    Args:
        platform_name (str): name of the platform.

    Returns:
        list: list of dictionaries containing property name and values of devices that belong to the platform.
    """
    platform = Platform.getPlatformByName(platform_name)
    platform_name = platform.getName()
    device_property_list = []
    if platform_name == "OpenCL":
        device_property_list = get_OpenCL_property_list(platform)
    elif platform_name == "CUDA":
        device_property_list = get_CUDA_property_list(platform)
    elif platform_name == "CPU":
        device_property_list = get_CPU_property_list(platform)
    if len(device_property_list) == 0 and platform_name != "Reference":
        print(
            f'{tc.TRESET}Could not find available devices on the "{tc.TBOLD}{platform_name}{tc.TRESET}" platform. Check list of available platforms with VCM or check OpenMM\'s installation and plugins.'
        )
        raise
    return device_property_list


def parse_platform_device_selection(device_index_choices, selected_device):
    """Parse user selected device.

    Args:
        device_index_choices (list): List of integers containing indices of platform's available devices.
        selected_device (str, int): Input by user referring to a platform device index.

    Raises:
        ValueError: If input is not an integer or not in the list of device choices.

    Returns:
        int: Index of the platform's device.
    """
    if selected_device is None:
        selected_device = get_user_input("Please choose a device (index): ")
    selected_device = check_index_type(selected_device)
    if selected_device in device_index_choices:
        platform_device = selected_device
    else:
        raise ValueError(
            f'{tc.TRESET}"{tc.TFILE}{selected_device}{tc.TRESET}" not in list of available device indices: {device_index_choices}'
        )
    return platform_device


def get_platform_device(platform_name, selected_device):
    """Get platform's device index from user selection.

    Args:
        platform_name (str): Platform name (CPU, OpenCL, etc)
        selected_device (str, None): Platform's device index.

    Returns:
        int: Index of the device on the platform.
    """
    device_property_list = get_platform_device_properties(platform_name)
    print_device_property_list(platform_name, device_property_list)
    num_of_devices = len(device_property_list)
    if num_of_devices > 1:
        platform_device = parse_platform_device_selection(
            list(range(num_of_devices)), selected_device
        )
    else:
        platform_device = 0
    return platform_device


def parse_platform_index(selected_platform):
    """Get platform name and index from input index.

    Args:
        selected_platform (str): selected platform index.

    Raises:
        ValueError: When index is out or range of platform indices.

    Returns:
        platform_index: int
            Index of the platform. 'None' if selected_platform is not an index.
        platform_name: str
            Name of the platform. 'None' if selected_platform is not an index.
    """
    platform_names = get_list_of_platform_names()
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
        raise ValueError(
            f'{tc.TRESET}"{tc.TFILE}{user_platform_id}{tc.TRESET}" not in list of available platforms indices. {list(range(num_platforms))}'
        )
    return platform_index, platform_name


def parse_platform_name(selected_platform):
    """Parse the input platform name.

    Args:
        selected_platform (str): Name of the selected platform.

    Raises:
        ValueError: If name is not in the list of available platforms.

    Returns:
        platform_index: int
            Index of the selected platform.
        platform_name: str
            Name of the selected platform.
    """
    platform_names = get_list_of_platform_names()
    if selected_platform in platform_names:
        platform_name = selected_platform
        platform_index = platform_names.index(selected_platform)
    else:
        raise ValueError(
            f'{tc.TRESET}"{tc.TFILE}{selected_platform}{tc.TRESET}" not in list of available platforms: {platform_names}'
        )
    return platform_index, platform_name


def parse_selected_platform(selected_platform):
    """Parse the selected platform input.

    Args:
        selected_platform (str): The index or name of the selected platform.

    Returns:
        platform_index: int
            Index of the selected platform.
        platform_name: str
            Name of the selected platform.
    """
    platform_index, platform_name = parse_platform_index(selected_platform)
    if platform_name is None:
        platform_index, platform_name = parse_platform_name(selected_platform)
    return platform_index, platform_name


def get_list_of_platform_names():
    """Get a list of available platform names.

    Returns:
        list: List of strings of platform names.
    """
    num_platforms = Platform.getNumPlatforms()
    platform_names = []
    for i in range(num_platforms):
        platform = Platform.getPlatform(i)
        platform_names.append(platform.getName())
    return platform_names


def print_available_platforms():
    """Print available platforms on the machine."""
    platform_names = get_list_of_platform_names()
    print(
        f"{tc.TOMM}\nOpenMM available platforms:\n{tc.TGRAY}Index Name \t  Speed (Estimated){tc.TRESET}"
    )
    for index, name in enumerate(platform_names):
        platform = Platform.getPlatform(index)
        print(f"({tc.TBOLD}{index}{tc.TRESET})  {name}\t   {platform.getSpeed():0.0f}")


def get_platform_index_and_name(selected_platform):
    """Get the platform name and index from user inputs.

    Args:
        selected_platform (str, None): User argparse input for the platform selection.

    Returns:
        platform_index: int
                Index of the platform in the list of available platforms.
        platform_name: str
                Name of the platform in the list of available platforms.

    """
    if selected_platform is None:
        print_available_platforms()
        selected_platform = get_user_input("Please choose a platform: ")
    platform_index, platform_name = parse_selected_platform(selected_platform)
    return platform_index, platform_name


def print_available_platforms_and_devices(user_selected_platform, user_selected_device):
    """Print the VCM input flags corresponding to user choice of platform and device.

    Args:
        user_selected_platform (str, None): User argparse input for the platform selection.
        user_selected_device (_type_, none): User argparse input for the platform device index selection.
    """
    platform_index, platform_name = get_platform_index_and_name(user_selected_platform)
    user_flags = f"\nSelected device flags:\n--platform {platform_name}"
    if platform_name != "Reference":
        platform_device_ID = get_platform_device(platform_name, user_selected_device)
        user_flags += f" --platform-device-ID {platform_device_ID}"
    print(user_flags)

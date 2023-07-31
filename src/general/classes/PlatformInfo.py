import sys
import openmm
from openmm.openmm import Platform
from openmm.openmm import System
from openmm.openmm import VerletIntegrator
from openmm.openmm import Context
import openmm.app.element as elem
import openmm.unit as units
from openmm.app.topology import Topology
from openmm.app.topology import Residue
from openmm.app.simulation import Simulation


from general.classes.print_colors import TerminalColors as tc


class PlatformInfo:
    def __init__(self):
        self.name = None
        self.index = None
        self.device_ID = 0
        self.openmm_plugin_path = None
        self.device_properties = []
        self.device_properties_report = []


def get_list_of_platforms():
    num_platforms = Platform.getNumPlatforms()
    platform_names = []
    for i in range(num_platforms):
        platform = Platform.getPlatform(i)
        platform_names.append(platform.getName())
    return platform_names


def parse_selected_platform(selected_platform):
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


def get_platform_device_properties(platform_name):
    platform = Platform.getPlatformByName(platform_name)

    device_property_list = []
    list_depth = 20
    system = System()
    system.addParticle(1.0)
    step_size_Ps = 0.001 * units.picosecond
    integrator = VerletIntegrator(step_size_Ps)
    if platform.getName() == "OpenCL":
        for plat_ind in range(list_depth):
            for dev_ind in range(list_depth):
                for precision in ["single", "double", "mixed"]:
                    properties = {
                        "OpenCLPlatformIndex": str(plat_ind),
                        "OpenCLDeviceIndex": str(dev_ind),
                        "Precision": precision,
                    }
                    try:
                        context = Context(system, integrator, platform, properties)
                        property_names = platform.getPropertyNames()
                        props = {}
                        for prop_name in property_names:
                            props[prop_name] = platform.getPropertyValue(
                                context, prop_name
                            )
                        device_property_list.append(props)
                    except:
                        pass
    elif platform.getName() == "CUDA":
        for plat_ind in range(list_depth):
            for precision in ["single", "double", "mixed"]:
                properties = {
                    "DeviceIndex": str(plat_ind),
                    "Precision": precision,
                }
                try:
                    context = Context(system, integrator, platform, properties)
                    property_names = platform.getPropertyNames()
                    props = {}
                    for prop_name in property_names:
                        props[prop_name] = platform.getPropertyValue(
                            context, prop_name
                        )
                    device_property_list.append(props)
                except:
                    pass
    elif platform.getName() == "CPU":
        context = Context(system, integrator, platform)
        property_names = platform.getPropertyNames()
        props = {}
        for prop_name in property_names:
            props[prop_name] = platform.getPropertyValue(context, prop_name)
        device_property_list.append(props)
    return device_property_list


def get_platform_device_from_user(device_property_list, selected_device):
    num_of_devices = len(device_property_list)
    if num_of_devices > 1:
        if selected_device is None:
            selected_device = input("Please choose a device (index): \n{tc.TFILE}")
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
    print(f'{tc.TRESET}-----------------------------------------')
    if platform_name == "OpenCL":
        tc_color = tc.TOCL
        print(f'Available devices on the {tc_color}{platform_name}{tc.TRESET} platform:')
    elif platform_name == "CUDA":
        tc_color = tc.TCUD
        print(f'Available devices on the {tc_color}{platform_name}{tc.TRESET} platform:')
    elif platform_name == "CPU":
        tc_color = tc.TCPU
        print(f'{tc_color}{platform_name}{tc.TRESET} platform properties:')
    for i, dict in enumerate(device_property_list):
        print(f'{tc.TBOLD}{tc_color}{i}{tc.TRESET} : ',end="")
        for key, val in dict.items():
            print(f'\t{key}\t{tc_color}{val}{tc.TRESET}')
        print(f'{tc.TRESET}-----------------------------------------')

    

def get_platform_device(platform_name, selected_device):
    device_property_list = get_platform_device_properties(platform_name)
    print_device_property_list(platform_name, device_property_list)
    platform_device = get_platform_device_from_user(
        device_property_list, selected_device
    )
    return platform_device


def print_available_platforms(user_inputs):
    platform_index, platform_name = get_platform_index_and_name(
        user_inputs.user_selected_platform
    )
    user_flags = f'Selected device flags:\n--platform {platform_name}'
    if platform_name != "Reference":
        platform_device_ID = get_platform_device(
            platform_name, user_inputs.user_selected_device
        )
        user_flags += f' --platform-device-ID {platform_device_ID}'


from openmm.openmm import Platform
from general.classes.print_colors import TerminalColors as tc


class PlatformInfo:
    def __init__(self):
        self.platform_name = "Reference"
        self.platform_device_ID = 0
        self.openmm_plugin_path = None


def print_platform_info(user_inputs):
    numPlatforms = Platform.getNumPlatforms()
    print(
        f"{tc.TOMM}\nOpenMM available platforms:\n{tc.TGRAY}Index Name \t  Speed (Estimated){tc.TRESET}"
    )
    for i in range(numPlatforms):
        platform = Platform.getPlatform(i)
        print(f"({tc.TBOLD}{i}{tc.TRESET})  {platform.getName()}\t   {platform.getSpeed():0.0f}")


from openmm.openmm import Platform

class PlatformInfo:
    def __init__(self):
        self.platform_name = "Reference"
        self.platform_device_ID = 0
        self.openmm_plugin_path = None


def print_platform_info(user_inputs):
    numPlatforms = Platform.getNumPlatforms()
    print("There are", numPlatforms, "Platforms available:")
    print()
    for i in range(numPlatforms):
        platform = Platform.getPlatform(i)
        print(i+1, platform.getName())
    
    
    import sys
    sys.exit()


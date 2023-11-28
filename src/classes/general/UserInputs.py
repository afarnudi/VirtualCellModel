from classes.general.PlatformInfo import PlatformInfo


class UserInputs:
    def __init__(self):
        self.user_selected_platform = None
        self.user_selected_device = None
        self.platform_info = PlatformInfo()
        self.resume_path = None
        self.config_file_path = None
        self.write_at_end = False

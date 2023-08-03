class Settings:
    def __init__(self, default_value=None, help=None):
        self._default_value = default_value
        self._help = Settings.check_description(help)

    @property
    def help(self):
        return self._help

    @property
    def default_value(self):
        return self._default_value

    @help.setter
    def help(self, string):
        self._help = Settings.check_description(string)

    @default_value.setter
    def default_value(self, value):
        self._default_value = value

    def check_description(help):
        if help.startswith("#") is False:
            help = "#"+help
        return help
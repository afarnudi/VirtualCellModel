class Settings:
    def __init__(self, default_value=None, help=None):
        self._value = default_value
        self._help = Settings.check_description(help)

    @property
    def help(self):
        return self._help

    @property
    def value(self):
        return self._value

    @help.setter
    def help(self, string):
        self._help = Settings.check_description(string)

    @value.setter
    def value(self, val):
        self._value = val

    def check_description(help):
        if help.startswith("#") is False:
            help = "#"+help
        return help
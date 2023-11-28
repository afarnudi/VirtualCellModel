"""VCM configuration and settings containers.

Provide storage units for VCM configurations and Section settings.
"""
from classes.general.print_colors import TerminalColors as tc


class Settings:
    """Class storing values and help descriptions."""

    def __init__(self, default_value=None, help=None):
        """Initilise and instance.

        Args:
            default_value (str, optional): Setting value. Defaults to None.
            help (str, optional): Description for the Setting value, options, and format. Defaults to None.
        """
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
        """Check description format.

        Add a '#' character to the beginning of the description if it does not exist.

        Args:
            help (str): Description of the format of the value.

        Returns:
            string: Value description.
        """
        if help.startswith("#") is False:
            help = "#" + help
        return help


class Configuration:
    """A dictionary based class to store user selected VCM configuration.

    Create an instance of configuration only if the section title is registered with VCM.
    """

    """Store VCM user selected configurations.

    Store and preprocess user configurations for VCM's different sections.

    Raises:
        TypeError: If section name does not follow the "-SectionName" format and it is not placed in the beginning of the line.
        KeyError: If the configuration is asked to return the value of a configuration it does not hold.

    Returns:
        _type_: _description_
    """

    SECTION_LIST = [
        "GeneralParameters",
        "Membrane",
        "Actin",
        "Chromatin",
        "InteractionTable",
    ]
    section_declaration_message = f'Section declarations must be at the beginning of the line and in the following format "-SectionName". SectionName choices are: {SECTION_LIST}'

    def __init__(self, section_title_string):
        """Initialise an instance.

        Args:
            section_title_string (str): One of VCM's section name.

        Raises:
            TypeError: If section name does not follow the "-SectionName" format and it is not placed in the beginning of the line.
        """
        if self.is_section_declaration(section_title_string) is False:
            raise TypeError(
                f'{tc.TRESET}{tc.TFAILED}Configuration file parsing: "{tc.TFILE}{section_title_string}{tc.TFAILED}" is not a section name.\n{Configuration.section_declaration_message}{tc.TRESET}'
            )
        self.title = self._get_title(section_title_string)
        self.configurations = {}

    def _get_title(self, title):
        """Extract the VCM section name.

        Args:
            title (str): VCM section name together with the character decorators.

        Returns:
            str: VCM section name
        """
        if self.is_section_declaration(title):
            return title.split()[0][1:]

    def is_section_declaration(self, line):
        """Find section deceleration pattern.

        Args:
            line (str): A line containing VCM section deceleration.

        Raises:
            TypeError: If section name is not in VCM's list of valid names.

        Returns:
            bool: True if the section deceleration is valid. False if the section name pattern is not recognised.
        """
        line = self.clean_line(line)
        if line.startswith("-"):
            start = line.split()[0]
            if start[1:] in Configuration.SECTION_LIST:
                return True
            else:
                raise TypeError(
                    f'{tc.TRESET}{tc.TFAILED}Configuration file parsing: "{tc.TFILE}{start}{tc.TFAILED}" is not a section name.\n{Configuration.section_declaration_message}{tc.TRESET}'
                )
        return False

    def add(self, line):
        """Add a section configuration.

        Extract configuration name and value

        Args:
            line (str): _description_

        Raises:
            ValueError: If the configuration does not have a value in front of it.
        """
        line = self.clean_line(line)
        if len(line) != 0:
            key = line.split()[0]
            if len(line.split()) == 1:
                raise ValueError(
                    f'{tc.TRESET}{tc.TFAILED}Configuration file parsing: "{tc.TFILE}{start}{tc.TFAILED}" does not state a value for.\n{key}. Consult the template and try agin.{tc.TRESET}'
                )
            setting = Settings(
                default_value=line.split(key)[-1].rstrip().lstrip(), help="#"
            )
            self.configurations[key] = setting

    def __str__(self):
        out = f"{self.title}\n"
        for k, v in self.configurations.items():
            out += f"{k}: {v}\n"
        return out

    def get_name(self):
        return self.title

    def get_value(self, key):
        """Get configuration value.

        Args:
            key (str): configuration option name.

        Raises:
            KeyError: If the instance does not hold this configuration.

        Returns:
            str: Value of the configuration.
        """
        keys = self.get_key_names()
        if key in keys:
            return self.configurations[key].value
        raise KeyError(
            f'The configuration does not contain a "{tc.TFILE}{key}{tc.TRESET}" key.'
        )

    def get_key_names(self):
        return self.configurations.keys()

    def clean_line(self, line):
        """Clean up the configuration line.

        Remove comments and trailing white space on the right and left of the line.

        Args:
            conf_lines (str):Configuration line.

        Returns:
            str: Cleaned up string.
        """
        line = self.strip_of_comments(line)
        line = line.rstrip()
        line = line.lstrip()

        return line

    def strip_of_comments(self, line, comment="#"):
        """Remove characters to the right of the comment symbol (including the comment symbol).

        Args:
            line (str): String containing statements and/or comments.
            comment (str, optional): String identified as the comment symbol. Defaults to "#".

        Returns:
            str: String content that is on the left of the comment symbol.
        """
        return line.split(comment)[0]


SECTION_LIST = Configuration.SECTION_LIST.copy()
# ----------------------------------------------------------------------
# Create one instance of the Configuration file, and export its methods
# as module-level functions.

_inst = Configuration("-GeneralParameters")
strip_of_comments = _inst.strip_of_comments
clean_line = _inst.clean_line
is_section_declaration = _inst.is_section_declaration

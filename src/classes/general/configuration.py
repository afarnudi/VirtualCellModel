from classes.general.print_colors import TerminalColors as tc


class Configuration:
    SECTION_LIST = [
        "GeneralParameters",
        "Membrane",
        "Actin",
        "Chromatin",
        "InteractionTable",
    ]
    section_declaration_message = f'Section declarations must be at the beginning of the line and in the following format "-SectionName". SectionName choices are: {SECTION_LIST}'

    def __init__(self, section_title_string):
        self.title = self._get_title(section_title_string)
        self.configurations = {}

    def _get_title(self, title):
        if Configuration.check_string_for_section_declaration(title):
            return title.split()[0][1:]
        else:
            raise TypeError(
                f"{tc.TFAILED}Configuration file parsing: {Configuration.section_declaration_message}{tc.TRESET}"
            )

    def check_string_for_section_declaration(line):
        start = line.split()[0]

        if start.startswith("-"):
            if start[1:] in Configuration.SECTION_LIST:
                return True
            else:
                raise TypeError(
                    f'{tc.TRESET}{tc.TFAILED}Configuration file parsing: "{tc.TFILE}{start}{tc.TFAILED}" is not a section name.\n{Configuration.section_declaration_message}{tc.TRESET}'
                )
        return False

    def add(self, line):
        key = line.split()[0]
        self.configurations[key] = line.split(key)[-1]

    def __str__(self):
        out = f"{self.title}\n"
        for k, v in self.configurations.items():
            out += f"{k}: {v}\n"
        return out

    def get_name(self):
        return self.title
    
    
    def get_value(self, key):
        keys = self.get_key_names()
        if key in keys:
            return self.configurations[key]
        raise KeyError(f'The configuration does not contain a "{tc.TFILE}{key}{tc.TRESET}" key.')
    
    def get_key_names(self):
        return self.configurations.keys()


from general.classes.print_colors import TerminalColors as tc


class Configuration:
    CLASS_LIST = [
        "GeneralParameters",
        "Membrane",
        "Actin",
        "Chromatin",
        "InteractionTable",
    ]

    def __init__(self, class_title_string):
        self.title = self._get_title(class_title_string)
        self.configurations = {}

    def _get_title(self, title):
        if Configuration.check_string_for_class_declaration(title):
            print(title.split()[0][1:])
            return title.split()[0][1:]
        else:
            raise TypeError(
                f'{tc.TFAILED}Configuration file parsing: Class declarations must be at the beginning of the line and in the following format "-ClassName". ClassName choices are: {self.CLASS_LIST}{tc.TRESET}'
            )

    def check_string_for_class_declaration(line):
        start = line.split()[0]
        if start.startswith("-"):
            if start[1:] in Configuration.CLASS_LIST:
                return True
            else:
                raise TypeError(
                    f'{tc.TRESET}{tc.TFAILED}Configuration file parsing: "{tc.TFILE}{start}{tc.TFAILED}" is not a class name.\nClass declarations must be at the beginning of the line and in the following format "-ClassName". ClassName choices are: {Configuration.CLASS_LIST}{tc.TRESET}'
                )
        return False

    def add(self, line):
        key = line.split()[0]
        value = line.split(key)[-1]
        self.configurations[key] = value.split()

    def __str__(self):
        out = f"{self.title}\n"
        for k, v in self.configurations.items():
            out += f"{k}: {v}\n"
        return out

    def get_class_name(self):
        return self.title
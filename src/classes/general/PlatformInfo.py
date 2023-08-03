from openmm.openmm import Context
import openmm.app.element as elem

from openmm.app.topology import Topology
from openmm.app.topology import Residue
from openmm.app.simulation import Simulation


from classes.general.print_colors import TerminalColors as tc


class PlatformInfo:
    def __init__(self):
        self.name = None
        self.index = None
        self.device_ID = 0
        self.openmm_plugin_path = None
        self.device_properties = []
        self.device_properties_report = []



























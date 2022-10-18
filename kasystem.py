# Walter Fontana, 2022
"""
This module defines the 'System' object.
'System' holds overall state variables and pointers to data structures of the system.
"""
import datetime
import uuid
import sys

system = None


class System:
    def __init__(self):
        """
        Initializes the 'System' object, which holds overall state variables and
        pointers to data structures.
        """
        self.date = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
        self.uuid = uuid.uuid1()

        self.cmdline = None

        # program behavior directives
        self.db_level = 0
        self.consolidate = True
        self.canonicalize = True

        self.kappa = None           # Kappa parser
        self.sgm = None             # Site Graph Matcher
        self.signature = None       # system signature
        self.parameters = None      # system parameters
        self.mixture = None         # mixture of kappa expressions
        self.sim = None             # ctmc simulator
        self.monitor = None         # observable monitor

        self.sim_limit_type = 'time'  # {time, event}
        self.sim_limit = 0.

        self.report_file = None
        self.parameter_file = None
        self.mixture_file = None
        self.signature_string = ''  # signature string

        # reaction rate constants
        self.rc_bond_formation_intra = 0.
        self.rc_bond_formation_inter = 0.
        self.rc_bond_dissociation = {}  # indexed by bond type

        self.inflow_rate = {}   # indexed by atom type
        self.outflow_rate = {}  # indexed by atom type

    def report(self):
        with open(self.report_file, "w") as report:
            if not self.parameters:
                print(f'No parameters!')
                sys.exit()
            else:
                param_info = str(self.parameters)

            if not self.signature:
                print(f'No signature!')
                sys.exit()
            else:
                sig_info = str(self.signature)

            if not self.mixture:
                print(f'No initial mixture!')
                sys.exit()
            else:
                mix_info = str(self.mixture)

            sim_info = str(self.sim)

            report.write(f'\n\n{"date and time (UTC)":>30}: {self.date}\n')
            report.write(f'{"uuid":>30}: {self.uuid}\n\n')
            info = f"\n{'COMMAND LINE '.ljust(70, '-')}\n\n"
            report.write(f'{info}{self.cmdline}\n')

            report.write(sig_info)
            report.write(param_info)
            report.write(sim_info)
            report.write(mix_info)


def init_system():
    global system
    system = System()  # global scope

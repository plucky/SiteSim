# Walter Fontana, 2022
"""
This module defines the 'System' object.
'System' holds overall state variables and pointers to data structures of the system.
"""
from datetime import datetime, timezone
from zoneinfo import ZoneInfo
import time
import psutil
import uuid
import sys
import kainit

system = None


class System:
    def __init__(self):
        """
        Initializes the 'System' object, which holds overall state variables and
        pointers to data structures.
        """
        self.date_utc = datetime.now(timezone.utc)
        self.date_eastern = self.date_utc.astimezone(tz=ZoneInfo("America/New_York"))
        self.uuid = uuid.uuid1()

        self.cmdline = None

        # program behavior directives
        self.db_level = 0
        self.consolidate = True
        self.canonicalize = True
        self.barcode = False

        self.kappa = None           # Kappa parser
        self.sgm = None             # Site Graph Matcher
        self.signature = None       # system signature
        self.parameters = None      # system parameters
        self.mixture = None         # mixture of kappa expressions
        self.sim = None             # ctmc simulator
        self.monitor = None         # monitor of observables
        self.alarm = None           # stopping conditions

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

            sys_info = f'\n\n{"initialized (UTC)":>30}: {self.date_utc.strftime("%Y-%m-%d %H:%M:%S")}\n'
            sys_info += f'{"initialized (Boston)":>30}: {self.date_eastern.strftime("%Y-%m-%d %H:%M:%S")}\n'
            sys_info += f'{"uuid":>30}: {self.uuid}\n\n'
            sys_info += f"\n{'COMMAND LINE '.ljust(70, '-')}\n\n"
            sys_info += f'{self.cmdline}\n'
            sys_info += f"\n{'SYSTEM SETTINGS '.ljust(70, '-')}\n\n"
            sys_info += f'{"consolidate":>20}: {self.consolidate}\n'
            sys_info += f'{"canonicalize":>20}: {self.canonicalize}\n'
            sys_info += f'{"barcode":>20}: {self.barcode}\n\n'

            report.write(sys_info)
            report.write(sig_info)
            report.write(param_info)
            report.write(sim_info)
            report.write(mix_info)

    def resources_report(self):
        end_utc = datetime.now(timezone.utc)
        end_eastern = end_utc.astimezone(tz=ZoneInfo("America/New_York"))
        with open(self.report_file, "a") as report:
            sys_info = f"\n{'RESOURCES '.ljust(70, '-')}\n\n"
            sys_info += f'{"uuid":>30}: {self.uuid}\n'
            sys_info += f'{"terminated (UTC)":>30}: {end_utc.strftime("%Y-%m-%d %H:%M:%S")}\n'
            sys_info += f'{"terminated (Boston)":>30}: {end_eastern.strftime("%Y-%m-%d %H:%M:%S")}\n'
            sys_info += f'{"real time since initializing":>30}: {end_eastern - self.date_eastern}\n'
            sys_info += f'{"events":>30}: {self.sim.event}\n'
            sys_info += f'{"cpu":>30}: {time.process_time():.3f} seconds\n'
            process = psutil.Process()
            sys_info += f'{"memory rss":>30}: {process.memory_info().rss / (1024 * 1024):.3f} Mb\n'
            sys_info += f'{"memory vms":>30}: {process.memory_info().vms / (1024 * 1024):.3f} Mb\n'
            report.write(sys_info)

def init_system():
    global system
    system = System()  # global scope

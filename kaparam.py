# Walter Fontana, 2022
"""
This module defines the 'Parameters' object.
'Parameters' holds the system parameters.
"""
# from scipy import constants
import os
import re
import sys

import kasystem as ka
import kasig
import math


def is_number(s):
    """
    Checks if s is a number. Returns True or False.
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


class Parameters:
    """
    Stores the system parameters.
    """
    def __init__(self, file=None):

        # constants

        self.Avogadro = 6.02214e+23
        self.GasConstant_SI = 8.31446261815324  #  J K^{-1} mol^{-1}
        self.GasConstant_cal = 1.98720425864083  #  cal K^{-1} mol^{-1} (1 J = 0.2390057361 cal, 1 cal = 4.184 J)
        self.Volume_choices = {'fibro': 2.25e-12, 'yeast': 4.2e-14}  # L

        self.referenceVol = self.Volume_choices['fibro']
        self.referenceTemp = 273.15 + 25.0  # K
        self.referenceTemp_C = 25.0  # C  (IUPAC "room temperature")

        self.Volume = self.referenceVol
        self.Temperature = self.referenceTemp  # K
        self.default_concentration = 100   # default initial agent concentration in nM

        # default parameters------------------------------------------------------

        self.Kd_weak   = 1000.e-9  #   1 mumol
        self.Kd_medium =  100.e-9  # 100 nmol
        self.Kd_strong =   10.e-9  #  10 nmol

        self.k_on = 1.e+9  # (M s)^-1; diffusion-controlled limit

        self.RingClosureFactor = 1.e+5  # ratio of binary Kd to unary Kd

        self.ResizeVolume = 1.
        self.RescaleTemperature = 1.

        # in/out rates -----------------------------------------------------------

        self.inflow = {}     # for each atom type
        self.outflow = {}    # for each atom type

        # stochastic rate constants ----------------------------------------------

        self.s_on = 0.
        self.s_ring_on = 0.
        self.s_off_weak = 0.
        self.s_off_medium = 0.
        self.s_off_strong = 0.

        # initial abundances (in molecules) --------------------------------------

        self.init_agents = {}

        # random number generator seed -------------------------------------------

        self.rng_seed = None

        # signature string -------------------------------------------------------

        self.signature_string = None

        # update and finalize the parameters
        if file:
            self.read_parameters(file)
        else:
            sys.exit("no parameter file.")

        # process SIGNATURE
        if not ka.system.signature_string:  # this would have come from the commandline
            if self.signature_string:
                ka.system.signature_string = self.signature_string
            else:
                sys.exit('signature information is missing')
        ka.system.signature = kasig.KappaSignature(ka.system.signature_string)

        # sanity check
        if self.inflow and not ka.system.canonicalize:
            sys.exit("in/out flow requires canonicalization.")

    def apply_parameters(self):
        """
        Apply the system volume and temperature scale factor to the parameters.
        """
        # sanity
        if ka.system.barcode:  # if we barcode, there's no point in consolidating
            ka.system.consolidate = False
        if not ka.system.consolidate:  # if we don't consolidate, we might as well not canonicalize
            ka.system.canonicalize = False

        if self.Volume != self.referenceVol:  # the inputted volume has precedence
            print(f'Using volume setting and adjusting scale factor relative to default reference')
            # self.Volume = self.Volume
            self.ResizeVolume = self.Volume / self.referenceVol
        else:
            self.Volume = self.referenceVol * self.ResizeVolume

        if self.Temperature != self.referenceTemp:  # the inputted temperature has precedence
            print(f'Using temperature setting and adjusting scale factor relative to default reference')
            # self.Temperature = self.Temperature
            self.RescaleTemperature = self.Temperature / self.referenceTemp
        else:
            self.Temperature = self.referenceTemp * self.RescaleTemperature

        # this assumes an ideal mono-atomic gas...
        self.RingClosureFactor = self.RingClosureFactor * self.ResizeVolume * math.pow(self.RescaleTemperature, 3./2.)

        # stochastic rate constants-----------------------------------------------

        # bi-molecular on-constant  (inter-binding)
        self.s_on = self.k_on / (self.Avogadro * self.Volume)
        # uni-molecular on-constant (intra-binding)
        self.s_ring_on = self.RingClosureFactor * self.s_on

        # off-constants are the same regardless of whether they lead to a fission
        self.s_off_weak = self.k_on * math.pow(self.Kd_weak, 1. / self.RescaleTemperature)
        self.s_off_medium = self.k_on * math.pow(self.Kd_medium, 1. / self.RescaleTemperature)
        self.s_off_strong = self.k_on * math.pow(self.Kd_strong, 1. / self.RescaleTemperature)

        # better mnemonic, stored at 'system'
        ka.system.rc_bond_formation_inter = self.s_on
        ka.system.rc_bond_formation_intra = self.s_ring_on

        ka.system.inflow_rate = {}
        ka.system.outflow_rate = {}
        for at in self.inflow:
            ka.system.inflow_rate[at] = self.inflow[at] * self.Avogadro * self.Volume
        for at in self.outflow:
            ka.system.outflow_rate[at] = self.outflow[at]

        # These are dissociation constants (Kd)!
        # However, terms like w(eak), m(edium), s(trong), def(ault) refer to affinities (1/Kd);
        # thus, "weak" means a high Kd
        ka.system.rc_bond_dissociation = {}
        for bt in ka.system.signature.bond_types:
            if is_number(ka.system.signature.bond_types[bt]):
                Kd = float(ka.system.signature.bond_types[bt]) * 1.e-9  # at reference temperature
                Kd = math.pow(Kd, 1. / self.RescaleTemperature)
                ka.system.rc_bond_dissociation[bt] = self.k_on * Kd
            else:
                if ka.system.signature.bond_types[bt] == 'w':  # weak affinity
                    ka.system.rc_bond_dissociation[bt] = self.s_off_weak
                elif ka.system.signature.bond_types[bt] == 'm':  # medium affinity
                    ka.system.rc_bond_dissociation[bt] = self.s_off_medium
                elif ka.system.signature.bond_types[bt] == 's':  # strong affinity
                    ka.system.rc_bond_dissociation[bt] = self.s_off_strong
                elif ka.system.signature.bond_types[bt] == 'def':  # default affinity
                    ka.system.rc_bond_dissociation[bt] = self.s_off_medium
                else:  # default affinity
                    ka.system.rc_bond_dissociation[bt] = self.s_off_medium

        # initial agent counts
        for a in ka.system.signature.init_agents:
            init = ka.system.signature.init_agents[a]
            if init == '*':  # default abundance in nM
                self.init_agents[a] = int(self.default_concentration * 1.e-9 * self.Avogadro * self.Volume)
            else:
                self.init_agents[a] = int(float(init) * 1.e-9 * self.Avogadro * self.Volume)

    def read_parameters(self, par_file):
        """
        Reads parameters from a file.

        par_file: file name containing parameter values
        Returns: nothing

        Current declaration keywords:
            %par, %sig, %rep, (%obs declarations are read in kamon.py)
        Current keywords:
            'Volume', 'Kd_weak', 'Kd_medium', 'Kd_strong', 'k_on', 'Resize', 'RingClosureFactor',
            'seed', 'inflow', 'outflow', 'sim_limit', 'obs_frequency', 'report_fn', 'snap_root',
            'output_fn', 'numbering'
        """
        if not os.path.isfile(par_file):
            sys.exit("Cannot find parameter file %s" % par_file)
        else:
            with open(par_file, "r", encoding='utf-8') as data:
                while True:
                    line = data.readline()
                    if not line:
                        break
                    # parse the line
                    match = re.match(r'^%(par:|sig:|rep:)\s?', line)
                    if match:
                        if match.group(1) == 'par:':
                            match = re.match(r'%par:\s*(.*)\s*=\s*(\S*)\s?', line)
                            if match:
                                name = match.group(1).strip()
                                value = match.group(2).strip()
                                if name == 'Volume':
                                    if value in self.Volume_choices:
                                        self.Volume = self.Volume_choices[value]
                                    else:
                                        if is_number(value):
                                            self.Volume = float(value)
                                        else:
                                            sys.exit(f'No such volume choice: {value}')
                                elif name == 'Temperature':  # in C
                                    self.Temperature = float(value) + 273.15
                                elif name == 'ReferenceVolume':
                                    if value in self.Volume_choices:
                                        self.referenceVol = self.Volume_choices[value]
                                    else:
                                        if is_number(value):
                                            self.referenceVol = float(value)
                                        else:
                                            sys.exit(f'No such volume choice: {value}')
                                elif name == 'ReferenceTemp':  # in C
                                    self.referenceTemp = float(value) + 273.15
                                elif name == 'Kd_weak':
                                    self.Kd_weak = float(value)
                                elif name == 'Kd_medium':
                                    self.Kd_medium = float(value)
                                elif name == 'Kd_strong':
                                    self.Kd_strong = float(value)
                                elif name == 'k_on':
                                    self.k_on = float(value)
                                elif name == 'ResizeVolume':
                                    self.ResizeVolume = float(value)
                                elif name == 'RescaleTemp':
                                    self.RescaleTemperature = float(value)
                                elif name == 'RingClosureFactor':
                                    self.RingClosureFactor = float(value)
                                elif name == 'initial_mixture':
                                    ka.system.mixture_file = value.strip()
                                elif name == 'reproducible':
                                    if "True" in value or 'true' in value:
                                        ka.system.monitor.reproducible = True
                                    else:
                                        ka.system.monitor.reproducible = False
                                elif name == 'canonicalize':
                                    if "True" in value or 'true' in value:
                                        ka.system.canonicalize = True
                                    else:
                                        ka.system.canonicalize = False
                                elif name == 'consolidate':
                                    if "True" in value or 'true' in value:
                                        ka.system.consolidate = True
                                    else:
                                        ka.system.consolidate = False
                                elif name == 'barcode':
                                    if "True" in value or 'true' in value:
                                        ka.system.barcode = True
                                    else:
                                        ka.system.barcode = False
                                elif name == 'sim_limit':
                                    match = re.match(r'%par:\s*sim_limit\s*=\s*(\S*)\s*(\S*)\s*', line)
                                    ka.system.sim_limit = float(match.group(1))
                                    ka.system.sim_limit_type = match.group(2)
                                elif name == 'obs_frequency':
                                    ka.system.monitor.obs_period = float(value)
                                elif name == 'snap_frequency':
                                    ka.system.monitor.snap_period = float(value)
                                elif name == 'seed':
                                    if value != 'None':
                                        self.rng_seed = int(value)
                                elif name == "memory":
                                    ka.system.monitor.memory = int(value)
                                elif name == 'inflow':
                                    match = re.match(r'%par:\s*inflow\s*=\s*(\S*)\s*(\S*)\s*', line)
                                    self.inflow[match.group(2)] = float(match.group(1))
                                elif name == 'outflow':
                                    match = re.match(r'%par:\s*outflow\s*=\s*(\S*)\s*(\S*)\s*', line)
                                    self.outflow[match.group(2)] = float(match.group(1))
                                else:
                                    sys.exit(f"unknown parameter file keyword in {line}")
                        elif match.group(1) == 'sig:':
                            match = re.match(r'%sig: ([^/]*)/?', line)
                            if match:
                                self.signature_string = match.group(1).strip()
                        elif match.group(1) == 'rep:':
                            match = re.match(r'%rep:\s*(.*)\s*=\s*(\S*)\s?', line)
                            if match:
                                name = match.group(1).strip()
                                value = match.group(2).strip()
                                if name == 'report_fn':
                                    ka.system.report_file = value
                                if name == 'output_fn':
                                    ka.system.monitor.obs_file_name = value
                                if name == 'snap_root':
                                    ka.system.monitor.snap_root_name = value
                                if name == 'numbering':
                                    ka.system.monitor.snap_numbering = value

    def report(self, pp_width=40):
        """
        Pretty print the system parameters
        """
        form = '1.5E'
        info = f"\n{'PARAMETERS '.ljust(70, '-')}\n\n"

        # now = datetime.datetime.now()
        # info += f'{"date":>30}: {now:%Y-%m-%d %H:%M}\n'
        # info += f'{"uuid":>30}: {self.uuid}\n'
        #
        # info += '\n'

        info += f'{"reference Vol":>{pp_width}}: {self.referenceVol:{form}} L\n'
        info += f'{"reference Temp":>{pp_width}}: {self.referenceTemp:{form}} K\n'
        info += f'{"Volume":>{pp_width}}: {self.Volume:{form}} L\n'
        info += f'{"Temperature":>{pp_width}}: {self.Temperature:{form}} K ({self.Temperature - 273.15:.3f} ºC)\n'

        info += f'{"Kd weak":>{pp_width}}: {self.Kd_weak}\n'
        info += f'{"Kd medium":>{pp_width}}: {self.Kd_medium}\n'
        info += f'{"Kd strong":>{pp_width}}: {self.Kd_strong}\n'
        info += f'{"k_on":>{pp_width}}: {self.k_on:{form}}\n'

        info += '\n'

        info += f'{"ResizeVolume":>{pp_width}}: {self.ResizeVolume}\n'
        info += f'{"RescaleTemperature":>{pp_width}}: {self.RescaleTemperature}\n'
        info += f'{"RingClosureFactor (adjusted)":>{pp_width}}: {self.RingClosureFactor:{form}}\n'

        # stochastic rate constants-----------------------------------------------

        info += '\n'

        info += f'{"inter-molecular on-rate (s_on)":>{pp_width}}: {self.s_on:{form}}\n'
        info += f'{"intra-molecular on-rate (s_ring_on)":>{pp_width}}: {self.s_ring_on:{form}}\n'
        for bt in ka.system.rc_bond_dissociation:
            text = f'off-rate ({bt[0]}--{bt[1]})'
            info += f'{text:>{pp_width}}: {ka.system.rc_bond_dissociation[bt]:{form}}\n'

        info += '\n'

        if ka.system.inflow_rate or ka.system.outflow_rate:
            for at in ka.system.inflow_rate:
                text = f'inflow rate ({at})'
                info += f'{text:>{pp_width}}: {ka.system.inflow_rate[at]:{form}}\n'
            for at in ka.system.outflow_rate:
                text = f'outflow rate ({at})'
                info += f'{text:>{pp_width}}: {ka.system.outflow_rate[at]:{form}}\n'
            info += '\n'

        info += f'{"defaults":>{pp_width}}\n'
        info += f'{"s_off_weak":>{pp_width}}: {self.s_off_weak:{form}}\n'
        info += f'{"s_off_medium":>{pp_width}}: {self.s_off_medium:{form}}\n'
        info += f'{"s_off_strong":>{pp_width}}: {self.s_off_strong:{form}}\n'

        info += '\n'

        for a in self.init_agents:
            text = f'initial agents {a}'
            if ka.system.signature.init_agents[a] == '*':
                molar = str(self.default_concentration)
            else:
                molar = ka.system.signature.init_agents[a]
            info += f'{text:>{pp_width}}: {self.init_agents[a]} ({molar} nM)\n'

        info += '\n'

        info += f'{"random number seed":>{pp_width}}: {self.rng_seed}\n'

        return info

    def __str__(self, pp_width=40):
        return self.report(pp_width=pp_width)

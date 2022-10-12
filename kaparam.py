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
    def __init__(self, file=None, pdict=None):

        # constants

        # self.Avogadro = constants.Avogadro  # temporarily avoid scipy for Pythonista's sake...
        self.Avogadro = 6.02214e+23
        self.Volume_choices = {'fibro': 2.25e-12}
        self.init_default = 100   # default initial agent concentration in nM

        # default parameters------------------------------------------------------

        self.Kd_weak   = 1000.e-9  #   1 mumol
        self.Kd_medium =  100.e-9  # 100 nmol
        self.Kd_strong =   10.e-9  #  10 nmol

        self.k_on = 1.e+9  # (M s)^-1; diffusion-controlled limit

        self.Volume = self.Volume_choices['fibro']
        self.RingClosureFactor = 1.e+5  # ratio of binary Kd to unary Kd
        self.Resize = 1.

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

        self.rng_seed = 4711  # Eau de Cologne

        # signature string -------------------------------------------------------

        self.signature_string = None

        # update and finalize the parameters
        if file:
            self.read_parameters(file)
        elif pdict:
            self.set_parameters(pdict)

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

        self.update_parameters()

    def update_parameters(self):
        """
        Apply the system (volume) scale factor to the parameters.
        """
        self.Volume = self.Volume * self.Resize
        self.RingClosureFactor = self.RingClosureFactor * self.Resize

        # stochastic rate constants-----------------------------------------------

        # bi-molecular on-constant  (inter-binding)
        self.s_on = self.k_on / (self.Avogadro * self.Volume)
        # uni-molecular on-constant (intra-binding)
        self.s_ring_on = self.RingClosureFactor * self.s_on

        # off-constants are the same regardless of whether they lead to a fission
        self.s_off_weak = self.k_on * self.Kd_weak
        self.s_off_medium = self.k_on * self.Kd_medium
        self.s_off_strong = self.k_on * self.Kd_strong

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
                ka.system.rc_bond_dissociation[bt] = self.k_on * float(ka.system.signature.bond_types[bt]) * 1.e-9
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
            if init == '*':  # default abundance 100 nM
                self.init_agents[a] = int(self.init_default * 1.e-9 * self.Avogadro * self.Volume)
            else:
                self.init_agents[a] = int(int(init) * 1.e-9 * self.Avogadro * self.Volume)

    def set_parameters(self, pdict):  # needs updating; currently not used
        """
        Set parameters using a dictionary of {keyword: value}.
        It calls update_parameters() and can be used for resetting parameters.
        """
        for keyword in pdict:
            if keyword == 'Volume':
                if is_number(pdict['Volume']):
                    self.Volume = pdict['Volume']
                else:
                    if pdict['Volume'] in self.Volume_choices:
                        self.Volume = self.Volume_choices[pdict['Volume']]
                    else:
                        sys.exit(f"No such volume choice: {pdict['Volume']}")
            elif keyword == 'Kd_weak':
                self.Kd_weak = pdict['Kd_weak']
            elif keyword == 'Kd_medium':
                self.Kd_medium = pdict['Kd_medium']
            elif keyword == 'Kd_strong':
                self.Kd_strong = pdict['Kd_strong']
            elif keyword == 'k_on':
                self.k_on = pdict['k_on']
            elif keyword == 'Resize':
                self.Resize = pdict['Resize']
            elif keyword == 'RingClosureFactor':
                self.RingClosureFactor = pdict['RingClosureFactor']
            elif keyword == 'seed':
                self.rng_seed = pdict['seed']
            elif keyword == 'Signature':
                self.signature_string = pdict['Signature']
            else:
                sys.exit(f"Parameters.set_parameters(): No such keyword: {keyword}")

    def read_parameters(self, par_file):
        """
        Reads parameters from a file.

        par_file: file name containing parameter values
        Returns: nothing

        Current declaration keywords:
            %par
            %sig
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
                    match = re.match(r'^%(par:|sig:|obs:|rep:)\s?', line)
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
                                elif name == 'Kd_weak':
                                    self.Kd_weak = float(value)
                                elif name == 'Kd_medium':
                                    self.Kd_medium = float(value)
                                elif name == 'Kd_strong':
                                    self.Kd_strong = float(value)
                                elif name == 'k_on':
                                    self.k_on = float(value)
                                elif name == 'Resize':
                                    self.Resize = float(value)
                                elif name == 'RingClosureFactor':
                                    self.RingClosureFactor = float(value)
                                elif name == 'sim_limit':
                                    match = re.match(r'%par:\s*sim_limit\s*=\s*(\S*)\s*(\S*)\s*', line)
                                    ka.system.sim_limit = float(match.group(1))
                                    ka.system.sim_limit_type = match.group(2)
                                elif name == 'obs_frequency':
                                    ka.system.obs_freq = float(value)
                                elif name == 'seed':
                                    self.rng_seed = int(value)
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
                        elif match.group(1) == 'obs:':
                            match = re.match(r'%obs: \s*(\S*)\s*([^/]*)/?', line)
                            ka.system.observable[match.group(1)] += [match.group(2).strip()]
                        elif match.group(1) == 'rep:':
                            match = re.match(r'%rep:\s*(.*)\s*=\s*(\S*)\s?', line)
                            if match:
                                name = match.group(1).strip()
                                value = match.group(2).strip()
                                if name == 'report_fn':
                                    ka.system.report_file = value
                                if name == 'output_fn':
                                    ka.system.obs_file = value
                                if name == 'snap_root':
                                    ka.system.snap_root = value
                                if name == 'numbering':
                                    ka.system.snap_numbering = value

    def __str__(self, pp_width=40):
        """
        Pretty print the system parameters
        """
        info = f"\n{'PARAMETERS '.ljust(70, '-')}\n\n"

        # now = datetime.datetime.now()
        # info += f'{"date":>30}: {now:%Y-%m-%d %H:%M}\n'
        # info += f'{"uuid":>30}: {self.uuid}\n'
        #
        # info += '\n'

        info += f'{"Avogadro":>{pp_width}}: {self.Avogadro:1.5E}\n'
        info += f'{"Kd weak":>{pp_width}}: {self.Kd_weak}\n'
        info += f'{"Kd medium":>{pp_width}}: {self.Kd_medium}\n'
        info += f'{"Kd strong":>{pp_width}}: {self.Kd_strong}\n'
        info += f'{"k_on":>{pp_width}}: {self.k_on:1.3E}\n'

        info += '\n'

        info += f'{"Resize":>{pp_width}}: {self.Resize}\n'
        info += f'{"Volume (resized)":>{pp_width}}: {self.Volume:1.3E}\n'
        info += f'{"RingClosureFactor (resized)":>{pp_width}}: {self.RingClosureFactor:1.3E}\n'

        # stochastic rate constants-----------------------------------------------

        info += '\n'

        info += f'{"inter-molecular on-rate (s_on)":>{pp_width}}: {self.s_on:1.5E}\n'
        info += f'{"intra-molecular on-rate (s_ring_on)":>{pp_width}}: {self.s_ring_on:1.5E}\n'
        for bt in ka.system.rc_bond_dissociation:
            text = f'off-rate ({bt[0]}--{bt[1]})'
            info += f'{text:>{pp_width}}: {ka.system.rc_bond_dissociation[bt]:1.5E}\n'

        info += '\n'

        if ka.system.inflow_rate or ka.system.outflow_rate:
            for at in ka.system.inflow_rate:
                text = f'inflow rate ({at})'
                info += f'{text:>{pp_width}}: {ka.system.inflow_rate[at]:1.5E}\n'
            for at in ka.system.outflow_rate:
                text = f'outflow rate ({at})'
                info += f'{text:>{pp_width}}: {ka.system.outflow_rate[at]:1.5E}\n'
            info += '\n'

        info += f'{"defaults":>{pp_width}}\n'
        info += f'{"s_off_weak":>{pp_width}}: {self.s_off_weak:1.5E}\n'
        info += f'{"s_off_medium":>{pp_width}}: {self.s_off_medium:1.5E}\n'
        info += f'{"s_off_strong":>{pp_width}}: {self.s_off_strong:1.5E}\n'

        info += '\n'

        for a in self.init_agents:
            text = f'initial agents {a}'
            if ka.system.signature.init_agents[a] == '*':
                molar = str(self.init_default)
            else:
                molar = ka.system.signature.init_agents[a]
            info += f'{text:>{pp_width}}: {self.init_agents[a]} ({molar} nM)\n'

        info += '\n'

        info += f'{"random number seed":>{pp_width}}: {self.rng_seed}\n'

        return info

# -----------------------------------------------------------------


if __name__ == '__main__':
    import kainit

    kainit.minit(signature=None, parameter_file='TestData/parameters.txt')
    print(str(ka.system.parameters))


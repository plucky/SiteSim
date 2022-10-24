# Walter Fontana at 10/11/22

import math
import re
import sys
from collections import defaultdict
import kamol
import kasystem as ka


class Monitor:
    def __init__(self):
        self.monitor_size_re = re.compile(r'\s*size\s*\[(\d*)\s*-\s*(\d*)\]')
        self.monitor_max_size_re = re.compile(r'\s*maxsize\s*\[(\d*)\]')

        self.obs_period = 0.
        # observable types: ! -> molecule, ? -> pattern, b -> bond, s -> free site, p -> property
        self.observable = {'!': [], '?': [], 'b': [], 's': [], 'mb': [], 'ms': [], 'p': []}

        self.obs_file_name = ''            # observation file (typically a csv)
        self.snap_root_name = ''           # root fn of snapshots
        self.snap_numbering = 'serial'     # snapshot numbering scheme: {serial, event}
        self.snap_period = 0.
        self.snap_counter = 0
        self.reproducible = 'False'

        self.observation_time = 0
        self.snap_time = 0

        self.min_size = 0
        self.max_size = 0
        self.max_size_ranks = 0
        self.name_form = ''

    def initialize(self):
        """
        Convert molecular observables into canonical form for fast retrieval. (Requires the mixture-wide
        local views.) Set up other observables. Write the column labels of the data file.
        """

        if ka.system.mixture_file and self.reproducible:
            rev = ka.system.mixture_file[::-1]   # reverse a string
            rev = rev[3:]  # get rid of ak.
            n = ''
            i = 0
            while rev[i:i+1].isdigit():
                n += rev[i:i+1]
                i += 1
            self.snap_counter = int(n[::-1])  # reverse again

        # set up file name format for snapshots
        if self.snap_numbering == 'serial':
            # compute the field size to pad numbers for files to be numerically sorted by OS
            if ka.system.sim_limit > 0:
                field_size = int(math.log10(ka.system.sim_limit / self.snap_period) + 1)
                self.name_form = '{' + f':0{field_size}d' + '}'

        if ka.system.sim_limit_type == 'time':
            self.observation_time = ka.system.sim.time
            self.snap_time = ka.system.sim.time
            info = f'time,'
        else:
            self.observation_time = ka.system.sim.event
            self.snap_time = ka.system.sim.event
            info = f'event,'

        # set up molecule observables
        internal = []
        for m_expr in self.observable['!']:
            # this also generates canonical forms for the molecular observables
            m = kamol.KappaComplex(m_expr, system=ka.system)
            internal.append(m)
            txt = re.sub(r'\),', ')', m.kappa_expression())
            info += f"{txt.strip()},"
        # replace observables with their internal representation
        self.observable['!'] = internal

        # set up pattern observables
        internal = []
        for pattern in self.observable['?']:
            m = kamol.KappaComplex(pattern, system=None, canon=False)
            internal.append(m)
            txt = re.sub(r'\),', ')', m.kappa_expression())
            info += f'?{txt.strip()},'
        self.observable['?'] = internal

        # set up bond type observables
        internal = []
        for item in self.observable['b']:
            temp = item.split('-')
            internal.append(tuple(sorted(temp)))
            info += f'{item},'
        self.observable['b'] = internal

        # site type observables
        internal = []
        for item in self.observable['s']:
            info += f'{item},'
        self.observable['s'] = internal

        # set up bond type observables for maximer
        internal = []
        for item in self.observable['mb']:
            temp = item.split('-')
            internal.append(tuple(sorted(temp)))
            info += f'm:{item},'
        self.observable['mb'] = internal

        # site type observables for maximer
        internal = []
        for item in self.observable['ms']:
            info += f'm:{item},'
        self.observable['ms'] = internal

        # size distribution observables
        for item in self.observable['p']:
            if 'size' in item:  # size distribution
                if 'maxsize' in item:
                    # syntax is
                    # %obs: p maxsize [n]
                    # where n is the number of the largest n complexes
                    match = self.monitor_max_size_re.match(item)
                    self.max_size_ranks = int(match.group(1))
                    for i in range(1, self.max_size_ranks + 1):
                        info += f'sz-rank {i},'
                else:
                    # syntax is
                    # %obs: p size [min-max]
                    # where min, max define the size range
                    match = self.monitor_size_re.match(item)
                    if not match:
                        sys.exit("could not parse request for size observation")
                    else:
                        self.min_size = int(match.group(1))
                        self.max_size = int(match.group(2))
                        for i in range(self.min_size, self.max_size+1):
                            info += f'size {i},'

        info = info[:-1]

        # initialize file with column labels
        with open(self.obs_file_name, "w") as fp:
            fp.write(info + '\n')

    def observe(self):
        """
        Execute the requested observations and add to the data file.
        """
        if ka.system.sim_limit_type == 'time':
            info = f'{ka.system.sim.time}, '
        else:
            info = f'{ka.system.sim.event}, '

        for m in self.observable['!']:
            if m.canonical in ka.system.mixture.canonical:
                molecule = ka.system.mixture.canonical[m.canonical]
                info += f'{molecule.count}, '
            else:
                info += '0, '

        for item in self.observable['?']:
            embeddings = 0
            for m in ka.system.mixture.complexes:
                embeddings += ka.system.sgm.number_of_all_embeddings(m, item) * m.count
            info += f'{embeddings}, '

        for item in self.observable['b']:
            info += f'{ka.system.mixture.total_bond_type[item]}, '

        for item in self.observable['s']:
            info += f'{ka.system.mixture.total_free_sites[item]}, '

        sorted_complexes = []
        if self.observable['p'] or self.observable['mb'] or self.observable['ms']:
            sorted_complexes = sorted(ka.system.mixture.complexes, key=lambda x: x.size, reverse=True)

        for item in self.observable['mb']:
            info += f'{sorted_complexes[0].bond_type[item]}, '

        for item in self.observable['ms']:
            info += f'{sorted_complexes[0].free_site[item]}, '

        for item in self.observable['p']:
            if 'size' in item:  # size distribution
                if 'maxsize' in item:
                    # determine max size count
                    i = 1
                    for m in sorted_complexes:
                        info += f'{m.size}, '
                        i += 1
                        if i > self.max_size_ranks:
                            break
                else:
                    obs = defaultdict(int)
                    for m in ka.system.mixture.complexes:
                        if self.min_size <= m.size <= self.max_size:
                            obs[m.size] += m.count
                    for i in range(self.min_size, self.max_size+1):
                        info += f'{obs[i]}, '

        info = info[:-2]
        with open(self.obs_file_name, "a") as fp:
            fp.write(info + '\n')

        self.observation_time += self.obs_period

    def snapshot(self):
        """
        Make a snapshot of the mixture.
        """
        # assemble filename
        if self.name_form:
            snap_fn = self.snap_root_name + self.name_form.format(self.snap_counter) + '.ka'
        else:
            snap_fn = self.snap_root_name + f'{self.snap_counter}' + '.ka'

        ka.system.mixture.make_snapshot(snap_fn)
        self.snap_counter += 1

        self.snap_time += self.snap_period


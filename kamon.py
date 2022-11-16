# Walter Fontana, 2022

import math
import re
import sys
from collections import defaultdict
import kamol
import kasystem as ka


class Monitor:
    """
    --> %obs: ! <Kappa expression>              // # of instances of fully specified molecule, given in Kappa syntax

        examples
        %obs: ! A(s[1]), B(s[2]), S(a[1] b[2])  // a molecule
        %obs: ! S(a[.] b[.])                    // another molecule

    --> %obs: ? <Kappa pattern>                 // # of instances of all molecules matching a pattern
    --> %obs: ? A(r[.]) size [min-max]]         // # of instances of pattern in size classes [min-max]

        examples
        %obs: ? S(a[#] b[#])                    // a molecular pattern ('#' means "don't care")
        %obs: ? S(a[.] b[_])                    // a pattern ('_' means 'bound')
        %obs: ? S(a[.] b[1]), B(s[1])           // a pattern is relative to a signature ('don't care, don't mention')
        %obs: ? A(r[.]) size [1-20]             // pattern occurrence in size classes

    --> %obs: b <bond type>                     // total # of instance of a given bond type...
    --> %obs: mb <bond type>                    // total # of instance of a given bond type in maximer...
                                                // ...with bond type: <atom-type>.<site-name> - <atom-type>.<site-name>
        examples
        %obs: b A.l-A.r                         // a bond type
        %obs: b P.d-P.d                         // a bond type
        %obs: mb A.p-P.a1                       // a bond type in the current maximer

    --> %obs: s <site-type>                     // # of instances of a free site type...
    --> %obs: ms <site-type>                    // # of instances of a free site type in maximer...
                                                // ... with free site type: <atom-type>.<site-name>
        examples
        %obs: s A.l                             // a free site type
        %obs: ms P.d                            // a free site type in the current maximer

    --> %obs: p size [min-max]]                 // report the size distribution in the size range [min-max]
    --> %obs: p maxsize [n]                     // # of particles of each of the n largest molecules

        example
        %obs: p size [1-20]                     // # particles of sizes 1 to 20 (reported for each size class)
    """
    def __init__(self):
        self.monitor_size_re = re.compile(r'\s*size\s*\[(\d*)\s*-\s*(\d*)\]')
        self.monitor_range_re = re.compile(r'\s*\[(\d*)\s*-\s*(\d*)\]')
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
            rev = ka.system.mixture_file[::-1]   # reverse string
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
            if ka.system.sim_limit > 0 and self.snap_period > 0:
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

        if self.snap_period == 0:
            self.snap_time = ka.system.sim_limit

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
            interval = None
            if 'size' in pattern:
                pat = pattern.split('size')
                match = self.monitor_range_re.match(pat[1])
                interval = (int(match.group(1)), int(match.group(2)))
            m = kamol.KappaComplex(pat[0].strip(), system=None, canon=False)
            internal.append((m, interval))
            txt = re.sub(r'\),', ')', m.kappa_expression())
            if interval:
                for i in range(interval[0], interval[1] + 1):
                    info += f'?{txt.strip()} in size {i},'
            else:
                info += f'?{txt.strip()},'
        self.observable['?'] = internal

        # set up bond type observables
        internal = []
        for item in self.observable['b']:
            temp = item.split('-')
            internal.append(tuple(sorted(temp)))
            info += f'{item},'
        self.observable['b'] = internal

        # set up site type observables
        for item in self.observable['s']:
            info += f'{item},'

        # set up bond type observables for maximer
        internal = []
        for item in self.observable['mb']:
            temp = item.split('-')
            internal.append(tuple(sorted(temp)))
            info += f'mb {item},'
        self.observable['mb'] = internal

        # site type observables for maximer
        internal = []
        for item in self.observable['ms']:
            info += f'ms {item},'
        self.observable['ms'] = internal

        # size distribution observables
        for item in self.observable['p']:
            if 'size' in item:  # size distribution
                if 'maxsize' in item:
                    # syntax is
                    # %obs: p maxsize [n]
                    # where n denotes the largest n complexes
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
                        for i in range(self.min_size, self.max_size + 1):
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
            pattern = item[0]
            embeddings = defaultdict(int)
            if item[1]:
                # user specified a size range stratifying the embeddings of the pattern
                from_, to_ = item[1]
                for m in ka.system.mixture.complexes:
                    if from_ <= m.size <= to_:
                        embeddings[m.size] += ka.system.sgm.number_of_all_embeddings(m, pattern) * m.count
                for i in range(from_, to_ + 1):
                    info += f'{embeddings[i]}, '
            else:
                for m in ka.system.mixture.complexes:
                    embeddings[0] += ka.system.sgm.number_of_all_embeddings(m, pattern) * m.count
                info += f'{embeddings[0]}, '

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
                    # determine the largest 'max_size_ranks' sizes
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

    def snapshot(self, flag=''):
        """
        Make a snapshot of the mixture.
        """
        # assemble filename
        if not flag:
            if self.name_form:
                snap_fn = self.snap_root_name + self.name_form.format(self.snap_counter) + '.ka'
            else:
                snap_fn = self.snap_root_name + f'{self.snap_counter}' + '.ka'
        elif flag == 'first':
            snap_fn = self.snap_root_name + "_start.ka"
        else:
            snap_fn = self.snap_root_name + "_end.ka"

        ka.system.mixture.make_snapshot(snap_fn)
        self.snap_counter += 1

        self.snap_time += self.snap_period


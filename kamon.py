# Walter Fontana, 2022

import math
import re
import sys
import os
from collections import defaultdict
from collections import deque
import kamol
import kasystem as ka


class Monitor:
    """
    assigning a [name] to an observable is optional, but not permitted for p-type observables

    --> %obs: ([name]) ! <Kappa expression>     // # of instances of fully specified molecule, given in Kappa syntax

        examples
        %obs: ! A(s[1]), B(s[2]), S(a[1] b[2])  // a molecule
        %obs: [free_S] ! S(a[.] b[.])           // another molecule

    --> %obs: ([name]) ? <Kappa pattern>        // # of instances of all molecules matching a pattern
    --> %obs: ? A(r[.]) size [min-max]]         // # of instances of pattern in size classes [min-max]

        examples
        %obs: ? S(a[#] b[#])                    // a molecular pattern ('#' means "don't care")
        %obs: ? S(a[.] b[_])                    // a pattern ('_' means 'bound')
        %obs: [SB] ? S(a[.] b[1]), B(s[1])      // a pattern is relative to a signature ('don't care, don't mention')
        %obs: ? A(r[.]) size [1-20]             // pattern occurrence in size classes 1 to 20

    --> %obs: ([name]) b <bond type>            // total # of instance of a given bond type...
    --> %obs: ([name]) mb <bond type>           // total # of instance of a given bond type in maximer...
                                                // ...with bond type: <atom-type>.<site-name> - <atom-type>.<site-name>
        examples
        %obs: [A_dimer] b A.l-A.r               // a bond type
        %obs: b P.d-P.d                         // a bond type
        %obs: mb A.p-P.a1                       // a bond type in the current maximer

    --> %obs: ([name]) s <site-type>            // # of instances of a free site type...
    --> %obs: ([name]) ms <site-type>           // # of instances of a free site type in maximer...
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
        self.monitor_range_re = re.compile(r'\s*\[(\d*)\s*-\s*(\d*)\]')
        self.monitor_max_size_re = re.compile(r'\s*\[(\d*)\]')

        self.obs_period = 0.
        # observable types: ! -> molecule, ? -> pattern, b -> bond, s -> free site, p -> property
        self.observable = {'!': {}, '?': {}, 'b': {}, 's': {}, 'mb': {}, 'ms': {}, 'p': {}}
        self.observable_by_name = {}
        self.memory = 1                    # number of events to remember values of

        self.obs_file_name = ''            # observation file (typically a csv)
        self.snap_root_name = ''           # root fn of snapshots
        self.snap_numbering = 'serial'     # snapshot numbering scheme: {serial, event}
        self.snap_period = 0.
        self.snap_counter = 0
        self.reproducible = 'False'

        self.observation_time = 0
        self.snap_time = 0

        self.name_form = ''

    def parse_observables(self, file=None):
        """
        Reads the %obs portion from the input file and constructs observables
        observable[type] =
                        {
                        name:
                            {
                            'pattern': pattern,
                            'label': column label,
                            'spec': {'type': type, 'value': value}
                            'value': {}  # dictionary of deques to hold time series
                            }
                        }
        """
        if not file:
            sys.exit("No parameter file specified.")

        if not os.path.isfile(file):
            sys.exit("Cannot find parameter file %s" % file)
        else:
            name_list = {}
            default_name = 1
            with open(file, "r", encoding='utf-8') as data:
                while True:
                    line = data.readline()
                    if not line:
                        break

                    # parse the line
                    match = re.match(r'%obs: \s*((?:".+"))?\s*(\S*)\s*([^/]*)/?', line)

                    if match:
                        # group 1: optional name of observable
                        # group 2: observable type
                        # group 3: observable
                        obs_name = match.group(1)
                        if obs_name:
                            # sanity check that name is well-formed
                            if not re.match(r'".+"', obs_name):
                                sys.exit(f'invalid observable name: {obs_name}')
                            # sanity check that name is unique
                            if obs_name in name_list:
                                sys.exit(f'observable name is duplicate: {obs_name}')
                        else:
                            obs_name = '*' + str(default_name)
                        default_name += 1
                        name_list[obs_name] = 1

                        obs_type = match.group(2).strip()
                        item = match.group(3).strip()

                        self.observable[obs_type][obs_name] = {}
                        self.observable[obs_type][obs_name]['pattern'] = None
                        self.observable[obs_type][obs_name]['label'] = None
                        self.observable[obs_type][obs_name]['spec'] = None
                        self.observable[obs_type][obs_name]['value'] = None

                        if obs_name[0] == '*':
                            label = item
                        else:
                            label = obs_name

                        match obs_type:
                            # molecule observables
                            case '!':
                                m = kamol.KappaComplex(item, system=ka.system)
                                if obs_name[0] == '*':
                                    text = re.sub(r'\),', ')', m.kappa_expression())
                                    label = f'!{text.strip()}'
                                else:
                                    label = obs_name
                                #------------------------------------------------------------
                                self.observable[obs_type][obs_name]['pattern'] = m
                                self.observable[obs_type][obs_name]['label'] = [label]
                                self.observable[obs_type][obs_name]['value'] = {0: deque([])}
                                # ------------------------------------------------------------
                            # pattern observables
                            case '?':
                                if 'size' in item:
                                    pattern, size_range = item.split('size')
                                    m = kamol.KappaComplex(pattern, system=ka.system)
                                    result = self.monitor_range_re.match(size_range)
                                    if not result:
                                        sys.exit("could not parse size range")
                                    interval = (int(result.group(1)), int(result.group(2)))
                                    spec = {'type': 'size range', 'value': interval}
                                    if obs_name[0] == '*':
                                        text = re.sub(r'\),', ')', m.kappa_expression())
                                        label = f'!{text.strip()}'
                                    else:
                                        label = obs_name
                                    labels = []
                                    values = {}
                                    for i in range(interval[0], interval[1] + 1):
                                        labels += [f'?{label} in size {i}']
                                        values[i] = deque([])
                                    # ------------------------------------------------------------
                                    self.observable[obs_type][obs_name]['pattern'] = m
                                    self.observable[obs_type][obs_name]['label'] = labels
                                    self.observable[obs_type][obs_name]['spec'] = spec
                                    self.observable[obs_type][obs_name]['value'] = values
                                    # ------------------------------------------------------------
                                else:
                                    m = kamol.KappaComplex(item, system=ka.system)
                                    if obs_name[0] == '*':
                                        text = re.sub(r'\),', ')', m.kappa_expression())
                                        label = f'!{text.strip()}'
                                    else:
                                        label = obs_name
                                    # ------------------------------------------------------------
                                    self.observable[obs_type][obs_name]['pattern'] = m
                                    self.observable[obs_type][obs_name]['label'] = [f'?{label}']
                                    self.observable[obs_type][obs_name]['value'] = {0: deque([])}
                                    # ------------------------------------------------------------
                            # bond type observables
                            case 'b':
                                ports = item.split('-')
                                # ------------------------------------------------------------
                                self.observable[obs_type][obs_name]['pattern'] = tuple(sorted(ports))
                                self.observable[obs_type][obs_name]['label'] = [label]
                                self.observable[obs_type][obs_name]['value'] = {0: deque([])}
                                # ------------------------------------------------------------
                            # site type observables
                            case 's':
                                # ------------------------------------------------------------
                                self.observable[obs_type][obs_name]['pattern'] = item
                                self.observable[obs_type][obs_name]['label'] = [label]
                                self.observable[obs_type][obs_name]['value'] = {0: deque([])}
                                # ------------------------------------------------------------
                            # bond type observables in maximer
                            case 'mb':
                                ports = item.split('-')
                                # ------------------------------------------------------------
                                self.observable[obs_type][obs_name]['pattern'] = tuple(sorted(ports))
                                self.observable[obs_type][obs_name]['label'] = [f'mb {label}']
                                self.observable[obs_type][obs_name]['value'] = {0: deque([])}
                                # ------------------------------------------------------------
                            # site type observables in maximer
                            case 'ms':
                                # ------------------------------------------------------------
                                self.observable[obs_type][obs_name]['pattern'] = item
                                self.observable[obs_type][obs_name]['label'] = [f'ms {label}']
                                self.observable[obs_type][obs_name]['value'] = {0: deque([])}
                                # ------------------------------------------------------------
                            # special (property) observables
                            case 'p':
                                if 'size' in item:  # size distribution
                                    if 'maxsize' in item:
                                        # syntax is /%obs: p maxsize [n]/, where n denotes the largest n complexes
                                        pattern, size_range = item.split('maxsize')
                                        result = self.monitor_max_size_re.match(size_range.strip())
                                        ranks = int(result.group(1))
                                        spec = {'type': 'top sizes', 'value': ranks}
                                        labels = []
                                        values = {}
                                        for i in range(1, ranks + 1):
                                            labels += [f'sz-rank {i}']
                                            values[i] = deque([])
                                        # ------------------------------------------------------------
                                        self.observable[obs_type][obs_name]['label'] = labels
                                        self.observable[obs_type][obs_name]['spec'] = spec
                                        self.observable[obs_type][obs_name]['value'] = values
                                        # ------------------------------------------------------------
                                    else:
                                        # syntax is /%obs: p size [min-max]/, where min, max define the size range
                                        pattern, size_range = item.split('size')
                                        result = self.monitor_range_re.match(size_range.strip())
                                        if not result:
                                            sys.exit("could not parse size range")
                                        interval = (int(result.group(1)), int(result.group(2)))
                                        spec = {'type': 'size range', 'value': interval}
                                        labels = []
                                        values = {}
                                        for i in range(interval[0], interval[1] + 1):
                                            labels += [f'size {i}']
                                            values[i] = deque([])
                                        # ------------------------------------------------------------
                                        self.observable[obs_type][obs_name]['label'] = labels
                                        self.observable[obs_type][obs_name]['spec'] = spec
                                        self.observable[obs_type][obs_name]['value'] = values
                                        # ------------------------------------------------------------

    def initialize(self, file=None):
        """
        Convert molecular observables to canonical form for fast retrieval. (Requires the mixture-wide
        local views.) Set up all other observables. Write the column labels of the data file.
        """

        self.parse_observables(file=file)

        # get the latest snapshot file count suffix
        if ka.system.mixture_file and self.reproducible:
            rev = ka.system.mixture_file[::-1]   # reverse string
            rev = rev[3:]  # get rid of ak.
            n = ''
            i = 0
            while rev[i:i+1].isdigit():
                n += rev[i:i+1]
                i += 1
            if n != '':
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
            info = f'time, '
        else:
            self.observation_time = ka.system.sim.event
            self.snap_time = ka.system.sim.event
            info = f'event, '

        if self.snap_period == 0:
            self.snap_time = ka.system.sim_limit

        # generate the column labels for the monitor file and
        # a "by name" indexed copy of the dictionary of observables
        for obs_type in self.observable:
            for name in self.observable[obs_type]:
                self.observable_by_name[name] = {obs_type: self.observable[obs_type][name]}
                for item in self.observable[obs_type][name]['label']:
                    info += f'{item}, '
        info = info[:-2]
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

        sorted_complexes = []
        if self.observable['p'] or self.observable['mb'] or self.observable['ms']:
            sorted_complexes = sorted(ka.system.mixture.complexes, key=lambda x: x.size, reverse=True)

        for obs_type in self.observable:
            match obs_type:
                case '!':
                    for name in self.observable[obs_type]:
                        obs = self.observable[obs_type][name]
                        m = obs['pattern']
                        if m.canonical in ka.system.mixture.canonical:
                            molecule = ka.system.mixture.canonical[m.canonical]
                            value = molecule.count
                        else:
                            value = 0
                        obs['value'][0].append(value)
                        if len(obs['value'][0]) > self.memory:
                            obs['value'][0].popleft()
                        info += f'{value}, '
                case '?':
                    for name in self.observable[obs_type]:
                        obs = self.observable[obs_type][name]
                        pattern = obs['pattern']
                        spec = obs['spec']
                        if spec:
                            if spec['type'] == 'size range':
                                # user specified a size range stratifying the embeddings of the pattern
                                from_, to_ = spec['value']
                                embed = defaultdict(int)
                                for m in ka.system.mixture.complexes:
                                    if from_ <= m.size <= to_:
                                        embed[m.size] += ka.system.sgm.number_of_all_embeddings(m, pattern) * m.count
                                for i in range(from_, to_ + 1):
                                    obs['value'][i].append(embed[i])
                                    if len(obs['value'][i]) > self.memory:
                                        obs['value'][i].popleft()
                                    info += f'{embed[i]}, '
                        else:
                            embed = 0
                            for m in ka.system.mixture.complexes:
                                embed += ka.system.sgm.number_of_all_embeddings(m, pattern) * m.count
                            obs['value'][0].append(embed)
                            if len(obs['value'][0]) > self.memory:
                                obs['value'][0].popleft()
                            info += f'{embed}, '
                case 'b':
                    for name in self.observable[obs_type]:
                        obs = self.observable[obs_type][name]
                        bt = obs['pattern']
                        value = ka.system.mixture.total_bond_type[bt]
                        obs['value'][0].append(value)
                        if len(obs['value'][0]) > self.memory:
                            obs['value'][0].popleft()
                        info += f'{value}, '
                case 's':
                    for name in self.observable[obs_type]:
                        obs = self.observable[obs_type][name]
                        st = obs['pattern']
                        value = ka.system.mixture.total_free_sites[st]
                        obs['value'][0].append(value)
                        if len(obs['value'][0]) > self.memory:
                            obs['value'][0].popleft()
                        info += f'{value}, '
                case 'mb':
                    for name in self.observable[obs_type]:
                        obs = self.observable[obs_type][name]
                        bt = obs['pattern']
                        value = sorted_complexes[0].bond_type[bt]
                        obs['value'][0].append(value)
                        if len(obs['value'][0]) > self.memory:
                            obs['value'][0].popleft()
                        info += f'{value}, '
                case 'ms':
                    for name in self.observable[obs_type]:
                        obs = self.observable[obs_type][name]
                        st = obs['pattern']
                        value = sorted_complexes[0].free_site[st]
                        obs['value'][0].append(value)
                        if len(obs['value'][0]) > self.memory:
                            obs['value'][0].popleft()
                        info += f'{value}, '
                case 'p':
                    for name in self.observable[obs_type]:
                        obs = self.observable[obs_type][name]
                        spec = obs['spec']
                        if spec['type'] == 'top sizes':
                            # determine the largest 'size_ranks' sizes
                            i = 1
                            for m in sorted_complexes:
                                obs['value'][i].append(m.size)
                                if len(obs['value'][i]) > self.memory:
                                    obs['value'][i].popleft()
                                info += f'{m.size}, '
                                i += 1
                                if i > spec['value']:
                                    break
                        elif spec['type'] == 'size range':
                            counts = defaultdict(int)
                            min_size, max_size = spec['value']
                            for m in ka.system.mixture.complexes:
                                if min_size <= m.size <= max_size:
                                    counts[m.size] += m.count
                            for i in range(min_size, max_size + 1):
                                obs['value'][i].append(counts[i])
                                if len(obs['value'][i]) > self.memory:
                                    obs['value'][i].popleft()
                                info += f'{counts[i]}, '

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


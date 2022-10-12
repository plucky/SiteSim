# Walter Fontana at 10/11/22

import math
import re
import kamol
import kasystem as ka


def make_snapshot_filename():
    """
    Construct a snapshot file name.
    """
    if ka.system.snap_numbering == 'serial':
        # compute the field size to pad numbers for files to be numerically sorted by OS
        field_size = int(math.log10(ka.system.sim_limit / ka.system.obs_freq) + 1)
        form = '{' + f':0{field_size}d' + '}'
        snap_fn = ka.system.snap_root + form.format(ka.system.obs_counter) + '.ka'
    else:
        snap_fn = ka.system.snap_root + f'{ka.system.obs_counter}' + '.ka'
    return snap_fn


def make_snapshot():
    pass


def set_up_monitoring():
    """
    Convert molecular observables into canonical form for fast retrieval.
    Requires the mixture-wide local views.
    """
    internal = []
    info = 'time,'
    for m_expr in ka.system.observable['!']:
        # this also generates canonical forms for the molecular observables
        m = kamol.KappaComplex(m_expr, system=ka.system)
        internal.append(m)
        txt = re.sub(r'\),', ')', m.kappa_expression())
        info += f"{txt.strip()},"
    # replace observables with their internal representation
    ka.system.observable['!'] = internal

    internal = []
    for pattern in ka.system.observable['?']:
        m = kamol.KappaComplex(pattern, system=None, canon=False)
        internal.append(m)
        txt = re.sub(r'\),', ')', m.kappa_expression())
        info += f'?{txt.strip()},'
    ka.system.observable['?'] = internal

    internal = []
    for item in ka.system.observable['b']:
        temp = item.split('-')
        internal.append(tuple(sorted(temp)))
        info += f'{item},'
    ka.system.observable['b'] = internal

    for item in ka.system.observable['s']:
        info += f'{item},'
    info = info[:-1]

    # initialize file with contents string
    with open(ka.system.obs_file, "w") as fp:
        fp.write(info + '\n')


def monitor():
    info = f'{ka.system.sim.time}, '
    for m in ka.system.observable['!']:
        if m.canonical in ka.system.mixture.canonical:
            molecule = ka.system.mixture.canonical[m.canonical]
            info += f'{molecule.count}, '
        else:
            info += '0, '

    for item in ka.system.observable['?']:
        embeddings = 0
        for m in ka.system.mixture.complexes:
            embeddings += ka.system.sgm.number_of_all_embeddings(m, item) * m.count
        info += f'{embeddings}, '
    for item in ka.system.observable['b']:
        info += f'{ka.system.mixture.total_bond_type[item]}, '
    for item in ka.system.observable['s']:
        info += f'{ka.system.mixture.total_free_sites[item]}, '
    info = info[:-2]
    with open(ka.system.obs_file, "a") as fp:
        fp.write(info + '\n')



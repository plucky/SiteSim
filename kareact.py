# Walter Fontana 2022

from collections import deque
import re

import kamol
import kasystem as ka


def convert(text): return int(text) if text.isdigit() else text


def alphanum_key(key): return [convert(c) for c in re.split('([0-9]+)', key)]


def components(kappaMol, traverse='bfs'):
    """
    Determines the graphical components of a KappaMolecule after bond dissociation.
    Returns an agent dictionary for each component (as requisite for creating a new KappaMolecule).

    This is standard bfs or dfs graph traversal. We use queues rather than recursion,
    because recursion seems very slow in Python.
    """
    visited = set()
    agents = {}
    nc = 0

    # BFS with queue
    if traverse == 'bfs':
        queue = deque()  # [1st, 2nd, 3rd, 4th, ...]  add via append(); remove bia popleft()
        for node in kappaMol.agents:
            if node not in visited:
                # if nc == 1 at this point, we can exit, since dissociation
                # cannot generate more than 2 components. (Check if this is efficient.)
                if nc == 1:
                    break
                nc += 1
                queue.append(node)
                visited.add(node)
                while queue:
                    current = queue.popleft()  # a node in this component
                    agents[current] = kappaMol.agents[current]
                    for neighbor in kappaMol.adjacency[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
    # DFS with stack
    elif traverse == 'dfs':
        stack = deque()  # [1st, 2nd, 3rd, 4th, ...]  add via append(); remove via pop()
        for node in kappaMol.agents:
            if node not in visited:
                if nc == 1:
                    break
                nc += 1
                stack.append(node)
                while stack:
                    current = stack.pop()  # a node in this component
                    if current not in visited:
                        visited.add(current)
                        agents[current] = kappaMol.agents[current]
                        for neighbor in kappaMol.adjacency[current]:
                            if neighbor not in visited:
                                stack.append(neighbor)
    # expressions = []
    # for i in range(0, nc):
    #     expressions += [kappa_expression(agents[i], kappaMol.bonds, kappaMol.bondsep)]

    if len(agents) == len(kappaMol.agents):
        return [kappaMol.agents]
    else:
        list(map(lambda k: kappaMol.agents.pop(k, None), agents))
        return [agents, kappaMol.agents]  # agent dictionaries, one for each component


def local_view(name, agents, bond_sep='@'):
    """
    Determine the local view of 'name' in 'agents'.
    """
    lv = []
    iface = agents[name]['iface']
    for s in iface:
        view = ''
        b = iface[s]['bond']
        if b != '.' and b != '#':
            other_name, other_s = b.split(bond_sep)
            other_type = agents[other_name]['info']['type']
            view += f'[{other_type}.{other_s}]'
        else:
            view += f'[{b}]'
        # skip the state in this specific context
        # view += '{' + f"{iface[s]['state']}" + '}'
        lv.append((s, view))
    local = ''
    for site_view in [f"{s}{view} " for (s, view) in sorted(lv)]:
        local += site_view
    # this is the local view of agent 'name'
    loc_view = agents[name]['info']['type'] + '(' + local[:-1] + ')'

    # update the local_views at the systems level, if needed
    mix = ka.system.mixture
    if loc_view not in mix.local_views:
        running_id = mix.local_views[next(reversed(mix.local_views))]
        mix.local_views[loc_view] = running_id + 1

    return loc_view


def dissociate_bond(A, port1, port2):
    """
    Remove the bond (port1, port2) within KappaMolecule 'A'. (port1, port2) is a standardized representation
    of the link between two agents at specific sites.
    """
    # (port1, port2) is standardized
    (a_agent, a_site) = port1
    (b_agent, b_site) = port2
    # break the bond
    A.agents[a_agent]['iface'][a_site]['bond'] = '.'
    A.agents[b_agent]['iface'][b_site]['bond'] = '.'
    # update adjacency
    # Using remove, eg A.adjacency[a_agent].remove(b_agent), is correct, but can result in a different
    # list sequence compared with a de novo construction, such as when reading in a snapshot.
    # To ensure reproducibility I avoid 'remove' here for now; it's also faster.
    iface = A.agents[a_agent]['iface']
    A.adjacency[a_agent] = [iface[s1]['bond'].split(A.bond_sep)[0] for s1 in iface if iface[s1]['bond'] != '.']
    iface = A.agents[b_agent]['iface']
    A.adjacency[b_agent] = [iface[s1]['bond'].split(A.bond_sep)[0] for s1 in iface if iface[s1]['bond'] != '.']

    if ka.system.canonicalize:
        # adjust the local views after bond loss
        A.agents[a_agent]['local_view'] = local_view(a_agent, A.agents)
        A.agents[b_agent]['local_view'] = local_view(b_agent, A.agents)

    # check for connectedness
    # dfs seems a tad bit faster than bfs when graphs are lightly connected;
    # need to check behavior for large densely connected aggregates
    component_agents = components(A, traverse='dfs')

    nc = len(component_agents)
    if nc == 2:
        # A fission has occurred. We create KappaMolecules reusing the fragments from the old;
        # id_shift = 0 by default.
        B = kamol.KappaMolecule(component_agents[0], count=0, system=ka.system, l_views=True)
        C = kamol.KappaMolecule(component_agents[1], count=0, system=ka.system, l_views=True)
        return 2, [B, C]

    else:
        # no fission has occurred
        # We need to update the KappaMolecule properties now that a bond has been lost.

        # remove from the bond dictionary
        del A.bonds[(port1, port2)]
        # standardize and update the bond *types*; labels don't matter
        b = sorted([port1, port2], key=lambda x: (alphanum_key(x[0]), alphanum_key(x[0])))
        ta, tb = kamol.bond2type(tuple(b))
        # Delete the bond from the bond list
        # This rigamarole is needed to make removal from a list O(1) rather than O(n).
        # We cannot use a dictionary, because in select_reaction() we need to randomly choose
        # from the available bonds, which is best done using a list... C'est la vie.
        remove = A.bond_list_idx[(ta, tb)][(port1, port2)]
        last = A.bond_list[(ta, tb)][-1]
        A.bond_list[(ta, tb)][remove] = last
        A.bond_list_idx[(ta, tb)][last] = remove
        A.bond_list[(ta, tb)].pop()
        A.bond_list_idx[(ta, tb)].pop((port1, port2))   # It's a dictionary, thus O(1)
        A.bond_type[(ta, tb)] -= 1
        # update degree
        A.agents[a_agent]['info']['degree'] -= 1
        A.agents[b_agent]['info']['degree'] -= 1
        # update free sites
        site1_type = ''.join([re.sub(r'.\d+.', '', a_agent), '.', a_site])
        site2_type = ''.join([re.sub(r'.\d+.', '', b_agent), '.', b_site])
        A.free_site_list[site1_type].append(port1)
        A.free_site_list[site2_type].append(port2)
        A.free_site_list_idx[site1_type][port1] = A.free_site[site1_type]
        A.free_site[site1_type] += 1
        A.free_site_list_idx[site2_type][port2] = A.free_site[site2_type]
        A.free_site[site2_type] += 1

        # kamol.sort_site_and_bond_lists(A)

        # update the agent self-binding counts
        iface_a = A.agents[a_agent]['iface']
        iface_b = A.agents[b_agent]['iface']
        for bt in A.signature.bond_types:
            st1, st2 = bt
            if st1 != st2:
                at1, s1 = st1.split('.')
                at2, s2 = st2.split('.')
                if st1 == site1_type and s2 in iface_a:
                    if iface_a[s2]['bond'] == '.':
                        A.agent_self_binding[bt] += 1
                if st2 == site1_type and s1 in iface_a:
                    if iface_a[s1]['bond'] == '.':
                        A.agent_self_binding[bt] += 1
                if st1 == site2_type and s2 in iface_b:
                    if iface_b[s2]['bond'] == '.':
                        A.agent_self_binding[bt] += 1
                if st2 == site2_type and s1 in iface_b:
                    if iface_b[s1]['bond'] == '.':
                        A.agent_self_binding[bt] += 1
        A.label_counter = int(kamol.get_identifier(next(reversed(A.agents)), delimiters=A.id_sep)[1])
        A.size = len(A.agents)

        if A.nav:
            # get the type lists for matching
            A.type_slice = []
            for at in A.composition:
                A.type_slice.extend([[name for name in A.agents if A.agents[name]['info']['type'] == at]])
            A.embedding_anchor = A.type_slice[0][0]
            # update navigation lists; we should not delete the entry if there is a second link
            A.make_navigation_list()

        if A.canon:
            lva = A.agents[a_agent]['local_view']
            if lva in A.local_views:
                A.local_views[lva].append(a_agent)
            else:
                A.local_views[lva] = [a_agent]
            lvb = A.agents[b_agent]['local_view']
            if lvb in A.local_views:
                A.local_views[lvb].append(b_agent)
            else:
                A.local_views[lvb] = [b_agent]
            A.canonical = A.canonicalize()

        # calculate reaction propensities
        A.internal_reactivity()
        return 1, [A, None]


def make_bond(A, B, A_port, B_port):
    """
    Binds KappaMolecules A and B by grafting B onto A's data structure. (Nothing is returned. The result is A.)
    B's agent labels must have been shifted by the number of agents in A. If the binding is intra-molecular,
    B is empty.
    """
    # A is the molecule onto which B is grafted;
    # B must have had its labels shifted by A.label_counter!

    # these are the binding ports
    (a_agent, a_site) = A_port
    (b_agent, b_site) = B_port

    # bimolecular case
    if B:
        # add agents of B to agents of A
        A.agents.update(B.agents)
        # update main lists and dictionaries
        for bt in B.bond_type:
            if B.bond_type[bt]:
                A.bond_type[bt] += B.bond_type[bt]
                n = len(A.bond_list[bt])
                for b in B.bond_list_idx[bt]:
                    B.bond_list_idx[bt][b] += n
                A.bond_list_idx[bt] = {**A.bond_list_idx[bt], **B.bond_list_idx[bt]}
                A.bond_list[bt].extend(B.bond_list[bt])
            # Further down, we will correct for a potential decrease
            # in self-binding due to the new bond.
            A.agent_self_binding[bt] += B.agent_self_binding[bt]
        for st in B.free_site:
            if B.free_site[st]:
                A.free_site[st] += B.free_site[st]
                n = len(A.free_site_list[st])
                for s in B.free_site_list_idx[st]:
                    B.free_site_list_idx[st][s] += n
                A.free_site_list_idx[st] = {**A.free_site_list_idx[st], **B.free_site_list_idx[st]}
                A.free_site_list[st].extend(B.free_site_list[st])
        A.bonds.update(B.bonds)
        # A.adjacency.update(B.adjacency)
        A.adjacency = {**A.adjacency, **B.adjacency}
        # get the composition
        A.get_composition()
        A.rarest_type = next(iter(A.composition))
        if A.nav:
            A.navigation.update(B.navigation)

    # make the bond
    iface_a = A.agents[a_agent]['iface']
    iface_b = A.agents[b_agent]['iface']
    iface_a[a_site]['bond'] = ''.join([b_agent, A.bond_sep, b_site])
    iface_b[b_site]['bond'] = ''.join([a_agent, A.bond_sep, a_site])
    # update adjacency
    A.adjacency[a_agent] = [iface_a[s1]['bond'].split(A.bond_sep)[0] for s1 in iface_a if iface_a[s1]['bond'] != '.']
    A.adjacency[b_agent] = [iface_b[s1]['bond'].split(A.bond_sep)[0] for s1 in iface_b if iface_b[s1]['bond'] != '.']
    if ka.system.canonicalize:
        # adjust the local views after bond loss
        A.agents[a_agent]['local_view'] = local_view(a_agent, A.agents)
        A.agents[b_agent]['local_view'] = local_view(b_agent, A.agents)
    # standardize the bond
    b = sorted([A_port, B_port], key=lambda x: (alphanum_key(x[0]), alphanum_key(x[1])))
    b = tuple(b)
    # update the bond list
    A.bonds[b] = 1
    (ta, tb) = kamol.bond2type(b)
    A.bond_list[(ta, tb)].append(b)
    A.bond_list_idx[(ta, tb)][b] = A.bond_type[(ta, tb)]
    A.bond_type[(ta, tb)] += 1
    # update degree
    A.agents[a_agent]['info']['degree'] += 1
    A.agents[b_agent]['info']['degree'] += 1
    # update free sites
    a_site_type = ''.join([re.sub(r'.\d+.', '', a_agent), '.', a_site])
    b_site_type = ''.join([re.sub(r'.\d+.', '', b_agent), '.', b_site])
    # Remove the newly freed sites from the site list
    # This rigamarole is needed to make removal from a list O(1) rather than O(n).
    # We cannot use a dictionary, because in select_reaction() we need to randomly choose
    # from the available sites, which is best done using a list... Mais oui.
    remove = A.free_site_list_idx[a_site_type][A_port]
    last = A.free_site_list[a_site_type][-1]
    A.free_site_list[a_site_type][remove] = last
    A.free_site_list_idx[a_site_type][last] = remove
    A.free_site_list[a_site_type].pop()
    A.free_site_list_idx[a_site_type].pop(A_port)  # It's a dictionary, thus O(1)
    A.free_site[a_site_type] -= 1
    remove = A.free_site_list_idx[b_site_type][B_port]
    last = A.free_site_list[b_site_type][-1]
    A.free_site_list[b_site_type][remove] = last
    A.free_site_list_idx[b_site_type][last] = remove
    A.free_site_list[b_site_type].pop()
    A.free_site_list_idx[b_site_type].pop(B_port)  # It's a dictionary, thus O(1)
    A.free_site[b_site_type] -= 1

    # kamol.sort_site_and_bond_lists(A)

    # agent self-binding correction
    for bt in A.signature.bond_types:
        st1, st2 = bt
        if st1 != st2:
            a1, s1 = st1.split('.')
            a2, s2 = st2.split('.')
            if st1 == a_site_type and s2 in iface_a:
                if iface_a[s2]['bond'] == '.':
                    A.agent_self_binding[bt] -= 1
            if st2 == a_site_type and s1 in iface_a:
                if iface_a[s1]['bond'] == '.':
                    A.agent_self_binding[bt] -= 1
            if st1 == b_site_type and s2 in iface_b:
                if iface_b[s2]['bond'] == '.':
                    A.agent_self_binding[bt] -= 1
            if st2 == b_site_type and s1 in iface_b:
                if iface_b[s1]['bond'] == '.':
                    A.agent_self_binding[bt] -= 1

    A.label_counter = int(kamol.get_identifier(next(reversed(A.agents)), delimiters=A.id_sep)[1])
    # size
    A.size = len(A.agents)

    if A.nav:
        # get the type lists for matching
        A.type_slice = []
        for at in A.composition:
            A.type_slice.extend([[name for name in A.agents if A.agents[name]['info']['type'] == at]])
        A.embedding_anchor = A.type_slice[0][0]
        # update the navigation list
        (a1, s1), (a2, s2) = b
        A.navigation[(a1, a2)] = s1
        A.navigation[(a2, a1)] = s2

    if A.canon:
        A.make_local_view_lists()
        A.canonical = A.canonicalize()

    # calculate reaction propensities
    A.internal_reactivity()


def bond_dissociation(reaction):
    """
    Organizes the execution of a bond dissolution.
    """
    choice, (molecule, _), (agent1, site1), (agent2, site2) = reaction
    # dissociation of the bond between (agent1, site1), (agent2, site2) in molecule
    # "molecule" refers to a species, not to an instance

    # keep in mind that we decreased the counts in the preceding negative update;
    # hence we correct counts by -1
    if molecule.count > 0:
        # we need to act on an instance, so let's make a copy
        A = kamol.copy_molecule(molecule, id_shift=0, system=ka.system)
    else:
        #  if molecule.count == 1, we modify in place
        A = molecule

    # execute the bond dissolution; 'A' might fragment.
    n_products, products = dissociate_bond(A, (agent1, site1), (agent2, site2))
    # if n_products == 1, 'A' was modified in place.
    # if n_products == 2, the agents of 'A' live on in the fragments...

    return n_products, products


def bimolecular_binding(reaction):
    """
    Organizes the execution of an inter-molecular binding event.
    """
    choice, (molecule1, molecule2), (agent1, site1), (agent2, site2) = reaction
    # In the cases below, we modify the recipient A in place; this is meant
    # to avoid (deep)copying a large molecule, since such a molecule is more likely
    # to be present as a single instance.
    # Notation: A is augmented by B

    # keep in mind that we decreased the counts in the preceding negative update;
    # hence we correct counts by -1
    if molecule1.count == 0 and molecule2.count == 0:
        # modify in place
        if molecule1.size > molecule2.size:
            # the recipient
            A = molecule1
            A_port = (agent1, site1)
            # the attachment
            B = kamol.copy_molecule(molecule2, id_shift=A.label_counter, system=ka.system)
            B_port = (agent2, site2)
            # molecule2 will be removed from the mixture; there might be equivalent molecules
            # in the mixture depending on consolidation status.
            # Note: the negative updates to reaction propensity were already done by calling negativeUpdate().
        else:
            A = molecule2
            A_port = (agent2, site2)
            B = kamol.copy_molecule(molecule1, id_shift=A.label_counter, system=ka.system)
            B_port = (agent1, site1)
            # molecule1 will be removed from the mixture; there might be equivalent molecules
            # in the mixture depending on compaction
            # Note: the negative updates to reaction propensity were already done by calling negativeUpdate().
    elif molecule1.count == 0:
        # modify in place
        A = molecule1
        A_port = (agent1, site1)
        # we could use copy.deepcopy(), but it seems slow; this one is fast
        B = kamol.copy_molecule(molecule2, id_shift=A.label_counter, system=ka.system)
        B_port = (agent2, site2)
    elif molecule2.count == 0:
        # modify in place
        A = molecule2
        A_port = (agent2, site2)
        B = kamol.copy_molecule(molecule1, id_shift=A.label_counter, system=ka.system)
        B_port = (agent1, site1)
    else:
        # construct a new KappaMolecule object to be added to the mixture
        A = kamol.copy_molecule(molecule1, system=ka.system)
        B = kamol.copy_molecule(molecule2, id_shift=A.label_counter, system=ka.system)
        A_port = (agent1, site1)
        B_port = (agent2, site2)
    # we first need to add the attachment B to the recipient's A data structure:
    # shift the identifiers in B by the largest identifier in A;
    # this affects B.agents{} dict
    # B.remap(id_shift=A.label_counter)
    # cognitive ergonomics...
    (b_agent, b_site) = B_port
    # remap the port name on the B side
    b_agent_type, b_agent_label = kamol.get_identifier(b_agent)
    b_agent_remapped = kamol.add_identifier(b_agent_type, str(int(b_agent_label) + A.label_counter))
    # execute the binding and update the data structure of A
    make_bond(A, B, A_port, (b_agent_remapped, b_site))

    return A


def unimolecular_binding(reaction):
    """
    Organizes the execution of an intra-molecular binding event.
    """
    choice, (molecule, _), (agent1, site1), (agent2, site2) = reaction
    # unimolecular binding between (agent1, site1), (agent2, site2) in molecule
    # "molecule" refers to a species, not to an instance

    # keep in mind that we decreased the counts in the preceding negative update;
    # hence we correct counts by -1
    if molecule.count == 0:
        # modify the only instance in place
        new_molecule = molecule  # new_molecule is just a new reference to the old object
    else:
        # create a copy and remove an old instance from the count
        # we'll modify the copy
        new_molecule = kamol.copy_molecule(molecule, system=ka.system)
        # can't have become zero, or we would be in the other branch
    # execute the binding
    make_bond(new_molecule, None, (agent1, site1), (agent2, site2))

    return new_molecule


def inflow(atom_type):
    """
    Injects an atom of atom_type into the mixture
    """
    atom = ka.system.kappa.parser(ka.system.signature.default_agent_state[atom_type])
    new_atom = kamol.KappaMolecule(atom, count=0, system=ka.system)

    return new_atom


def outflow(atom_type):
    """
    Ejects an atom of atom_type from the mixture
    """
    # Instance of the atom type must exist in the mixture or else we wouldn't be here.
    atom_canon = ka.system.mixture.atom_canonical[atom_type]
    m = ka.system.mixture.canonical[atom_canon]

    return m

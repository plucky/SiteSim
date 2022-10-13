# Walter Fontana 2022

from collections import deque
import sys

import kamol
import kasystem as ka


def components(kappaMol, traverse='bfs'):
    """
    Determines the graphical components of a KappaMolecule after bond dissociation.
    Returns an agent dictionary for each component (as requisite for creating a new KappaMolecule).

    This is standard bfs or dfs graph traversal. We use queues rather than recursion,
    because recursion seems very slow in Python.
    """
    visited = set()
    agents = []
    nc = 0

    # BFS with queue
    if traverse == 'bfs':
        queue = deque()  # [1st, 2nd, 3rd, 4th, ...]  add via append(); remove bia popleft()
        for node in kappaMol.agents:
            if node not in visited:
                # if nc == 1 at this point, we could exit, since dissociation
                # cannot generate more than 2 components. However, we may pay a higher cost
                # when recovering the components. To think about.
                nc += 1
                agents += [{}]
                queue.append(node)
                visited.add(node)
                while queue:
                    current = queue.popleft()  # a node in this component
                    agents[nc - 1][current] = kappaMol.agents[current]
                    for neighbor in kappaMol.adjacency[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
    # DFS with stack
    # Maybe use disjoint-sets data approach?
    elif traverse == 'dfs':
        stack = deque()  # [1st, 2nd, 3rd, 4th, ...]  add via append(); remove bia pop()
        for node in kappaMol.agents:
            if node not in visited:
                nc += 1
                agents += [{}]
                stack.append(node)
                while stack:
                    current = stack.pop()  # a node in this component
                    if current not in visited:
                        visited.add(current)
                        agents[nc - 1][current] = kappaMol.agents[current]
                        for neighbor in kappaMol.adjacency[current]:
                            if neighbor not in visited:
                                stack.append(neighbor)
    # expressions = []
    # for i in range(0, nc):
    #     expressions += [kappa_expression(agents[i], kappaMol.bonds, kappaMol.bondsep)]

    return agents  # a list of agent dictionaries, one for each component


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
    # update adjacency (will only remove one instance in case there are multiple)
    A.adjacency[a_agent].remove(b_agent)
    A.adjacency[b_agent].remove(a_agent)

    # check for connectedness
    component_agents = components(A, traverse='bfs')

    nc = len(component_agents)
    if nc == 2:
        # a fission has occurred
        # We create KappaMolecules from the debris.

        # We create new KappaMolecules reusing the fragments from the old; id_shift = 0 by default.
        B = kamol.KappaMolecule(component_agents[0], count=0, system=ka.system)
        C = kamol.KappaMolecule(component_agents[1], count=0, system=ka.system)
        return 2, [B, C]

    elif nc == 1:
        # no fission has occurred
        # We need to update the KappaMolecule properties now that a bond has been lost.

        # remove from the bond dictionary
        del A.bonds[(port1, port2)]
        # standardize and update the bond *types*; labels don't matter
        n1, l1 = kamol.get_identifier(a_agent)
        n2, l2 = kamol.get_identifier(b_agent)
        # possible exchange when n1 == n2, but doesn't matter in that case
        (t1, s1), (t2, s2) = sorted([(n1, a_site), (n2, b_site)])  # sorting not needed... (already sorted)
        ta, tb = (''.join([t1, '.', s1]), ''.join([t2, '.', s2]))
        A.bond_type[(ta, tb)] -= 1
        A.bond_type_list[(ta, tb)].remove((port1, port2))
        # update degree
        A.agents[a_agent]['info']['degree'] -= 1
        A.agents[b_agent]['info']['degree'] -= 1
        # update free sites
        site1_type = ''.join([n1, '.', a_site])  # this could be simply ta, if ports are guaranteed to be sorted...
        site2_type = ''.join([n2, '.', b_site])  # this could be simply tb, if ports are guaranteed to be sorted...
        A.free_site[site1_type] += 1
        A.free_site_list[site1_type] += [port1]
        A.free_site[site2_type] += 1
        A.free_site_list[site2_type] += [port2]
        # update the agent self-binding counts
        for bt in A.signature.bond_types:
            st1, st2 = bt
            if st1 != st2:
                at1, s1 = st1.split('.')
                at2, s2 = st2.split('.')
                if st1 == site1_type and s2 in A.agents[a_agent]['iface']:
                    if A.agents[a_agent]['iface'][s2]['bond'] == '.':
                        A.agent_self_binding[bt] += 1
                if st2 == site1_type and s1 in A.agents[a_agent]['iface']:
                    if A.agents[a_agent]['iface'][s1]['bond'] == '.':
                        A.agent_self_binding[bt] += 1
                if st1 == site2_type and s2 in A.agents[b_agent]['iface']:
                    if A.agents[b_agent]['iface'][s2]['bond'] == '.':
                        A.agent_self_binding[bt] += 1
                if st2 == site2_type and s1 in A.agents[b_agent]['iface']:
                    if A.agents[b_agent]['iface'][s1]['bond'] == '.':
                        A.agent_self_binding[bt] += 1
        A.label_counter = int(kamol.get_identifier(next(reversed(A.agents)), delimiters=A.id_sep)[1])
        # size
        A.size = len(A.agents)
        # get the composition
        A.get_composition()
        A.rarest_type = next(iter(A.composition))
        # make the adjacency lists. This should be updated more efficiently than a de novo construction...
        A.make_adjacency_lists()

        if A.nav:
            # get the type lists for matching
            A.type_slice = []
            for at in A.composition:
                A.type_slice.extend([[name for name in A.agents if A.agents[name]['info']['type'] == at]])
            A.embedding_anchor = A.type_slice[0][0]
            # update navigation lists
            # we should not delete the entry if there is a second link
            A.make_navigation_list()

        if A.canon:
            # requires adjacency list
            A.get_local_views()
            A.canonical = A.canonicalize()

        # calculate reaction propensities
        A.internal_reactivity()
        return 1, [A, None]
    else:
        sys.exit("Fission produced more than 2 components.")


def bind_molecules(A, B, A_port, B_port):
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
            A.bond_type[bt] += B.bond_type[bt]
            A.bond_type_list[bt].extend(B.bond_type_list[bt])
            # Further down, we will correct for a potential decrease
            # in self-binding due to the new bond.
            A.agent_self_binding[bt] += B.agent_self_binding[bt]
        for st in B.free_site:
            A.free_site[st] += B.free_site[st]
            A.free_site_list[st].extend(B.free_site_list[st])
        A.bonds.update(B.bonds)
        # A.adjacency.update(B.adjacency)
        if A.nav:
            A.navigation.update(B.navigation)
        # at this point we can delete the reference B
        # del B

    # make the bond
    A.agents[a_agent]['iface'][a_site]['bond'] = ''.join([b_agent, A.bond_sep, b_site])
    A.agents[b_agent]['iface'][b_site]['bond'] = ''.join([a_agent, A.bond_sep, a_site])

    # standardize the bond
    n1, l1 = kamol.get_identifier(a_agent)
    n2, l2 = kamol.get_identifier(b_agent)
    (t1, l1, s1), (t2, l2, s2) = sorted([(n1, int(l1), a_site), (n2, int(l2), b_site)])
    b = (kamol.add_identifier(t1, str(l1)), s1), (kamol.add_identifier(t2, str(l2)), s2)
    # update the bond dict
    A.bonds[b] = 1
    # update the adjacency list with the new bond
    # A.adjacency[a_agent].append(b_agent)
    # A.adjacency[b_agent].append(a_agent)
    # standardize and update the bond *types*; labels don't matter
    (t1, s1), (t2, s2) = sorted([(n1, a_site), (n2, b_site)])
    ta, tb = (''.join([t1, '.', s1]), ''.join([t2, '.', s2]))
    A.bond_type[(ta, tb)] += 1
    A.bond_type_list[(ta, tb)] += [b]
    # update degree
    A.agents[a_agent]['info']['degree'] += 1
    A.agents[b_agent]['info']['degree'] += 1
    # update free sites
    a_site_type = ''.join([n1, '.', a_site])
    b_site_type = ''.join([n2, '.', b_site])
    A.free_site[a_site_type] -= 1
    A.free_site_list[a_site_type].remove(A_port)
    A.free_site[b_site_type] -= 1
    A.free_site_list[b_site_type].remove(B_port)

    # agent self-binding correction
    for bt in A.signature.bond_types:
        st1, st2 = bt
        if st1 != st2:
            a1, s1 = st1.split('.')
            a2, s2 = st2.split('.')
            if st1 == a_site_type and s2 in A.agents[a_agent]['iface']:
                if A.agents[a_agent]['iface'][s2]['bond'] == '.':
                    A.agent_self_binding[bt] -= 1
            if st2 == a_site_type and s1 in A.agents[a_agent]['iface']:
                if A.agents[a_agent]['iface'][s1]['bond'] == '.':
                    A.agent_self_binding[bt] -= 1
            if st1 == b_site_type and s2 in A.agents[b_agent]['iface']:
                if A.agents[b_agent]['iface'][s2]['bond'] == '.':
                    A.agent_self_binding[bt] -= 1
            if st2 == b_site_type and s1 in A.agents[b_agent]['iface']:
                if A.agents[b_agent]['iface'][s1]['bond'] == '.':
                    A.agent_self_binding[bt] -= 1

    A.label_counter = int(kamol.get_identifier(next(reversed(A.agents)), delimiters=A.id_sep)[1])
    # size
    A.size = len(A.agents)
    # get the composition
    A.get_composition()
    A.rarest_type = next(iter(A.composition))
    # make the adjacency lists. This should be updated more efficiently than a de novo construction...
    A.make_adjacency_lists()

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
        # requires adjacency list
        A.get_local_views()
        A.canonical = A.canonicalize()

    # calculate reaction propensities
    A.internal_reactivity()


def bond_dissociation(reaction):
    """
    Organizes the execution of a bond dissolution. Care is taken (or so I believe) to free up memory...
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
    Organizes the execution of an inter-molecular binding event. Care is taken to free up memory...
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
    bind_molecules(A, B, A_port, (b_agent_remapped, b_site))

    return A


def unimolecular_binding(reaction):
    """
    Organizes the execution of an intra-molecular binding event. Care is taken to free up memory...
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
    bind_molecules(new_molecule, None, (agent1, site1), (agent2, site2))

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

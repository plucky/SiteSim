# Walter Fontana, 2022
"""
This module manages the 'mixture' of KappaMolecules.
"""
# population of Kappa complexes

import pprint
import kasnap as snap
import kasystem as ka
import kamol


def create_atoms(init_counts, signature=None, system=None):
    """
    Generates a simple mixture of monomeric agents in accordance with the signature.
    """
    kappa = kamol.Kappa()
    complexes = []
    local_views = {}
    for agent in signature.default_agent_state:
        count = init_counts[agent]
        if count:
            atom = kappa.parser(signature.default_agent_state[agent])
            kappa_atom = kamol.KappaMolecule(atom, count=count, system=system, sig=signature, views=local_views)
            complexes.append(kappa_atom)
    return 0, None, 0., complexes, local_views


class Mixture(snap.SnapShot):
    """
    Reads a snapshot from a .ka file and generates internal representations of molecules (complexes).
    The complexes are stored in the list self.complexes as KappaMolecule objects.
    Alternatively, import an already-made list of KappaMolecules (complexes) into the Mixture.
    """
    # To merge mixtures, we need to change the code a bit to ensure local views are inherited.

    def __init__(self, file=None, system=None, complexes=[]):
        super().__init__(file=file, complexes=complexes, system=system, signature=system.signature)

        # This is a dictionary version of self.complexes, but keys are canonicalized kappa expressions for fast
        # identification. It adds another "footprint of self.complexes" to the memory requirements. Not a problem.
        self.canonical = {}
        # This dictionary is for locating an 'atomic' species in 'complexes':
        # complexes[self.index[self.canonical[atom_canonical[atom_type]]]
        # We fill the dictionary over time.
        self.atom_canonical = {}
        # 'self.index' is a dictionary {molecule: index} indexing the list self.complexes. All operations
        # on 'self.complexes' must be mediated via 'self.index'.This allows us to delete an arbitrary molecule
        # entry of 'self,complexes' in O(1) _and_ we can use self.index to organize the heap for
        # efficient O(log n) reaction selection by the simulator.
        self.index = {}
        # note: self.local_views is declared in snap.Snapshot.

        if not file:
            value = create_atoms(ka.system.parameters.init_agents, signature=system.signature, system=system)
            self.event, self.origin_uuid, self.time, self.complexes, self.local_views = value
            self.origin_uuid = ka.system.uuid
            self.number_of_species = len(self.complexes)

        # make the heap index and the mapping from atom types to their canonical form
        for (i, mol) in enumerate(self.complexes):
            self.index[mol] = i  # indices start with 0
            self.canonical[mol.canonical] = mol

        # ------ only tracking
        self.total_bond_type = {}

        # ------ unimolecular reaction activities
        self.activity_unimolecular_binding = {}  # unimolecular activity: intra-molecular bond formation
        self.activity_bond_dissociation = {}     # unimolecular activity: bond dissociation

        self.unimolecular_reactivity_of_mixture()

        # ------ bimolecular reaction activities
        self.total_free_sites = {}
        self.activity_bimolecular_binding = {}  # bimolecular activity: inter-molecular bond formation

        self.bimolecular_reactivity_of_mixture()

        # ------ in/out flows (per atom type)
        self.activity_inflow = {}
        self.activity_outflow = {}

        self.flow_activity_of_mixture()

        # ------ total system activities
        self.total_activity = 0
        self.unimolecular_binding_activity = 0.
        self.bond_dissociation_activity = 0.
        self.bimolecular_binding_activity = 0.

        self.total_inflow = 0.
        self.total_outflow = 0.

        self.update_overall_activities()

    def unimolecular_reactivity_of_mixture(self):
        """
        Computes the unimolecular reactivity (dissociation and intra-molecular binding) of the mixture.
        """
        for bt in ka.system.signature.bond_types:
            self.activity_unimolecular_binding[bt] = 0.  # activity intra-molecular bond formation
            self.activity_bond_dissociation[bt] = 0.     # activity bond dissociation
            self.total_bond_type[bt] = 0                 # only tracking

        for m in self.complexes:
            for bt in ka.system.signature.bond_types:
                # multiplication with rate constant occurred when creating m
                self.activity_unimolecular_binding[bt] += (m.binding[bt] * m.count)
                # multiplication with rate constant occurred when creating m
                self.activity_bond_dissociation[bt] += (m.unbinding[bt] * m.count)
                # tracking
                self.total_bond_type[bt] += m.bond_type[bt]

    def bimolecular_reactivity_of_mixture(self):
        """
        Computes the bimolecular binding reactivity of the mixture.
        """
        for bt in ka.system.signature.bond_types:
            self.activity_bimolecular_binding[bt] = 0.  # activity inter-molecular bond formation
        # accumulate free_site counts across mixture
        for st in ka.system.signature.site_types:
            self.total_free_sites[st] = 0
        for m in self.complexes:
            for st in m.free_site:
                self.total_free_sites[st] += (m.free_site[st] * m.count)

        for m in self.complexes:   ### is this correct??
            for bt in ka.system.signature.bond_types:
                st1, st2 = bt
                factor = 1.
                if st1 == st2:  # symmetry correction
                    factor = 0.5
                # (Note to self: the below is array exo[] in the jupyter notebook)
                # inter-molecular bond formation *within* molecular species m
                a = m.free_site[st1] * m.free_site[st2] * (m.count - 1) * m.count
                # summing over all inter-molecular bond formation *between* molecular
                # species m and all others yields (see Overleaf notes for the calculation):
                a += m.free_site[st1] * m.count * (self.total_free_sites[st2] - m.free_site[st2] * m.count)
                self.activity_bimolecular_binding[bt] += (a * factor * ka.system.rc_bond_formation_inter)

    def flow_activity_of_mixture(self):
        """
        Tally agent in- and out-flow activity.
        """
        for a in ka.system.inflow_rate:
            # inflows are zero-molecular
            self.activity_inflow[a] = ka.system.inflow_rate[a]
        if ka.system.outflow_rate:
            for m in self.complexes:  # this is done only once, so we're good with looping
                if m.size == 1:
                    atom_type = m.agents[next(iter(m.agents))]['info']['type']
                    self.atom_canonical[atom_type] = m
                    if atom_type in ka.system.outflow_rate:
                        # outflows are unimolecular
                        self.activity_outflow[atom_type] = m.count * ka.system.outflow_rate[atom_type]

    def negativeUpdate(self, m):
        """
        This updates the aggregate binding and unbinding activities in the mixture stratified by bond or site type,
        as affected by the loss of a single molecule.
        """
        for bt in m.signature.bond_types:
            # unimolecular channels
            self.activity_unimolecular_binding[bt] -= m.binding[bt]
            self.activity_bond_dissociation[bt] -= m.unbinding[bt]
            # bimolecular channels
            st1, st2 = bt
            a = m.free_site[st1] * m.free_site[st2] * (m.count - 1)
            a += m.free_site[st1] * (self.total_free_sites[st2] - m.free_site[st2] * m.count)
            if st1 != st2:
                a += m.free_site[st2] * m.free_site[st1] * (m.count - 1)
                a += m.free_site[st2] * (self.total_free_sites[st1] - m.free_site[st1] * m.count)
            self.activity_bimolecular_binding[bt] -= a * ka.system.rc_bond_formation_inter
            # only tracking
            self.total_bond_type[bt] -= m.bond_type[bt]
        # update the total number of free sites per type
        for st in m.free_site:
            self.total_free_sites[st] -= m.free_site[st]

        # outflow is restricted to atoms (for now)
        if ka.system.outflow_rate and m.size == 1:
            agent_type = m.agents[next(iter(m.agents))]['info']['type']
            self.activity_outflow[agent_type] -= ka.system.outflow_rate[agent_type]
            self.total_outflow -= ka.system.outflow_rate[agent_type]

    def positiveUpdate(self, m):
        """
        This updates the aggregate binding and unbinding activities in the mixture stratified by bond or site type,
        as affected by the gain of a single molecule.
        """
        for st in m.free_site:
            self.total_free_sites[st] += m.free_site[st]
        for bt in m.signature.bond_types:
            # unimolecular channels
            self.activity_unimolecular_binding[bt] += m.binding[bt]
            self.activity_bond_dissociation[bt] += m.unbinding[bt]
            # bimolecular channels
            st1, st2 = bt
            a = m.free_site[st1] * m.free_site[st2] * (m.count - 1)
            a += m.free_site[st1] * (self.total_free_sites[st2] - m.free_site[st2] * m.count)
            if st1 != st2:
                a += m.free_site[st2] * m.free_site[st1] * (m.count - 1)
                a += m.free_site[st2] * (self.total_free_sites[st1] - m.free_site[st1] * m.count)
            self.activity_bimolecular_binding[bt] += a * ka.system.rc_bond_formation_inter
            # only tracking
            self.total_bond_type[bt] += m.bond_type[bt]

        # outflow is restricted to atoms (for now)
        if ka.system.outflow_rate and m.size == 1:
            agent_type = m.agents[next(iter(m.agents))]['info']['type']
            self.activity_outflow[agent_type] += ka.system.outflow_rate[agent_type]
            self.total_outflow += ka.system.outflow_rate[agent_type]

    def update_overall_activities(self):
        """
        Collects all activities of the three top channels:
            unimolecular binding
            bimolecular binding
            dissociation
        """
        self.total_activity = 0
        self.unimolecular_binding_activity = 0.
        self.bond_dissociation_activity = 0.
        self.bimolecular_binding_activity = 0.
        for bt in ka.system.signature.bond_types:
            self.unimolecular_binding_activity += self.activity_unimolecular_binding[bt]
            self.bond_dissociation_activity += self.activity_bond_dissociation[bt]
            self.bimolecular_binding_activity += self.activity_bimolecular_binding[bt]
        self.total_activity += self.unimolecular_binding_activity
        self.total_activity += self.bond_dissociation_activity
        self.total_activity += self.bimolecular_binding_activity

        self.total_inflow = 0.
        self.total_outflow = 0.
        for a in ka.system.inflow_rate:
            # inflows are zero-molecular
            self.total_inflow += self.activity_inflow[a]
        for a in ka.system.outflow_rate:
            # outflows are unimolecular
            self.total_outflow += self.activity_outflow[a]

        self.total_activity += self.total_inflow
        self.total_activity += self.total_outflow

    def remove_molecular_species(self, m):
        """
        Removes a molecular species from the mixture and syncs with heap.
        """
        # Note: only removes m from the list; m may still be referenced by other molecules.

        # The below amounts to self.complexes.remove(m). However, it executes _on average_ in O(1) time,
        # whereas self.complexes.remove(m) executes in O(n) with n the length of the list. We also use this
        # structure to build and maintain a heap on top.

        # save the index of the molecule to be deleted; we need it below to update the heap.
        old_index = self.index[m]
        # move the last molecular species to the location of the one being deleted.
        m_last = self.complexes[-1]
        self.complexes[self.index[m]] = m_last
        # update the index of the species we moved.
        self.index[m_last] = self.index[m]

        if ka.system.canonicalize:
            del self.canonical[m.canonical]
            if m.size == 1:
                atom_type = m.agents[next(iter(m.agents))]['info']['type']
                self.atom_canonical[atom_type] = None

        # finally, pop the last element from the list, as it has become redundant...
        self.complexes.pop()
        # ... and delete the index entry of the removed molecular species.
        self.index.pop(m)
        # update
        self.number_of_species -= 1
        # update the heap (using the index prior to deletion in complexes[])
        for t in ['bt+', 'bt-', 'st']:
            for k in ka.system.sim.heap[t]:
                ka.system.sim.heap[t][k].delete(old_index)

    def add_molecular_species(self, m):
        """
        Adds a molecular species to the mixture and syncs with heap.
        """
        self.index[m] = self.number_of_species
        self.complexes.append(m)

        if ka.system.canonicalize:
            self.canonical[m.canonical] = m
            if m.size == 1:
                atom_type = m.agents[next(iter(m.agents))]['info']['type']
                self.atom_canonical[atom_type] = m

        self.number_of_species += 1
        # update the heap
        for st in ka.system.sim.heap['st']:
            data_item = m.free_site[st] * m.count
            ka.system.sim.heap['st'][st].insert(data_item)
        for bt in ka.system.sim.heap['bt+']:
            data_item = m.binding[bt] * m.count
            ka.system.sim.heap['bt+'][bt].insert(data_item)
        for bt in ka.system.sim.heap['bt-']:
            data_item = m.unbinding[bt] * m.count
            ka.system.sim.heap['bt-'][bt].insert(data_item)

    def change_count(self, m, num, remove=True):
        """
        Updates the count and affected heaps.
        """
        m.count = m.count + num

        if m.count == 0 and remove:
            self.remove_molecular_species(m)
        else:
            # update the heap
            for st in ka.system.sim.heap['st']:
                data_item = m.free_site[st] * m.count
                ka.system.sim.heap['st'][st].modify(data_item, self.index[m])
            for bt in ka.system.sim.heap['bt+']:
                data_item = m.binding[bt] * m.count
                ka.system.sim.heap['bt+'][bt].modify(data_item, self.index[m])
            for bt in ka.system.sim.heap['bt-']:
                data_item = m.unbinding[bt] * m.count
                ka.system.sim.heap['bt-'][bt].modify(data_item, self.index[m])

    def __str__(self):
        self.report()

    def update_mixture(self, new):
        # add the reaction product 'm' to the mixture
        found = False
        if ka.system.canonicalize:
            if new.canonical in self.canonical:
                m = self.canonical[new.canonical]
                self.change_count(m, 1)
                found = True
        else:
            if ka.system.consolidate:
                # check if the new molecule belongs to an already existing species
                for m in self.complexes:
                    if m.size == m.size:
                        if m.sum_formula == new.sum_formula:
                            if ka.system.sgm.isomorphic(m, new):
                                # the molecular species exists; add the new molecule instance to it
                                self.change_count(m, 1)
                                found = True
                                break
        if not found:
            new.count = 1
            # append; multiple instances of the same species are not consolidated if consolidate=False
            self.add_molecular_species(new)

    def consolidate_mixture(self):
        """
        Determines equality between complexes and consolidates.
        """
        to_delete = []
        n_complexes = len(self.complexes)
        for i in range(0, n_complexes):
            m1 = self.complexes[i]
            for j in range(i+1, n_complexes):
                m2 = self.complexes[j]
                if m2.count == 0:
                    continue
                if m1.size == m2.size:
                    if m1.sum_formula == m2.sum_formula:
                        exists = False
                        if ka.system.canonicalize:
                            if m1.canonical == m2.canonical:
                                # the molecular species exists
                                exists = True
                        else:
                            if ka.system.sim.SGM.isomorphic(m1, m2):
                                # the molecular species exists
                                exists = True
                        if exists:
                            m1.count += m2.count
                            m2.count = 0
                            to_delete.append(m2)

        for m in to_delete:
            self.remove_molecular_species(m)

    def make_snapshot(self, file, label=False, pretty=False, sort=True):
        """
        Produces a snapshot of the mixture.
        'label' indicates whether kappa agents are labeled.
        'pretty' pretty-prints expressions to constant width
        """
        # consolidate temp_complexes? (use kasnap.consolidate())
        # don't mess up self.complexes, since it is indexed by self.index.
        if sort:
            temp_complexes = sorted(self.complexes, key=lambda c: c.size, reverse=True)
        else:
            temp_complexes = self.complexes
        with open(file, "w") as fp:
            s = f"// Snapshot [Event: {ka.system.sim.event}]\n"
            s += f'// "uuid" : "{ka.system.uuid}"\n'
            s += f'%def: "T0" "{ka.system.sim.time}"\n'
            s += '\n'
            fp.write(s)
            for m in temp_complexes:
                s = f"%init: {m.count} /*{m.size} agents*/ "
                if pretty:
                    temp = pprint.pformat(m.kappa_expression(), indent=1, width=60)
                    s += temp[1:-1].replace("'", "") + '\n'
                else:
                    s += m.kappa_expression(label=label) + '\n'
                fp.write(s)

    def report(self, pp_width=40):
        """
        Summarizes the mixture.
        """
        info = f"\n{'MIXTURE '.ljust(70, '-')}\n\n"
        info += f'{"initial mixture file":>20}: {self.file}\n'
        info += f'{"molecular species":>20}: {self.number_of_species}\n'
        info += f'{"agents, molecules":>20}: {self.count_agents_and_molecules()}\n'

        if ka.system.db_level == 1:
            info += f'{"free sites":>20}:\n'
            for st in self.total_free_sites:
                info += f'{st:>20}: {self.total_free_sites[st]}\n'

        info += f'{"size distribution":>20}:\n'
        size_dist = self.get_size_distribution()
        d1 = f'{size_dist[:4]}'
        d2 = f'{size_dist[-4:]}'
        info += f'{d1[1:-1]} ... {d2[1:-1]}'
        info += '\n'

        info += f'{"system activities ":>{pp_width}}\n'
        info += f'{"total system activity":>{pp_width}}: {self.total_activity:1.5E}\n'
        info += f'{"unimolecular binding activity":>{pp_width}}: {self.unimolecular_binding_activity:1.5E}\n'
        info += f'{"bond dissociation activity":>{pp_width}}: {self.bond_dissociation_activity:1.5E}\n'
        info += f'{"bimolecular binding activity":>{pp_width}}: {self.bimolecular_binding_activity:1.5E}\n'
        info += '\n'
        info += f'{"system activities by bond type ":>{pp_width}}\n'
        info += f'{"unimolecular binding activity ":>{pp_width}}\n'
        for bt in ka.system.signature.bond_types:
            s = f'{bt}'
            info += f'{s:>{pp_width}}: {self.activity_unimolecular_binding[bt]:1.5E}\n'
        info += f'{"bond dissociation activity ":>{pp_width}}\n'
        for bt in ka.system.signature.bond_types:
            s = f'{bt}'
            info += f'{s:>{pp_width}}: {self.activity_bond_dissociation[bt]:1.5E}\n'
        info += f'{"bimolecular binding activity ":>{pp_width}}\n'
        for bt in ka.system.signature.bond_types:
            s = f'{bt}'
            info += f'{s:>{pp_width}}: {self.activity_bimolecular_binding[bt]:1.5E}\n'
        info += f'{"inflow activity ":>{pp_width}}\n'
        for a in ka.system.inflow_rate:
            s = f'atom type {a}'
            info += f'{s:>{pp_width}}: {self.activity_inflow[a]:1.5E}\n'
        info += f'{"outflow activity ":>{pp_width}}\n'
        for a in ka.system.outflow_rate:
            s = f'atom type {a}'
            info += f'{s:>{pp_width}}: {self.activity_outflow[a]:1.5E}\n'
        info += '\n'

        reactivity = False
        show_bonds = False
        if ka.system.db_level > 0:
            if ka.system.db_level == 1:
                reactivity = True
                show_bonds = False
            elif ka.system.db_level == 2:
                reactivity = True
                show_bonds = True
            temp = sorted(self.complexes, key=lambda c: c.size, reverse=True)
            info += f"\n\n{'MIXTURE CONTENTS'.ljust(70, '-')}\n"
            info += '\n'
            n = 0
            for m in temp:
                s = '#' + str(n)
                info += (s + m.summary(reactivity=reactivity, show_bonds=show_bonds)[len(s)+1:])
                n += 1
        return info
        
# -------------------------------------------------------------------------------------------


if __name__ == '__main__':
    import kainit

    signature = 'A(p[a1.P$m a2.P$m a3.P$m], l[r.A$w], r[l.A]), P(a1[p.A], a2[p.A], a3[p.A], d[d.P$m])'
    parameters = 'TestData/parameters.txt'
    kainit.minit(signature=signature, parameter_file=parameters)

    mix = Mixture('TestData/snap__0181.ka')
    print(f'{mix.number_of_species} complexes')
    # size distribution
    dist = mix.get_size_distribution()
    for item in dist:
        print(f'size: {item[0]}   count: {item[1]}')
    print(' ')
    print(mix.report())
    print(' ')

    # for c in mix.complexes:
    #     print(c.kappa_expression())
    #     print(' ')
    #     print("free sites")
    #     for s in c.free_site:
    #         print(f'{s} -> {c.free_site[s]}')
    #     print(' ')
    #     print("bond types")
    #     n = 0
    #     for b in c.bond_type:
    #         print(f'{b} -> {c.bond_type[b]}')
    #         n += c.bond_type[b]
    #     if n != len(c.bonds):
    #         print('ERROR !!')
    #         sys.exit(-1)
    #     else:
    #         print(f'total bonds: {n} of {len(c.bonds)}')
    #     print(' ')

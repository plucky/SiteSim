# Walter Fontana, 2022
""""
This module defines the simulation engine.
"""

import kasystem as ka
import kaheap
import kareact as react

import numpy as np


class Random:
    """
    Sets up random number generators.
    (The 'pythonista' flag is for working with Pythonista on the iPad.)
    """
    def __init__(self, seed, pythonista=False):
        self.pythonista = pythonista
        # self.rng = np.random.default_rng(seed=ka.system.parameters.rng_seed)  # same as below
        # use np.random.MT19937(seed) for RandomState compatibility
        if self.pythonista:
            self.rng = np.random.RandomState(seed=seed)
        else:
            self.rng = np.random.Generator(np.random.PCG64(seed=seed))

    def integers(self, low=0, high=1):
        if self.pythonista:
            return self.rng.randint(low=low, high=high)
        else:
            return self.rng.integers(low=low, high=high)

    def uniform(self, low=0., high=1.):
        return self.rng.uniform(low=low, high=high)

    def exponential(self, scale=1.):
        return self.rng.exponential(scale=scale)


class CTMC:
    """
    Defines the continuous-time Monte Carlo simulator.
    """
    def __init__(self):

        self.time = ka.system.mixture.time  # we inherit the initial time from the mixture
        self.event = ka.system.mixture.event  # we inherit the initial event number from the mixture

        self.current_reaction = None

        # initialize random number generator
        self.rng = np.random.default_rng(seed=ka.system.parameters.rng_seed)
        # self.rng = Random(seed=ka.system.parameters.rng_seed, pythonista=True)

        # collection of heaps to speed up reaction selection
        self.heap = {'bt+': {}, 'bt-': {}, 'st': {}}
        self.initialize_heaps()

    def initialize_heaps(self):
        """
        Initialize the collection of heaps.
        """
        sites = set()
        for bt in ka.system.signature.bond_types:
            s1, s2 = bt
            sites.update([s1, s2])

            # heaps for handling reaction selection based on binding (bt+) and unbinding (bt-)
            # stratified by binding type
            data = [m.binding[bt] * m.count for m in ka.system.mixture.complexes]
            self.heap['bt+'][bt] = kaheap.Heap(data, ident=f'bt+ | {bt}')
            data = [m.unbinding[bt] * m.count for m in ka.system.mixture.complexes]
            self.heap['bt-'][bt] = kaheap.Heap(data, ident=f'bt- | {bt}')

        # heaps for handling reaction selection based on bimolecular binding stratified by site type
        for s in sites:
            data = [m.free_site[s] * m.count for m in ka.system.mixture.complexes]
            self.heap['st'][s] = kaheap.Heap(data, ident=f'st | {s}')

    def advance_time(self):
        """
        Simulates time.
        """
        self.time += self.rng.exponential(scale=1. / ka.system.mixture.total_activity)

    def execute_reaction(self):
        """
        Executes a reaction and updates.
        """
        choice, (molecule1, molecule2), (agent1, site1), (agent2, site2) = self.current_reaction

        if choice == 'inflow':
            # Inflow of an atom of type declared in 'molecule1'. (In this case, molecule1 is a string.)

            new = react.inflow(molecule1)
            ka.system.mixture.update_mixture(new)
            ka.system.mixture.positiveUpdate(new)
            ka.system.mixture.update_overall_activities()

        elif choice == 'outflow':
            # Outflow of an atom of type declared in 'molecule1'. (In this case, molecule1 is a string.)

            out = react.outflow(molecule1)
            ka.system.mixture.negativeUpdate(out)
            ka.system.mixture.change_count(out, -1)
            ka.system.mixture.update_overall_activities()

        elif choice == 'ub':
            # Unimolecular binding between site1 of agent1 and site2 of agent2 in molecule1.
            # molecule2 is None.

            # execute the reaction by creating the new molecule(s)
            new = react.unimolecular_binding(self.current_reaction)
            # negative update of propensities (before updating reactant counts!)
            ka.system.mixture.negativeUpdate(molecule1)
            # adjust actual counts (and remove species if count drops to zero)
            ka.system.mixture.change_count(molecule1, -1)
            # add the product to the mixture (updates product counts)
            ka.system.mixture.update_mixture(new)
            # positive update of propensities (before updating product counts!)
            ka.system.mixture.positiveUpdate(new)
            # update overall propensities
            ka.system.mixture.update_overall_activities()

        elif choice == 'bd':
            # Bond dissociation between (agent1, site1), (agent2, site2) in molecule1

            # update logic as for case 'ub'
            n_products, product = react.bond_dissociation(self.current_reaction)
            ka.system.mixture.negativeUpdate(molecule1)
            ka.system.mixture.change_count(molecule1, -1)
            for i in range(0, n_products):
                ka.system.mixture.update_mixture(product[i])
                ka.system.mixture.positiveUpdate(product[i])
            ka.system.mixture.update_overall_activities()

        elif choice == 'bb':
            # Bimolecular binding between (agent1, site1) in molecule1 and (agent2, site2) in molecule2

            # update logic as for case 'ub'
            new = react.bimolecular_binding(self.current_reaction)
            ka.system.mixture.negativeUpdate(molecule1)
            ka.system.mixture.change_count(molecule1, -1)
            ka.system.mixture.negativeUpdate(molecule2)
            ka.system.mixture.change_count(molecule2, -1)
            ka.system.mixture.update_mixture(new)
            ka.system.mixture.positiveUpdate(new)
            ka.system.mixture.update_overall_activities()

        else:
            return
            # sys.exit("execute_reaction(): Unknown reaction")

    # def select_reaction_basic(self):
    #     """
    #     Obsolete. Supplanted by the version with heaps. Kept for testing.
    #     """
    #     # The first few choices relate to in/outflow of atoms.
    #     # Since there are only a few types, looping is OK.
    #     rv = self.rng.uniform(low=0.0, high=ka.system.mixture.total_activity)
    #     if rv < ka.system.mixture.total_inflow:
    #         select = 'inflow'
    #         for a in ka.system.mixture.activity_inflow:
    #             if rv < ka.system.mixture.activity_inflow[a]:
    #                 self.current_reaction = select, (a, None), (None, None), (None, None)
    #                 return
    #             else:
    #                 rv -= ka.system.mixture.activity_inflow[a]
    #
    #     rv -= ka.system.mixture.total_inflow
    #     if rv < ka.system.mixture.total_outflow:
    #         select = 'outflow'
    #         for a in ka.system.mixture.activity_outflow:
    #             if rv < ka.system.mixture.activity_outflow[a]:
    #                 self.current_reaction =  select, (a, None), (None, None), (None, None)
    #                 return
    #             else:
    #                 rv -= ka.system.mixture.activity_outflow[a]
    #
    #     rv = self.rng.uniform(low=0.0, high=ka.system.mixture.total_activity)
    #     if rv < ka.system.mixture.unimolecular_binding_activity:
    #         # channel is a unimolecular binding event
    #         select = 'ub'
    #         for bt in ka.system.signature.bond_types:
    #             if rv < ka.system.mixture.activity_unimolecular_binding[bt]:
    #                 # the internal bond to be formed is of type bt
    #                 # refine search to the molecular level
    #                 for m in ka.system.mixture.complexes:
    #                     segment = m.binding[bt] * m.count
    #                     if rv < segment:
    #                         # the event is within molecular species m;
    #                         # now we uniformly choose the instance of bt in m
    #                         # bt is (agent_type1.site1), (agent_type2.site2)
    #                         # NOTE: I don't think we need to randomize which site we choose first.
    #                         s1, s2 = bt
    #                         r1 = self.rng.integers(low=0, high=m.free_site[s1])
    #                         # this is our choice of site1; it belongs to agent name1 in molecule m
    #                         name1, site1 = m.free_site_list[s1][r1]
    #                         if s1 == s2:
    #                             temp_list = [p for p in m.free_site_list[s2] if p != (name1, site1)]
    #                             r2 = self.rng.integers(low=0, high=m.free_site[s2] - 1)
    #                             name2, site2 = temp_list[r2]
    #                         else:
    #                             # exclude possibility of self-binding
    #                             site2 = s2.split('.')[1]
    #                             if (name1, site2) in m.free_site_list[s2]:
    #                                 temp_list = [p for p in m.free_site_list[s2] if p != (name1, site2)]
    #                                 r2 = self.rng.integers(low=0, high=m.free_site[s2] - 1)
    #                                 name2, site2 = temp_list[r2]
    #                                 self.current_reaction =  select, (m, None), (name1, site1), (name2, site2)
    #                                 return
    #                             else:
    #                                 r2 = self.rng.integers(low=0, high=m.free_site[s2])
    #                                 name2, site2 = m.free_site_list[s2][r2]
    #                         self.current_reaction =  select, (m, None), (name1, site1), (name2, site2)
    #                         return
    #                     else:
    #                         rv -= segment
    #             else:
    #                 rv -= ka.system.mixture.activity_unimolecular_binding[bt]
    #
    #     # Note to self: navigation across intervals may need to be changed to accommodate size-dependent alpha
    #     rv -= ka.system.mixture.unimolecular_binding_activity
    #     if rv < ka.system.mixture.bond_dissociation_activity:
    #         # channel is a bond dissociation
    #         select = 'bd'
    #         for bt in ka.system.signature.bond_types:
    #             if rv < ka.system.mixture.activity_bond_dissociation[bt]:
    #                 # refine search to molecular level
    #                 for m in ka.system.mixture.complexes:
    #                     segment = m.unbinding[bt] * m.count
    #                     if rv < segment:
    #                         # it's molecule of species m; now we need to uniformly choose the instance of bt in m
    #                         # bt is (agent_type_1.site_1), (agent_type_2.site_2)
    #                         # choose a bond of this type
    #                         r = self.rng.integers(low=0, high=m.bond_type[bt])
    #                         x, y = m.bond_type_list[bt][r]
    #                         self.current_reaction =  select, (m, None), x, y
    #                         return
    #                     else:
    #                         rv -= segment
    #             else:
    #                 rv -= ka.system.mixture.activity_bond_dissociation[bt]
    #
    #     rv -= ka.system.mixture.bond_dissociation_activity
    #     if rv < ka.system.mixture.bimolecular_binding_activity:
    #         # channel is a bimolecular binding
    #         select = 'bb'
    #         for bt in ka.system.signature.bond_types:
    #             if rv < ka.system.mixture.activity_bimolecular_binding[bt]:
    #                 s1, s2 = bt
    #                 # choose at random (uniformly) an s1, i.e. an agent and free site of required type
    #                 r1 = self.rng.integers(low=0, high=ka.system.mixture.total_free_sites[s1])
    #                 for m1 in ka.system.mixture.complexes:
    #                     segment = m1.free_site[s1] * m1.count
    #                     if r1 < segment:
    #                         # it's m1
    #                         break
    #                     else:
    #                         r1 -= segment
    #                 r2 = self.rng.integers(low=0, high=ka.system.mixture.total_free_sites[s2] - m1.free_site[s2])
    #                 m1.count -= 1  # only temporary!
    #                 for m2 in ka.system.mixture.complexes:
    #                     segment = m2.free_site[s2] * m2.count
    #                     if r2 < segment:
    #                         # it's m2
    #                         break
    #                     else:
    #                         r2 -= segment
    #                 m1.count += 1  # undo
    #                 r1 = self.rng.integers(low=0, high=m1.free_site[s1])
    #                 r2 = self.rng.integers(low=0, high=m2.free_site[s2])
    #                 name1, site1 = m1.free_site_list[s1][r1]
    #                 name2, site2 = m2.free_site_list[s2][r2]
    #                 self.current_reaction =  select, (m1, m2), (name1, site1), (name2, site2)
    #                 return
    #             else:
    #                 rv -= ka.system.mixture.activity_bimolecular_binding[bt]

    def select_reaction(self):
        """
        Fast reaction selection using heaps.
        """
        # The first few choices relate to in/outflow of atoms.
        # Since there are only a few types, looping is OK.
        rv = self.rng.uniform(low=0.0, high=ka.system.mixture.total_activity)
        if rv < ka.system.mixture.total_inflow:
            select = 'inflow'
            for a in ka.system.mixture.activity_inflow:
                if rv < ka.system.mixture.activity_inflow[a]:
                    self.current_reaction =  select, (a, None), (None, None), (None, None)
                    return
                else:
                    rv -= ka.system.mixture.activity_inflow[a]

        rv -= ka.system.mixture.total_inflow
        if rv < ka.system.mixture.total_outflow:
            select = 'outflow'
            for a in ka.system.mixture.activity_outflow:
                if rv < ka.system.mixture.activity_outflow[a]:
                    self.current_reaction =  select, (a, None), (None, None), (None, None)
                    return
                else:
                    rv -= ka.system.mixture.activity_outflow[a]

        rv -= ka.system.mixture.total_outflow
        if rv < ka.system.mixture.unimolecular_binding_activity:
            # channel is a unimolecular binding event
            select = 'ub'
            for bt in ka.system.signature.bond_types:
                if rv < ka.system.mixture.activity_unimolecular_binding[bt]:
                    # the internal bond to be formed is of type bt
                    # refine search to the molecular level
                    m = ka.system.mixture.complexes[self.heap['bt+'][bt].draw_node(rv)]
                    # the event is within molecular species m;
                    # now we uniformly choose the instance of bt in m
                    # bt is (agent_type1.site1), (agent_type2.site2)
                    # NOTE: I don't think we need to randomize which site we choose first.
                    s1, s2 = bt
                    r1 = self.rng.integers(low=0, high=m.free_site[s1])
                    # this is our choice of site1; it belongs to agent name1 in molecule m
                    name1, site1 = m.free_site_list[s1][r1]
                    if s1 == s2:
                        temp_list = [p for p in m.free_site_list[s2] if p != (name1, site1)]
                        r2 = self.rng.integers(low=0, high=m.free_site[s2] - 1)
                        name2, site2 = temp_list[r2]
                    else:
                        # exclude possibility of self-binding
                        site2 = s2.split('.')[1]
                        if (name1, site2) in m.free_site_list[s2]:
                            temp_list = [p for p in m.free_site_list[s2] if p != (name1, site2)]
                            r2 = self.rng.integers(low=0, high=m.free_site[s2] - 1)
                            name2, site2 = temp_list[r2]
                            self.current_reaction =  select, (m, None), (name1, site1), (name2, site2)
                            return
                        else:
                            r2 = self.rng.integers(low=0, high=m.free_site[s2])
                            name2, site2 = m.free_site_list[s2][r2]
                    self.current_reaction =  select, (m, None), (name1, site1), (name2, site2)
                    return
                else:
                    rv -= ka.system.mixture.activity_unimolecular_binding[bt]

        # Note to self: navigation across intervals may need to be changed to accommodate size-dependent alpha
        rv -= ka.system.mixture.unimolecular_binding_activity
        if rv < ka.system.mixture.bond_dissociation_activity:
            # channel is a bond dissociation
            select = 'bd'
            for bt in ka.system.signature.bond_types:
                if rv < ka.system.mixture.activity_bond_dissociation[bt]:
                    # refine search to molecular level
                    m = ka.system.mixture.complexes[self.heap['bt-'][bt].draw_node(rv)]
                    # it's molecule of species m; now we need to uniformly choose the instance of bt in m
                    # bt is (agent_type_1.site_1), (agent_type_2.site_2)
                    # choose a bond of this type
                    r = self.rng.integers(low=0, high=m.bond_type[bt])
                    x, y = m.bond_type_list[bt][r]
                    self.current_reaction =  select, (m, None), x, y
                    return
                else:
                    rv -= ka.system.mixture.activity_bond_dissociation[bt]

        rv -= ka.system.mixture.bond_dissociation_activity
        if rv < ka.system.mixture.bimolecular_binding_activity:
            # channel is a bimolecular binding
            select = 'bb'
            for bt in ka.system.signature.bond_types:
                if rv < ka.system.mixture.activity_bimolecular_binding[bt]:
                    s1, s2 = bt
                    # choose at random (uniformly) an s1, i.e. an agent and free site of required type
                    r1 = self.rng.integers(low=0, high=ka.system.mixture.total_free_sites[s1])
                    m1 = ka.system.mixture.complexes[self.heap['st'][s1].draw_node(r1)]
                    r2 = self.rng.integers(low=0, high=ka.system.mixture.total_free_sites[s2] - m1.free_site[s2])
                    ka.system.mixture.change_count(m1, -1, remove=False)  # only temporary!
                    m2 = ka.system.mixture.complexes[self.heap['st'][s2].draw_node(r2)]
                    ka.system.mixture.change_count(m1, 1)  # undo!
                    r1 = self.rng.integers(low=0, high=m1.free_site[s1])
                    r2 = self.rng.integers(low=0, high=m2.free_site[s2])
                    name1, site1 = m1.free_site_list[s1][r1]
                    name2, site2 = m2.free_site_list[s2][r2]
                    self.current_reaction =  select, (m1, m2), (name1, site1), (name2, site2)
                    return
                else:
                    rv -= ka.system.mixture.activity_bimolecular_binding[bt]

    def report(self, pp_width=40):
        info = f"\n{'SIMULATOR STATE '.ljust(70, '-')}\n\n"
        # info += f'{"system initialized on:":>{pp_width}} {ka.system.date}\n'
        # info += f'{"system uuid:":>{pp_width}} {ka.system.uuid}\n'
        # info += '\n'
        info += f'{"simulator status at time t=":>{pp_width}} {self.time}\n'
        info += '\n'
        info += f'{"total system activity":>{pp_width}}: {ka.system.mixture.total_activity:1.5E}\n'
        info += '\n'
        n_heaps = len(self.heap["bt+"]) + len(self.heap["bt-"]) + len(self.heap["st"])
        if n_heaps > 0:
            info += f'{"heaps":>{pp_width}}: {n_heaps} x ['
            a = next(iter(self.heap))
            b = next(iter(self.heap[a]))
            n_nodes = self.heap[a][b].n_internal_nodes + self.heap[a][b].n_leaves
            info += f"height: {self.heap[a][b].height} nodes: {n_nodes} "
            info += f"occ: {self.heap[a][b].n_entries / self.heap[a][b].n_leaves:.2f}]\n"
            for t in ['bt+', 'bt-', 'st']:
                for k in self.heap[t]:
                    info += f"HEAP {self.heap[t][k].id} -> root value: {self.heap[t][k].tree[0]:1.2E}\n"
        return info

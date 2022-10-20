# Walter Fontana, 2022
""""
This module defines the simulation engine.
"""

import kasystem as ka
import kaheap
import kareact as react

import numpy as np
import sys
import json


class CTMC:
    """
    Defines the continuous-time Monte Carlo simulator.
    """
    def __init__(self):

        self.sys = ka.system
        self.mix = ka.system.mixture
        self.sig = ka.system.signature

        self.time = ka.system.mixture.time  # we inherit the initial time from the mixture
        self.event = ka.system.mixture.event  # we inherit the initial event number from the mixture

        self.current_reaction = None

        # initialize random number generator
        self.rng = np.random.Generator(np.random.PCG64(seed=ka.system.parameters.rng_seed))
        # self.rng = np.random.default_rng(seed=ka.system.parameters.rng_seed)
        # [use np.random.MT19937(seed) for compatibility with RandomState]
        # restore state, if desired (mainly for continuation of a simulation)
        if ka.system.mixture.rg_state:
            self.rng.bit_generator.state = ka.system.mixture.rg_state

        # collection of heaps to speed up reaction selection
        self.heap = {'bt+': {}, 'bt-': {}, 'st': {}}
        self.initialize_heaps()

    def initialize_heaps(self):
        """
        Initialize the collection of heaps.
        """
        sites = set()
        for bt in self.sig.bond_types:
            s1, s2 = bt
            sites.update([s1, s2])

            # heaps for handling reaction selection based on binding (bt+) and unbinding (bt-)
            # stratified by binding type
            data = [m.binding[bt] * m.count for m in self.mix.complexes]
            self.heap['bt+'][bt] = kaheap.Heap(data, ident=f'bt+ | {bt}')
            data = [m.unbinding[bt] * m.count for m in self.mix.complexes]
            self.heap['bt-'][bt] = kaheap.Heap(data, ident=f'bt- | {bt}')

        # heaps for handling reaction selection based on bimolecular binding stratified by site type
        for s in sites:
            data = [m.free_site[s] * m.count for m in self.mix.complexes]
            self.heap['st'][s] = kaheap.Heap(data, ident=f'st | {s}')

    def advance_time(self):
        """
        Simulates time.
        """
        self.time += self.rng.exponential(scale=1. / self.mix.total_activity)

    def execute_reaction(self):
        """
        Executes a reaction and updates.
        """
        choice, (molecule1, molecule2), (agent1, site1), (agent2, site2) = self.current_reaction

        if choice == 'ub':
            # Unimolecular binding between site1 of agent1 and site2 of agent2 in molecule1.
            # molecule2 is None.

            # negative update of propensities (before updating reactant counts!)
            self.mix.negativeUpdate(molecule1)
            # adjust actual counts (and remove species if count drops to zero)
            self.mix.change_count(molecule1, -1)
            # execute the reaction by creating the new molecule(s)
            new = react.unimolecular_binding(self.current_reaction)
            # add the product to the mixture (updates product counts)
            self.mix.update_mixture(new)
            # positive update of propensities (before updating product counts!)
            self.mix.positiveUpdate(new)

        elif choice == 'bd':
            # Bond dissociation between (agent1, site1), (agent2, site2) in molecule1

            self.mix.negativeUpdate(molecule1)
            self.mix.change_count(molecule1, -1)
            n_products, product = react.bond_dissociation(self.current_reaction)
            for i in range(0, n_products):
                self.mix.update_mixture(product[i])
                self.mix.positiveUpdate(product[i])

        elif choice == 'bb':
            # Bimolecular binding between (agent1, site1) in molecule1 and (agent2, site2) in molecule2

            # update logic as for case 'ub'
            self.mix.negativeUpdate(molecule1)
            self.mix.change_count(molecule1, -1)
            self.mix.negativeUpdate(molecule2)
            self.mix.change_count(molecule2, -1)
            new = react.bimolecular_binding(self.current_reaction)
            self.mix.update_mixture(new)
            self.mix.positiveUpdate(new)

        elif choice == 'inflow':
            # Inflow of an atom of type declared in 'molecule1'. (In this case, molecule1 is a string.)

            new = react.inflow(molecule1)
            self.mix.update_mixture(new)
            self.mix.positiveUpdate(new)

        elif choice == 'outflow':
            # Outflow of an atom of type declared in 'molecule1'. (In this case, molecule1 is a string.)

            out = react.outflow(molecule1)
            self.mix.negativeUpdate(out)
            self.mix.change_count(out, -1)

        else:
            sys.exit("execute_reaction: Unknown reaction")

        # update overall propensities
        self.mix.update_overall_activities()

    def select_reaction(self):
        """
        Fast reaction selection using heaps.
        """

        rv = self.rng.uniform(low=0.0, high=self.mix.total_activity)

        if rv < self.mix.unimolecular_binding_activity:
            # channel is a unimolecular binding event
            select = 'ub'
            for bt in self.sig.bond_types:
                if rv < self.mix.activity_unimolecular_binding[bt]:
                    # the internal bond to be formed is of type bt
                    # refine search to the molecular level
                    m = self.mix.complexes[self.heap['bt+'][bt].draw_node(rv)]
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
                    self.current_reaction = select, (m, None), (name1, site1), (name2, site2)
                    return
                else:
                    rv -= self.mix.activity_unimolecular_binding[bt]

        # Note to self: navigation across intervals may need to be changed to accommodate size-dependent alpha
        rv -= self.mix.unimolecular_binding_activity
        if rv < self.mix.bond_dissociation_activity:
            # channel is a bond dissociation
            select = 'bd'
            for bt in self.sig.bond_types:
                if rv < self.mix.activity_bond_dissociation[bt]:
                    # refine search to molecular level
                    m = self.mix.complexes[self.heap['bt-'][bt].draw_node(rv)]
                    # it's molecule of species m; now we need to uniformly choose the instance of bt in m
                    # bt is (agent_type_1.site_1), (agent_type_2.site_2)
                    # choose a bond of this type
                    r = self.rng.integers(low=0, high=m.bond_type[bt])
                    x, y = m.bond_type_list[bt][r]
                    self.current_reaction =  select, (m, None), x, y
                    return
                else:
                    rv -= self.mix.activity_bond_dissociation[bt]

        rv -= self.mix.bond_dissociation_activity
        if rv < self.mix.bimolecular_binding_activity:
            # channel is a bimolecular binding
            select = 'bb'
            for bt in self.sig.bond_types:
                if rv < self.mix.activity_bimolecular_binding[bt]:
                    s1, s2 = bt
                    # choose at random (uniformly) an s1, i.e. an agent and free site of required type
                    r1 = self.rng.integers(low=0, high=self.mix.total_free_sites[s1])
                    m1 = self.mix.complexes[self.heap['st'][s1].draw_node(r1)]
                    r2 = self.rng.integers(low=0, high=self.mix.total_free_sites[s2] - m1.free_site[s2])
                    # temporarily modify the specific heap
                    self.heap['st'][s2].modify(m1.free_site[s2] * (m1.count - 1), self.mix.index[m1])
                    # draw the molecule
                    m2 = self.mix.complexes[self.heap['st'][s2].draw_node(r2)]
                    # undo the mod
                    self.heap['st'][s2].modify(m1.free_site[s2] * m1.count, self.mix.index[m1])
                    r1 = self.rng.integers(low=0, high=m1.free_site[s1])
                    r2 = self.rng.integers(low=0, high=m2.free_site[s2])
                    name1, site1 = m1.free_site_list[s1][r1]
                    name2, site2 = m2.free_site_list[s2][r2]
                    self.current_reaction =  select, (m1, m2), (name1, site1), (name2, site2)
                    return
                else:
                    rv -= self.mix.activity_bimolecular_binding[bt]

        rv -= self.mix.bimolecular_binding_activity
        if rv < self.mix.total_inflow:
            # Inflow of atoms. Since there are only a few types, looping is OK.
            select = 'inflow'
            for a in self.mix.activity_inflow:
                if rv < self.mix.activity_inflow[a]:
                    self.current_reaction =  select, (a, None), (None, None), (None, None)
                    return
                else:
                    rv -= self.mix.activity_inflow[a]

        rv -= self.mix.total_inflow
        if rv < self.mix.total_outflow:
            # Outflow of atoms. Since there are only a few types, looping is OK.
            select = 'outflow'
            for a in self.mix.activity_outflow:
                if rv < self.mix.activity_outflow[a]:
                    self.current_reaction =  select, (a, None), (None, None), (None, None)
                    return
                else:
                    rv -= self.mix.activity_outflow[a]

    def report(self, pp_width=40):
        form = '1.5E'
        info = f"\n{'SIMULATOR STATE '.ljust(70, '-')}\n\n"
        # info += f'{"system initialized on:":>{pp_width}} {ka.system.date}\n'
        # info += f'{"system uuid:":>{pp_width}} {ka.system.uuid}\n'
        # info += '\n'
        info += f'{"simulator status at time t=":>{pp_width}} {self.time}\n'
        info += '\n'
        info += f'{"total system activity":>{pp_width}}: {self.mix.total_activity:{form}}\n'
        n_heaps = len(self.heap["bt+"]) + len(self.heap["bt-"]) + len(self.heap["st"])
        if n_heaps > 0:
            info += f'{"heaps":>{pp_width}}: {n_heaps} x ['
            a = next(iter(self.heap))
            b = next(iter(self.heap[a]))
            n_nodes = self.heap[a][b].n_internal_nodes + self.heap[a][b].n_leaves
            info += f"height: {self.heap[a][b].height} nodes: {n_nodes} "
            info += f"occ: {self.heap[a][b].n_entries / self.heap[a][b].n_leaves:.2f}]\n\n"

            # for t in ['bt+', 'bt-', 'st']:
            #     for k in self.heap[t]:
            #         info += f"HEAP {self.heap[t][k].id} -> root value: {self.heap[t][k].tree[0]:1.2E}\n"

            info += f'{"random number generator state":>{pp_width}}\n'
            info += json.dumps(self.rng.bit_generator.state, sort_keys=False, indent=4)
            info += '\n'
        return info

    def __str__(self, pp_width=40):
        return self.report(pp_width=pp_width)

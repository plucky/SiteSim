# Walter Fontana, 2022
""""
This module defines debugging routines.
"""

import kaheap
import kasystem as ka
import kainit


def select_reaction_basic(self):
    """
    Obsolete. Supplanted by the version with heaps. Kept for testing.
    """
    # The first few choices relate to in/outflow of atoms.
    # Since there are only a few types, looping is OK.
    rv = ka.system.sim.rng.uniform(low=0.0, high=ka.system.mixture.total_activity)
    if rv < ka.system.mixture.total_inflow:
        select = 'inflow'
        for a in ka.system.mixture.activity_inflow:
            if rv < ka.system.mixture.activity_inflow[a]:
                self.current_reaction = select, (a, None), (None, None), (None, None)
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
                for m in ka.system.mixture.complexes:
                    segment = m.binding[bt] * m.count
                    if rv < segment:
                        # the event is within molecular species m;
                        # now we uniformly choose the instance of bt in m
                        # bt is (agent_type1.site1), (agent_type2.site2)
                        # NOTE: I don't think we need to randomize which site we choose first.
                        s1, s2 = bt
                        r1 = ka.system.sim.rng.integers(low=0, high=m.free_site[s1])
                        # this is our choice of site1; it belongs to agent name1 in molecule m
                        name1, site1 = m.free_site_list[s1][r1]
                        if s1 == s2:
                            temp_list = [p for p in m.free_site_list[s2] if p != (name1, site1)]
                            r2 = ka.system.sim.rng.integers(low=0, high=m.free_site[s2] - 1)
                            name2, site2 = temp_list[r2]
                        else:
                            # exclude possibility of self-binding
                            site2 = s2.split('.')[1]
                            if (name1, site2) in m.free_site_list[s2]:
                                temp_list = [p for p in m.free_site_list[s2] if p != (name1, site2)]
                                r2 = ka.system.sim.rng.integers(low=0, high=m.free_site[s2] - 1)
                                name2, site2 = temp_list[r2]
                                self.current_reaction =  select, (m, None), (name1, site1), (name2, site2)
                                return
                            else:
                                r2 = ka.system.sim.rng.integers(low=0, high=m.free_site[s2])
                                name2, site2 = m.free_site_list[s2][r2]
                        self.current_reaction =  select, (m, None), (name1, site1), (name2, site2)
                        return
                    else:
                        rv -= segment
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
                for m in ka.system.mixture.complexes:
                    segment = m.unbinding[bt] * m.count
                    if rv < segment:
                        # it's molecule of species m; now we need to uniformly choose the instance of bt in m
                        # bt is (agent_type_1.site_1), (agent_type_2.site_2)
                        # choose a bond of this type
                        r = ka.system.sim.rng.integers(low=0, high=m.bond_type[bt])
                        x, y = m.bond_type_list[bt][r]
                        self.current_reaction =  select, (m, None), x, y
                        return
                    else:
                        rv -= segment
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
                r1 = ka.system.sim.rng.integers(low=0, high=ka.system.mixture.total_free_sites[s1])
                for m1 in ka.system.mixture.complexes:
                    segment = m1.free_site[s1] * m1.count
                    if r1 < segment:
                        # it's m1
                        break
                    else:
                        r1 -= segment
                r2 = ka.system.sim.rng.integers(low=0, high=ka.system.mixture.total_free_sites[s2] - m1.free_site[s2])
                m1.count -= 1  # only temporary!
                for m2 in ka.system.mixture.complexes:
                    segment = m2.free_site[s2] * m2.count
                    if r2 < segment:
                        # it's m2
                        break
                    else:
                        r2 -= segment
                m1.count += 1  # undo
                r1 = ka.system.sim.rng.integers(low=0, high=m1.free_site[s1])
                r2 = ka.system.sim.rng.integers(low=0, high=m2.free_site[s2])
                name1, site1 = m1.free_site_list[s1][r1]
                name2, site2 = m2.free_site_list[s2][r2]
                self.current_reaction =  select, (m1, m2), (name1, site1), (name2, site2)
                return
            else:
                rv -= ka.system.mixture.activity_bimolecular_binding[bt]


def test_heap():
    """
    Heap tester
    """
    system = kainit.initialize()
    print(system.mixture.report())
    var = 'A.p'
    prob = {}
    for m in system.mixture.complexes:
        prob[system.mixture.index[m]] = float(m.free_site[var] * m.count) / float(system.mixture.total_free_sites[var])

    data = [m.free_site[var] * m.count for m in ka.system.mixture.complexes]

    heap = kaheap.Heap(data)
    print(str(heap))

    sample = [0] * len(system.mixture.complexes)
    n_total = 0
    for i in range(0, 1000000):
        rv = system.sim.rng.integers(low=0, high=system.mixture.total_free_sites[var])
        idx = system.heap.draw_node(rv)
        sample[idx] += 1
        n_total += 1

    for i in range(0, len(sample)):
        sample[i] /= float(n_total)
        info = f'{i} prob = {prob[i]}  heap = {sample[i]}'
        print(info)


def check_free_sites():
    totalFS = {}
    for st in ka.system.signature.site_types:
        totalFS[st] = 0
    for m in ka.system.mixture.complexes:
        for st in m.free_site:
            totalFS[st] += (m.free_site[st] * m.count)
    print(f'tallied free sites from <complexes>: {totalFS}')
    print(f'registered total free sites in <mixture>: {ka.system.mixture.total_free_sites}\n')


def check_heap(key):   # key = {'st', 'bt+', bt-'}
    for k in ka.system.sim.heap[key]:
        info = f"heap {key} {k} root value: {ka.system.sim.heap[key][k].tree[0]}"
        print(info)


def test_overall():
    # various test arrangements

    # for i in range(0, 100):
    #     print(system.sim.test_reaction_selection())

    # lprofiler = LineProfiler()
    # lprofiler.add_function(system.sim.select_reaction)
    # lp_wrapper = lprofiler(test)
    # lp_wrapper()
    # lprofiler.print_stats()

    # start_time = time.time()

    system = kainit.initialize()

    print(system.mixture.report())
    print(system.sim.report())
    #
    # check_free_sites()
    # check_heap()
    # ---------------------------------------------------------------
    for i in range(0, 10001):
        if (i % 1000) == 0:
            print(f'\nevent {i+1}')
            print(len(ka.system.mixture.complexes), len(ka.system.mixture.canonical), len(ka.system.mixture.local_views))
            # print(system.sim.report())
            # check_heap()
            # check_free_sites()
        # doubles = False
        # for x in range(0, len(ka.system.mixture.complexes)):
        #     mx = ka.system.mixture.complexes[x]
        #     for y in range(x+1, len(ka.system.mixture.complexes)):
        #         my = ka.system.mixture.complexes[y]
        #         if mx.canonical == my.canonical:
        #             doubles = True
        #             print(f'{x} {y}')
        # if not doubles:
        #     print("All good")

        # reaction = system.sim.select_reaction_basic()
        system.sim.select_reaction()
        system.sim.execute_reaction()
    #
    print(system.mixture.report())
    print(system.sim.report())
    check_free_sites()
    for k in ['st', 'bt+', 'bt-']:
        check_heap(k)


if __name__ == '__main__':
    # various test arrangements

    # test_heap()
    test_overall()

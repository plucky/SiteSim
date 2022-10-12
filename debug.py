# Walter Fontana, 2022
""""
This module defines debugging routines.
"""

import kasim
import kaheap
import kasystem as ka
import kainit

import sys


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

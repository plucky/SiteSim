# Walter Fontana at 10/5/22

import kasystem as ka
import kainit

import time


def main_loop():

    system = kainit.initialize(parameter_file='TestData/parameters_AP.txt')

    # if you read in from a snapshot file, you might want to zero these
    # system.sim.time = 0.
    # system.sim.event = 0

    system.report()
    system.monitor.initialize()
    system.monitor.observe()
    system.monitor.snapshot()

    # The core loop is slightly different for time-based vs event-based observations.
    # In the time-based case, the observation is a non-reactive event, whereas in the
    # event-based case, a reaction event is carried out in addition to the observation.
    # A slight amount of code duplication makes things more readable...

    if ka.system.sim_limit_type == 'time':
        while system.sim.time < system.sim_limit:
            system.sim.advance_time()
            skip = False
            # observation and snapshot are "null reactions"
            if system.sim.time >= system.monitor.observation_time:
                system.sim.time = system.monitor.observation_time
                system.monitor.observe()
                skip = True
            if system.sim.time >= system.monitor.snap_time:
                system.sim.time = system.monitor.snap_time
                system.monitor.snapshot()
                skip = True
            if not skip:
                system.sim.event += 1
                system.sim.select_reaction()
                system.sim.execute_reaction()
    else:
        while system.sim.event < system.sim_limit:
            system.sim.advance_time()
            if system.sim.event == system.monitor.observation_time:
                system.monitor.observe()
            if system.sim.event == system.monitor.snap_time:
                system.monitor.snapshot()
            system.sim.event += 1
            system.sim.select_reaction()
            system.sim.execute_reaction()

    system.monitor.snapshot()
    system.report()

    print("\nDone!")
    print(f'events: {system.sim.event}')
    print(f'cpu: {time.process_time()}\n')


if __name__ == '__main__':
    main_loop()

#!/opt/local/bin/python

# Walter Fontana at 10/5/22

import kasystem as ka
import kainit

import time


def main_loop():

    system = kainit.initialize(parameter_file='TestData/parameters_AP.txt')
    simulator = system.sim
    monitor = system.monitor
    
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

    if system.sim_limit_type == 'time':
        while simulator.time < system.sim_limit:
            simulator.advance_time()
            skip = False
            # observation and snapshot are "null reactions"
            if simulator.time >= monitor.observation_time:
                simulator.time = monitor.observation_time
                monitor.observe()
                skip = True
            if simulator.time >= monitor.snap_time:
                simulator.time = monitor.snap_time
                monitor.snapshot()
                print(f'events: {simulator.event}')
                print(f'cpu: {time.process_time()}\n')
                skip = True
            if not skip:
                simulator.event += 1
                simulator.select_reaction()
                # print(simulator.event)
                simulator.execute_reaction()
    else:
        while simulator.event < system.sim_limit:
            simulator.advance_time()
            if simulator.event == monitor.observation_time:
                monitor.observe()
            if simulator.event == monitor.snap_time:
                monitor.snapshot()
            simulator.event += 1
            simulator.select_reaction()
            simulator.execute_reaction()

    monitor.snapshot()
    system.report()

    print("\nDone!")
    print(f'events: {simulator.event}')
    print(f'cpu: {time.process_time()}\n')


if __name__ == '__main__':
    main_loop()

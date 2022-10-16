# Walter Fontana at 10/5/22

import kasystem as ka
import kainit
import kamon

import time


def main_loop():

    system = kainit.initialize(parameter_file='TestData/parameters_AP.txt')

    # if you read in from a snapshot file, you might want to zero these
    system.sim.time = 0.
    system.sim.event = 0

    system.report()

    observation_alert = 0
    monitor = kamon.Monitor()

    # The core loop is slightly different for time-based vs event-based observations.
    # In the time-based case, the observation is a non-reactive event, whereas in the
    # event-based case, a reaction event is carried out in addition to the observation.
    # A slight amount of code duplication makes things more readable...

    monitor.snapshot()
    # monitor.observe()

    if ka.system.sim_limit_type == 'time':
        while system.sim.time < system.sim_limit:
            system.sim.advance_time()
            if system.sim.time >= observation_alert:
                system.sim.time = observation_alert
                monitor.observe()
                observation_alert += system.obs_freq
            else:
                system.sim.event += 1
                system.sim.select_reaction()
                system.sim.execute_reaction()
    else:
        while system.sim.event < system.sim_limit:
            system.sim.advance_time()
            if system.sim.event == observation_alert:
                monitor.observe()
                observation_alert += system.obs_freq
            system.sim.event += 1
            system.sim.select_reaction()
            system.sim.execute_reaction()

    monitor.snapshot()
    system.report()

    print("\nDone!")
    print(f'events: {system.sim.event}')
    print(f'cpu: {time.process_time()}\n')


if __name__ == '__main__':
    main_loop()

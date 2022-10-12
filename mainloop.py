# Walter Fontana at 10/5/22

import kasystem as ka
import kainit
import kamon


def main_loop():

    system = kainit.initialize(parameter_file='TestData/parameters_test2.txt')

    # if you read in from a snapshot file, you might want to zero these
    system.sim.time = 0.
    system.sim.event = 0

    print(system.mixture.report())
    print(system.sim.report())

    observation_alert = 0
    kamon.set_up_monitoring()

    # The core loop is slightly different for time-based vs event-based observations.
    # In the time-based case, the observation is a non-reactive event, whereas in the
    # event-based case, a reaction event is carried out in addition to the observation.
    # A slight amount of code duplication makes things more readable...

    if ka.system.sim_limit_type == 'time':
        while system.sim.time < system.sim_limit:
            # print(system.sim.time)
            if system.sim.time >= observation_alert:
                system.sim.time = observation_alert
                kamon.monitor()
                observation_alert += system.obs_freq
            else:
                system.sim.advance_time()
                system.sim.select_reaction()
                system.sim.execute_reaction()
                system.sim.event += 1
    else:
        while system.sim.event < system.sim_limit:
            if system.sim.event == observation_alert:
                ka.system.monitor()
                observation_alert += system.obs_freq
            system.sim.advance_time()
            system.sim.select_reaction()
            system.sim.execute_reaction()
            system.sim.event += 1

    print(system.mixture.report())
    print(system.sim.report())


if __name__ == '__main__':
    main_loop()

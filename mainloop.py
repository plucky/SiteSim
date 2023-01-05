#!/opt/local/bin/python

# Walter Fontana, 2022

import kainit
import datetime

def in_notebook():
    """
    Determine whether we are in a notebook environment
    """
    try:
        from IPython import get_ipython
        if 'IPKernelApp' not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True


def initialize(parameter_file=None, invocation= None, modifier_fun=None, **kwargs):
    """
    Wrapper for initialization
    """
    if not in_notebook():
        kainit.commandline(invocation=invocation)
    # note that parameter_file may be overridden by the command line argument
    system = kainit.initialize(parameter_file=parameter_file, modifier_fun=modifier_fun, **kwargs)

    return system


def loop(system=None):
    """
    Main simulator loop
    """
    simulator = system.sim
    monitor = system.monitor
    
    # if you read in from a snapshot file, you might want to zero these
    # system.sim.time = 0.
    # system.sim.event = 0

    system.report()
    # system.monitor.initialize(file=parameter_file)
    system.monitor.observe()
    system.monitor.snapshot(flag='first')

    local_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f'\nSimulation <{system.uuid}> started at {local_time}')

    # ====================================================================================================
    # The core loop is slightly different for time-based vs event-based observations.
    # In the time-based case, the observation is a non-reactive event, whereas in the
    # event-based case, a reaction event is carried out in addition to the observation.
    # A slight amount of code duplication makes things more readable...

    if system.sim_limit_type == 'time':
        while simulator.time < system.sim_limit:
            simulator.advance_time()
            skip = False
            # future: add time-specific interventions here...
            if simulator.time >= monitor.observation_time:
                simulator.time = monitor.observation_time
                monitor.observe()
                # check for stopping conditions
                if system.alarm.trigger():
                    # a stopping condition was triggered
                    break
                # interim report
                system.report()
                # an observation (or snapshot or intervention) at a specified time is a "null reaction"
                skip = True
            if simulator.time >= monitor.snap_time:
                simulator.time = monitor.snap_time
                monitor.snapshot()
                skip = True
            # if we had an observation, skip the reaction
            if not skip:
                simulator.event += 1
                simulator.select_reaction()
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

    # ====================================================================================================

    monitor.snapshot(flag='last')
    # final reports
    system.report()
    system.resources_report()

    local_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f'Simulation <{system.uuid}> terminated at {local_time}\n')


def SiteSim_loop(parameter_file='TestData/parameters_AP.txt'):
    """
    Plain simulation loop
    """
    system = initialize(parameter_file=parameter_file, invocation=None)
    loop(system=system)
    # clear objects
    kainit.clear_SiteSim()


if __name__ == '__main__':
    SiteSim_loop()

#!/opt/local/bin/python

# Walter Fontana, 2022

import kainit


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


def SiteSim_loop(parameter_file='TestData/parameters.txt', modifier_fun=None):

    if in_notebook():
        # avoids sys.argv (command line parsing)
        system = kainit.init(parameter_file=parameter_file, parameter_mod_fun=modifier_fun)
    else:
        # note that parameter_file may be overridden by command line argument
        system = kainit.initialize(parameter_file=parameter_file, parameter_mod_fun=modifier_fun)

    simulator = system.sim
    monitor = system.monitor
    
    # if you read in from a snapshot file, you might want to zero these
    # system.sim.time = 0.
    # system.sim.event = 0

    system.report()
    # system.monitor.initialize(file=parameter_file)
    system.monitor.observe()
    system.monitor.snapshot(flag='first')

    print(f'\nSimulation <{system.uuid}> started')

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
                # print(f'events: {simulator.event}')
                # print(f'cpu: {time.process_time()}\n')
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
    # clear objects
    kainit.clear_SiteSim()

    print(f'Simulation <{system.uuid}> terminated')


if __name__ == '__main__':
    SiteSim_loop()

# Here is an example of an external parameter modifying function.
#
# def modify_parameters(system):
#
#     # file names
#     system.report_file = report_file
#     system.monitor.obs_file_name = csv_datafile
#     system.monitor.snap_root_name = snap_root
#
#     # physical parameters
#     system.parameters.RescaleTemperature = 1.
#     system.parameters.ResizeVolume = 0.1
#
#     # initial agent counts
#     system.signature.init_agents['P'] = 100
#
#     print(f'parameters updated')

# To loop over input parameters, you have to inline the initialization (kainit.init) into a version of SiteSim_loop,
# invoke the parameter modifications instead of calling the parameter modifying function. Then loop over everything
# varying the values of parameters.
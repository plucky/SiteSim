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


def initialize(parameter_file=None, modifier_fun=None, **kwargs):
    """
    Wrapper for initialization
    """
    if in_notebook():
        # avoids sys.argv (command line parsing)
        system = kainit.init(parameter_file=parameter_file, modifier_fun=modifier_fun, **kwargs)
    else:
        # note that parameter_file may be overridden by command line argument
        system = kainit.initialize(parameter_file=parameter_file, modifier_fun=modifier_fun, **kwargs)

    return system


def modify_parameters(system, **kwargs):
    for key, value in kwargs.items():
        match key:
            # file names
            case 'report_fn':
                system.report_file = value
                print(f'"report_file" updated to {system.report_file}')
            case 'output_fn':
                system.monitor.obs_file_name = value
                print(f'"csv_datafile" updated to {system.monitor.obs_file_name}')
            case 'snap_root':
                system.monitor.snap_root_name = value
                print(f'"snap_root" updated to {system.monitor.snap_root_name}')
            # physical parameters
            case 'RescaleTemperature':
                system.parameters.RescaleTemperature = value
                print(f'"RescaleTemperature" updated to {system.parameters.RescaleTemperature}')
            case 'Temperature':  # in C
                system.parameters.Temperature = float(value) + 273.15  # in K
                print(f'"Temperature" updated to {value} ºC ({system.parameters.Temperature} K)')
            case 'ResizeVolume':
                system.parameters.ResizeVolume = value
                print(f'"ResizeVolume" updated to {system.parameters.ResizeVolume}')
            # initial agent counts
            case 'initial_agent_counts':
                for agent_type, count in value.items():
                    system.signature.init_agents[agent_type] = count
                    print(f'Initial agent count for {agent_type} updated to {count}')


def main_loop(system=None):
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

    print(f'Simulation <{system.uuid}> terminated\n')


def SiteSim_loop(parameter_file='TestData/parameters_AP.txt'):
    """
    Plain simulation loop
    """
    system = initialize(parameter_file=parameter_file)
    main_loop(system=system)
    # clear objects
    kainit.clear_SiteSim()


def pd_loop(parameter_file='TestData/parameters_AP.txt'):
    """
    Simulation loop with modification of parameters
    """
    # modify below as appropriate
    for i in range(0, 11):
        temp = 25 - i * 2
        kwargs = dict()
        kwargs['report_fn'] = f'report_T{temp}.txt'
        kwargs['output_fn'] = f'output_T{temp}.csv'
        kwargs['snap_root'] = f'snap_T{temp}_'
        kwargs['Temperature'] = float(temp)
        kwargs['ResizeVolume'] = 0.1
        # kwargs['initial_agent_counts'] = {'A': 100, 'P': 100}

        system = initialize(parameter_file=parameter_file, modifier_fun=modify_parameters, **kwargs)
        main_loop(system=system)
        # clear objects
        kainit.clear_SiteSim()


if __name__ == '__main__':
    # SiteSim_loop()
    pd_loop()

#!/opt/local/bin/python

# Walter Fontana, 2023

import mainloop as main
import kainit
from pathlib import Path
import sys

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
            case 'referenceRingClosureFactor':
                system.parameters.referenceRingClosureFactor = value
                print(f'"referenceRingClosureFactor" updated to {system.parameters.referenceRingClosureFactor}')
            case 'initial_agent_counts':
                for agent_type, count in value.items():
                    system.signature.init_agents[agent_type] = count
                    print(f'Initial agent count for {agent_type} updated to {count}')

# ========================================================================================================
# Specialized loops
#
def TC_loop():
    """
    Simulation loop with modification of temperature and concentration of one agent type
    """
    # defaults
    dir_root = './'
    T_start = 25
    T_delta = 0
    agent_type = ''
    C_start = 100
    C_delta = 0
    T_loop_start = 0
    T_loop_end = 1
    C_loop_start = 0
    C_loop_end = 1

    system = kainit.commandline(invocation=None)
    # save the parameter file name from the command line for repeated initialization below
    parameter_file = system.parameter_file

    # process "extra" arguments
    if system.xargs:
        for key, value in system.xargs.arg_dict.items():
            if key == 'T':  # temperature scan
                for item in value:
                    var, val = item.split('=')
                    if var == 'initial':
                        T_start = float(val)
                    elif var == 'delta':
                        T_delta = float(val)
                    elif var == 'start':
                        T_loop_start = int(val)
                    elif var == 'end':
                        T_loop_end = int(val)
                    else:
                        sys.exit(f"X argument: Temperature parameter {var} not recognized.")
            elif key == 'C':  # concentration scan
                for item in value:
                    var, val = item.split('=')
                    if var == 'agent':
                        agent_type = val
                    elif var == 'initial':
                        C_start = float(val)
                    elif var == 'delta':
                        C_delta = float(val)
                    elif var == 'start':
                        C_loop_start = int(val)
                    elif var == 'end':
                        C_loop_end = int(val)
                    else:
                        sys.exit(f"X argument: Concentration parameter {var} not recognized.")
            elif key == 'dir':
                dir_root = value[0].strip()
                if dir_root[-1] != '/':
                    dir_root += '/'

    # preview
    print("These directories will be eventually created:")
    for t in range(T_loop_start, T_loop_end):
        temp = T_start + t * T_delta
        for c in range(C_loop_start, C_loop_end):
            concentration = C_start + c * C_delta
            if agent_type:
                run_name = f'T{temp}_{agent_type}{concentration}'
            else:
                run_name = f'T{temp}'
            dir_name = dir_root + run_name
            print(dir_name)
    ans = input("Continue? [N/y] ")
    if ans != 'y':
        sys.exit("Exiting.")

    # ====================================================================================================
    mod_args = {}
    agent_dict = {}
    for t in range(T_loop_start, T_loop_end):
        temp = T_start + t * T_delta
        for c in range(C_loop_start, C_loop_end):
            concentration = C_start + c * C_delta

            if agent_type:
                agent_dict[agent_type] = concentration
                run_name = f'T{temp}_{agent_type}{concentration}'
            else:
                run_name = f'T{temp}'

            mod_args['initial_agent_counts'] = agent_dict

            # make directory
            dir_name = dir_root + run_name
            try:
                Path(dir_name).mkdir(parents=True, exist_ok=False)
            except FileExistsError:
                while True:
                    print(f'directory {dir_name} already exists. Trying to add identifier.')
                    dir_split = dir_name.split('--')
                    if len(dir_split) == 1:
                        dir_name = dir_split[0] + f'--{1}'
                    else:
                        n = dir_split[-1]
                        if n.isdigit():
                            i = int(n)
                            if i > 100:
                                sys.exit(f'Refuse to add identifier to {dir_name}. Exiting.')
                        else:
                            i = 0
                        i = str(i + 1)
                        dir_name = ''.join(dir_split[:-1]) + f'--{i}'
                    try:
                        Path(dir_name).mkdir(parents=True, exist_ok=False)
                        break
                    except FileExistsError:
                        continue

            mod_args['report_fn'] = dir_name + f'/report_{run_name}.txt'
            mod_args['output_fn'] = dir_name + f'/output_{run_name}.csv'
            mod_args['snap_root'] = dir_name + f'/snap_{run_name}_'
            mod_args['Temperature'] = temp
            # mod_args['ResizeVolume'] = 0.1
            # mod_args['RescaleTemperature'] = 1.
            # mod_args['referenceRingClosureFactor'] = 1.e+5

            system = kainit.initialize(parameter_file=parameter_file, modifier_fun=modify_parameters, **mod_args)
            main.loop(system=system)

            # clear objects
            kainit.clear_SiteSim()


if __name__ == '__main__':
    TC_loop()

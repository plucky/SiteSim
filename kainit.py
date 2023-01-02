# Walter Fontana, 2022
"""
This module initializes the system.
"""
import os
import sys
import re
import argparse

import kasystem as ka
import kaparam
import kamol
import kamix
import kasim
import kamatch
import kamon
import kalarm


def test_invocation():
    """
    construct a commandline invocation for testing purposes

    Returns: a commandline invocation to be parsed by commandline()
    """
    test_line = []
    # test_line += ['-s', 'A(p[a1.P$w a2.P a3.P], l[r.A], r[l.A]), P@100(a1[p.A], a2[p.A], a3[p.A], d[d.P])']
    # test_line += ['-s', 'A(p[a1.P$m a2.P$m a3.P$m], l[r.A$w], r[l.A]), P(a1[p.A], a2[p.A], a3[p.A], d[d.P$m])']
    # test_line += ['--sig_file', 'TestData/signature.txt']
    test_line += ['-p', 'TestData/parameters_test2.txt']
    # test_line += ['-r', 'TestOutput/report.txt']
    # test_line += ['-m', 'TestData/snap__big.ka']
    # test_line += ['-m', 'TestData/snap__0015.ka']
    # test_line += ['-d', '0']
    # test_line += ['--seed', '1234']

    return test_line


def commandline(invocation=None):
    """
    parse the command line

    Returns: a pair consisting of (1) a dictionary {argument: value} and (2) the command line used
    """
    # define the command line arguments
    cmd = argparse.ArgumentParser(description='A simulator for arbitrary context-less binding interactions')
    cmd.add_argument('-s', '--signature', type=str,
                     help='signature of the interaction system')
    cmd.add_argument('-p', '--parameters', type=str, default='parameters.txt',
                     help='parameter file name')
    cmd.add_argument('-r', '--report', type=str, default='report.txt',
                     help='report file name')
    cmd.add_argument('-m', '--mixture', type=str, default=None,
                     help='initial mixture file')
    cmd.add_argument('--seed', type=int, default=-1,
                     help='random number seed')
    cmd.add_argument('-d', '--db', type=int, default=0,
                     help='reporting level')

    shortcuts = {'s': 'signature', 'r': 'report', 'p': 'parameters', 'm': 'mixture', 'd': 'db'}

    args = cmd.parse_args(invocation)

    arg_values = vars(args)  # command line {argument: value} dictionary

    # reconstruct the command line used so that it can be copy/pasted and rerun
    cmdline = os.path.split(sys.argv[0])[1]  # this is the program name
    for item in sys.argv[1:]:
        if re.match('--', item):
            key = item[2:]
            val = arg_values[key]
            if type(val) is str:
                cmdline += f' --{key} "{val}"'
            else:
                cmdline += f' --{key} {val}'
        elif re.match('-', item):
            key = item[1:]
            val = arg_values[shortcuts[key]]
            if type(val) is str:
                cmdline += f' -{key} "{val}"'
            else:
                cmdline += f' -{key} {val}'
        elif item in arg_values:
            val = arg_values[item]
            if type(val) is str:
                cmdline += f' {item} "{val}"'
            else:
                cmdline += f' {item} {val}'

    return arg_values, cmdline


def initialize(parameter_file=None, invocation=None, modifier_fun=None, **kwargs):
    """
    Initializes the system; includes commandline invocation
    Returns: a system object
    """
    ka.init_system()  # this creates, but does not fully initialize, the global object "system"

    ka.system.parameter_file = parameter_file  # will be overridden by command line, if any
    seed = -1  # flag to determine whether the command line contains a seed

    if len(sys.argv) > 1 or invocation is not None:
        # use command line or pass invocation
        arg_values, cmdline = commandline(invocation)  # returns {argument: value} dictionary and the invocation

        # store the COMMAND LINE
        ka.system.cmdline = cmdline

        for key in arg_values:

            # SIGNATURE
            if key == 'signature':
                ka.system.signature_string = arg_values[key]

            # PARAMETER file name
            elif key == 'parameters':
                ka.system.parameter_file = arg_values[key]

            # REPORT file name
            elif key == 'report':
                ka.system.report_file = arg_values[key]

            # MIXTURE file name
            elif key == 'mixture':
                ka.system.mixture_file = arg_values[key]

            # REPORTING level
            elif key == 'db':
                ka.system.db_level = arg_values[key]

            # random number SEED
            elif key == 'seed':
                seed = arg_values[key]

    # create monitor (needs further initializing before starting the main loop)
    ka.system.monitor = kamon.Monitor()

    # initialize Kappa PARSER
    ka.system.kappa = kamol.Kappa()

    # initialize site-graph MATCHER
    ka.system.sgm = kamatch.SiteGraphMatcher()

    # initialize PARAMETERS. This also initializes the SIGNATURE object.
    ka.system.parameters = kaparam.Parameters(file=ka.system.parameter_file)

    # change parameters using an external function
    if modifier_fun:
        modifier_fun(ka.system, **kwargs)

    # apply settings and compute rate constants
    ka.system.parameters.apply_parameters()

    # use the rng seed provided on the commandline, if any
    if seed != -1:
        ka.system.parameters.rng_seed = seed

    # initialize MIXTURE, which also generates all activities
    ka.system.mixture = kamix.Mixture(file=ka.system.mixture_file, system=ka.system)

    # initialize CTMC
    ka.system.sim = kasim.CTMC()

    # initialize MONITOR
    ka.system.monitor.initialize(file=ka.system.parameter_file)

    # initialize ALARMS
    ka.system.alarm = kalarm.Alarm()

    # initial REPORT
    ka.system.report()

    # ka.system.snapshot("TestOutput/snap0.ka", label=True)

    return ka.system


def init(parameter_file=None, modifier_fun=None, **kwargs):
    """
    Initializes a system without commandline and allows for parameter modifications.
    This is useful for scripting in jupyter.
    Returns: a system object
    """
    ka.init_system()  # this creates, but does not fully initialize, the global object "system"
    ka.system.parameter_file = parameter_file
    # create monitor (needs further initializing before starting the main loop)
    ka.system.monitor = kamon.Monitor()
    # initialize Kappa PARSER
    ka.system.kappa = kamol.Kappa()
    # initialize site-graph MATCHER
    ka.system.sgm = kamatch.SiteGraphMatcher()
    # initialize PARAMETERS. This also initializes the SIGNATURE object.
    ka.system.parameters = kaparam.Parameters(file=ka.system.parameter_file)
    # initialize MIXTURE, which also generates all activities

    # change parameters using an external function
    if modifier_fun:
        modifier_fun(ka.system, **kwargs)

    # apply settings and compute rate constants
    ka.system.parameters.apply_parameters()

    ka.system.mixture = kamix.Mixture(file=ka.system.mixture_file, system=ka.system)
    # initialize CTMC
    ka.system.sim = kasim.CTMC()
    # initialize MONITOR
    ka.system.monitor.initialize(file=ka.system.parameter_file)
    # initialize ALARMS
    ka.system.alarm = kalarm.Alarm()
    # initial REPORT
    ka.system.report()

    return ka.system


def clear_SiteSim():
    del ka.system.kappa
    del ka.system.sgm
    del ka.system.signature
    del ka.system.parameters
    del ka.system.mixture
    del ka.system.sim
    del ka.system.monitor
    del ka.system.alarm
    del ka.system


if __name__ == '__main__':
    pass

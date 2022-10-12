# Walter Fontana, 2022
"""
This module initializes the system.
"""
import os
import sys
import re
import argparse

import kasig
import kasystem as ka
import kaparam
import kamol
import kamix
import kasim
import kamatch


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


def initialize(parameter_file=None, invocation=None):
    """
    Initializes the system

    Returns: a system object
    """
    ka.init_system()  # this creates, but does not fully initialize, the global object "system"

    ka.system.parameter_file = parameter_file
    seed = -1

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

    # initialize PARAMETERS
    ka.system.parameters = kaparam.Parameters(file=ka.system.parameter_file)
    # use the rng seed provided on the commandline, if any
    if seed != -1:
        ka.system.parameters.rng_seed = seed

    # initialize Kappa PARSER
    ka.system.kappa = kamol.Kappa()

    # initialize site-graph matcher
    ka.system.sgm = kamatch.SiteGraphMatcher()

    # initialize MIXTURE, which also generates all activities
    ka.system.mixture = kamix.Mixture(file=ka.system.mixture_file, system=ka.system)

    # initialize CTMC and finish initializing SYSTEM (dependency on CTMC)
    ka.system.sim = kasim.CTMC()

    # initial REPORT
    ka.system.reporter()

    # ka.system.snapshot("TestOutput/snap0.ka", label=True)

    return ka.system


def minit(signature=None, parameter_file=None):
    """
    Initializes a minimal system.
    Returns: a system object
    """

    if not signature and not parameter_file:
        sys.exit("No signature or parameter file given.")

    ka.init_system()  # this creates, but does not fully initialize, the global object "system"
    # store the COMMAND LINE
    ka.system.cmdline = 'minit invocation'
    ka.system.signature_string = signature
    ka.system.parameters = kaparam.Parameters(parameter_file)
    ka.system.parameters.rng_seed = 42
    ka.system.kappa = kamol.Kappa()

    return ka.system


if __name__ == '__main__':
    pass

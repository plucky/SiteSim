# Walter Fontana, 2022

import re
import os
import sys
import kasystem as ka


class Watermark:  # not used for now
    def __init__(self):
        sorted_complexes = sorted(ka.system.mixture.complexes, key=lambda x: x.size, reverse=True)
        self.level = sorted_complexes[0].size
        self.threshold = 0

    def check(self):
        if self.level > self.threshold:
            return True
        else:
            return False

class Observable:
    def __init__(self, name, value, index=0, threshold=0):
        self.name = name
        self.value = value  # this should be the list of deques with values
        self.index = index
        self.threshold = threshold

    def check(self):
        # this is the most recently observed value
        if self.value[self.index][-1] > self.threshold:
            return True
        else:
            return False

class Alarm:
    """
    Manage and check stopping conditions.
    """
    def __init__(self):
        self.alarm = {}

        self.read_alarms()

    def trigger(self):
        """
        Check if stopping conditions are satisfied.
        """
        for alarm in self.alarm:
            if self.alarm[alarm].check():
                print(f'{alarm} hit threshold at time {ka.system.sim.time:.5f} (event {ka.system.sim.event})')
                return True
        return False

    def read_alarms(self):
        """
        Scan the parameter file for stopping conditions.
        """
        if not os.path.isfile(ka.system.parameter_file):
            sys.exit("Cannot find parameter file %s" % ka.system.parameter_file)
        else:
            with open(ka.system.parameter_file, "r", encoding='utf-8') as data:
                while True:
                    line = data.readline()
                    if not line:
                        break
                    # parse the line
                    match = re.match(r'^%stp:\s*(.*)\s*>\s*([0-9]+)\s?', line)
                    if match:
                        name = match.group(1).strip()
                        threshold = match.group(2).strip()
                        if name == 'size_watermark':  # not used for now
                            self.alarm['size_watermark'] = Watermark()
                            self.alarm['size_watermark'].threshold = int(threshold)
                        else:
                            match = re.match(r'(.*)\[(.*)\]$', name)
                            index = 0
                            if match:
                                _name = match.group(1)
                                index = int(match.group(2))
                            else:
                                _name = name
                            if _name in ka.system.monitor.observable_by_name:
                                # the next() construction is because a name has a unique obs_type, so
                                # no need to search for it.
                                obs_type = next(iter(ka.system.monitor.observable_by_name[_name]))
                                value = ka.system.monitor.observable_by_name[_name][obs_type]['value']
                                self.alarm[name] = Observable(name, value, index=index, threshold=int(threshold))

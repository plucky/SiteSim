# Walter Fontana, 2022

import re
import os
import sys
import kasystem as ka


class Watermark:
    def __init__(self):
        sorted_complexes = sorted(ka.system.mixture.complexes, key=lambda x: x.size, reverse=True)
        self.level = sorted_complexes[0].size
        self.threshold = 0

    def check(self):
        if self.level > self.threshold:
            return True
        else:
            return False


class Alarm:
    """
    Manage and check stopping conditions.
    """
    def __init__(self):
        # when scripting, external functions can be appended to self.fuse
        self.fuse = []
        self.watermark = Watermark()

        self.read_alarms()

        if self.watermark.threshold != 0:
            self.fuse.append(self.watermark.check)

    def check(self):
        """
        Check if stopping conditions are satisfied.
        """
        for fun in self.fuse:
            if fun():
                if fun == self.watermark.check:
                    s = f'{fun}'.split()[2]
                    print(f'{s} hit threshold')
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
                        value = match.group(2).strip()
                        if name == 'size_watermark':
                            self.watermark.threshold = int(value)

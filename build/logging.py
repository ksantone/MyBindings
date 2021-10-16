#!/usr/bin/env python

import denovo_sequencing
import sys

class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    
    def flush(self):
        pass

sys.stdout = Logger("/Users/kassimsantone/Desktop/fragment.txt")

denovo_sequencing.main()

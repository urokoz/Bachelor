#!/usr/bin/env python
import sys

infile = open(sys.argv[1], "r")

for i, line in enumerate(infile):
    sys.stdout.write(">seq_{:06d}\n".format(i) + line)

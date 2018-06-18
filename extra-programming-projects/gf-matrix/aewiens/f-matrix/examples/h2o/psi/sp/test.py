#!/usr/bin/env python

import re

f = open("output.dat","r").read()

test = r"\s+Total Energy\s+\=\s+(\-\d+\.\d+)"

print(re.findall(test,f))

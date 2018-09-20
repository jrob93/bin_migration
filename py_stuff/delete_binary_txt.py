#!/usr/bin/python
import subprocess
import os

dirs=next(os.walk('.'))[1]
dirs = [d for d in dirs if 'prone' in d]

for d in dirs:
    print subprocess.check_output(["rm {}/binary*txt".format(d)],shell=True)

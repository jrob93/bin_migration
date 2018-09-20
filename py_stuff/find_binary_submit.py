'''
read submit file, and see how the binaries were input
'''
import numpy
import os
from natsort import natsorted, ns
import subprocess
# import py_func as pf
import pandas as pd

path='/Volumes/Mirkwood/bin_migration/data'
dirs=next(os.walk(path))[1]
dirs = [d for d in dirs if 'prone' in d]
dirs=natsorted(dirs, key=lambda y: y.lower()) # sort unpadded numbers

print dirs

of=open('find_binary_submit.txt','w')
of.write('{}\t{}\t{}\n'.format('run','binary_file','input_file'))

for run in dirs:
    with open("{}/{}/submit.{}".format(path,run,run.split('_')[1])) as f:
        _dat = [line.rstrip() for line in f]
    args = _dat[9].split()
    print run,args[-3],args[-1]
    of.write('{}\t{}\t{}\n'.format(run,args[-3],args[-1]))

of.close()

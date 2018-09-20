#!/usr/bin/python
import os

files=next(os.walk('.'))[2]
files = [ fi for fi in files if fi.endswith(".tar")]
files.sort()

# files=["prone_1006.tar","prone_1007.tar","prone_1008.tar","prone_1009.tar","prone_1010.tar",
# "prone_1011.tar","prone_1012.tar","prone_1013.tar","prone_1014.tar","prone_1015.tar",
# "prone_1016.tar","prone_1017.tar","prone_1018.tar","prone_1019.tar"]

dirs=next(os.walk('.'))[1]

for f in files:
    d=f.split('.')[0]
    n=d.split('_')[-1]
    print f,n,d
    if (d in dirs) or (int(n)<1193):
        print "already unpacked"
        continue
    else:
        os.system ("tar -xvf {}".format(f)) # untar the main directory
        os.system("tar -xvf {}/{}_dump.tgz --directory {}".format(d,n,d)) # untar the dump files
    # break

import numpy
import matplotlib.pyplot as pyplot
import os
from natsort import natsorted, ns
import subprocess
import matplotlib.gridspec as gridspec
import sys
sys.path.insert(0, '../../../grav_cloud/python_stuff/')
import py_func as pf

path='/Volumes/Mirkwood/bin_migration/data'
dirs=["prone_1006","prone_1007","prone_1008","prone_1009","prone_1010",
"prone_1011","prone_1012","prone_1013","prone_1014","prone_1015",
"prone_1016","prone_1017","prone_1018","prone_1019"]#,"prone_1020",
# "prone_1021","prone_1022","prone_1023","prone_1024","prone_1025",
# "prone_1026","prone_1027","prone_1028","prone_1029","prone_1030",
# "prone_1031","prone_1032","prone_1033","prone_1034","prone_1035",
# "prone_1036","prone_1037","prone_1038","prone_1039","prone_1040",
# "prone_1041","prone_1042","prone_1043","prone_1044","prone_1045",
# "prone_1046","prone_1047","prone_1048","prone_1049","prone_1050",
# "prone_1051","prone_1052","prone_1053","prone_1054","prone_1055",
# "prone_1056","prone_1057","prone_1058","prone_1059","prone_1060",
# "prone_1061","prone_1062","prone_1063","prone_1064","prone_1065",
# "prone_1066","prone_1067","prone_1068","prone_1069","prone_1070",
# "prone_1071","prone_1072","prone_1073","prone_1074","prone_1075",
# "prone_1076","prone_1077","prone_1078","prone_1079","prone_1080",
# "prone_1081","prone_1082","prone_1083","prone_1084","prone_1085",
# "prone_1086","prone_1087","prone_1088","prone_1089","prone_1090",
# "prone_1091","prone_1092","prone_1093","prone_1094","prone_1095",
# "prone_1096","prone_1097","prone_1098","prone_1099","prone_1100"]

dirs=next(os.walk(path))[1]
# dirs=['prone_1335_timesteptest0.75']
# dirs=['prone_1335']
dirs=['prone_1000_timesteptest0.25','prone_1000_timesteptest0.75',
'prone_1040_timesteptest0.25','prone_1040_timesteptest0.75',
'prone_1335_timesteptest0.25','prone_1335_timesteptest0.75']
# dirs=['prone_1040_timesteptest0.25','prone_1335_timesteptest0.25']

# dirs=["prone_122","prone_123","prone_124","prone_125","prone_126",
# "prone_127","prone_128","prone_129","prone_130","prone_131","prone_132","prone_133"]
N_planet=4

for run in dirs:

    # if int(run.split('_')[1])<1006:
    #     continue

    print run
    subprocess.Popen(["mkdir","../{}_analysis".format(run)]) # make directory to store images

    run_path='{}/{}'.format(path,run)
    files=next(os.walk(run_path))[2]
    files = [ fi for fi in files if fi.endswith('.dump')]
    files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers

    time=[]
    a_planet=numpy.zeros((len(files),N_planet))
    e_planet=numpy.zeros((len(files),N_planet))
    I_planet=numpy.zeros((len(files),N_planet))
    for i in range(len(files)):
        filename='{}/{}'.format(run_path,files[i])
        print filename
        with open(filename) as f:
            dat = [line.rstrip() for line in f]
        dat = [line.split() for line in dat]
        print dat[0] # header line
        t=float(dat[0][3]) # timestamp
        dat = dat[1:] # strip header line
        dat = numpy.array(dat,dtype=float) # convert data list into array
        time.append(t)
        # data format:
        # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0

        for j in range(1,N_planet+1):
            planet_record_file="../{}_analysis/planet_{:03d}.txt".format(run,j)
            if i==0: # wipe file to avoid appending errors
                with open(planet_record_file, 'w') as file:
                    file.write("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n".format(t,dat[j,1],dat[j,2],dat[j,3],dat[j,4],dat[j,5],dat[j,6]))
            else:
                with open(planet_record_file, 'a') as file:
                    file.write("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n".format(t,dat[j,1],dat[j,2],dat[j,3],dat[j,4],dat[j,5],dat[j,6]))

        a_planet[i,:]=dat[1:N_planet+1,1]
        e_planet[i,:]=dat[1:N_planet+1,2]
        I_planet[i,:]=dat[1:N_planet+1,3]

    # # LOWER RESOLUTION
    # time=[]
    # print len(files)
    # file_nums=range(0,len(files),1000)
    # print file_nums
    # a_planet=numpy.zeros((len(file_nums),N_planet))
    # e_planet=numpy.zeros((len(file_nums),N_planet))
    # I_planet=numpy.zeros((len(file_nums),N_planet))
    # j=0
    # for i in file_nums:
    #     filename='{}/{}'.format(run_path,files[i])
    #     print filename
    #     with open(filename) as f:
    #         dat = [line.rstrip() for line in f]
    #     dat = [line.split() for line in dat]
    #     print dat[0] # header line
    #     t=float(dat[0][3]) # timestamp
    #     dat = dat[1:] # strip header line
    #     dat = numpy.array(dat,dtype=float) # convert data list into array
    #     time.append(t)
    #     # data format:
    #     # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0
    #
    #     for _j in range(1,N_planet+1):
    #         planet_record_file="../{}_analysis/planet_{:03d}.txt".format(run,_j)
    #         if i==0: # wipe file to avoid appending errors
    #             with open(planet_record_file, 'w') as file:
    #                 file.write("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n".format(t,dat[_j,1],dat[_j,2],dat[_j,3],dat[_j,4],dat[_j,5],dat[_j,6]))
    #         else:
    #             with open(planet_record_file, 'a') as file:
    #                 file.write("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n".format(t,dat[_j,1],dat[_j,2],dat[_j,3],dat[_j,4],dat[_j,5],dat[_j,6]))
    #
    #     a_planet[j,:]=dat[1:N_planet+1,1]
    #     e_planet[j,:]=dat[1:N_planet+1,2]
    #     I_planet[j,:]=dat[1:N_planet+1,3]
    #     j+=1

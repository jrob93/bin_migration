'''
Search for close encounters between binary centre of mass and the giant planets
'''
import numpy
import matplotlib.pyplot as pyplot
import os
from natsort import natsorted, ns
import subprocess
import matplotlib.gridspec as gridspec
import sys
sys.path.insert(0, '../../../grav_cloud/python_stuff/')
import py_func as pf
import rebound
import time

path='/Users/jrobinson/bin_migration/'
run='200_dump'

# path='/Users/jrobinson/bin_migration/migration'
# run='binSim'
# # run='binSim_20_0'

subprocess.Popen(["mkdir",run]) # make directory to store files for run

N_planet=4
bin_mass=0.5e18
planet_mass=numpy.array([1.8981e27,5.6831e26,8.6810e25,1.024e26])
r_hill_fac=10

run_path='{}/{}'.format(path,run)
files=next(os.walk(run_path))[2]
files = [ fi for fi in files if fi.endswith('.dump')]
files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers

primary_fname='{}/{}_primary_list.txt'.format(run,run)
single_fname='{}/{}_single_list.txt'.format(run,run)

dist=[]
times=[]
encounter_index=[]
dist=[]

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

    # data format:
    # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0

    pos=dat[:,7:10]*pf.AU
    vel=dat[:,10:13]*pf.AU/((24.0*60*60*365.25)/(2.0*numpy.pi)) #convert units!

    # search for close encounters
    i=N_planet+1
    while i < len(pos)-1:
        # Binary centre of mass
        bin_com=pf.centre_of_mass(pos[i:i+1,:],numpy.zeros(2)+bin_mass)
        # print bin_com
        # search against all planets
        for k in range(1,N_planet+1):
            # print k, pos[k,:]
            r_hill=pow(planet_mass[k-1]/(3.0*pf.M_sun),1.0/3.0)
            sep=bin_com-pos[k,:]
            _dist=numpy.linalg.norm(sep)
            if _dist/r_hill<r_hill_fac:
                print 'close encounter'
                times.append(t)
                encounter_index.append(i)
            # print _dist
            dist.append(_dist)
        i+=2

    # break

print times
print encounter_index
print numpy.amin(dist)

# # Make a plot of the distribution of close encounters
# fig = pyplot.figure()
# gs = gridspec.GridSpec(1, 1)
# ax1 = pyplot.subplot(gs[0,0])
# ax1.set_xlabel('separation/10*r_hill')
# ax1.set_ylabel('N')
# fig.suptitle('run: {}'.format(run))
#
# n_bin=100
# ax1.hist(dist/(r_hill*r_hill_fac),n_bin)
#
# save_path = '{}/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
# print save_path
# pyplot.savefig(save_path, bbox_inches='tight')
# pyplot.show()
# # pyplot.close()

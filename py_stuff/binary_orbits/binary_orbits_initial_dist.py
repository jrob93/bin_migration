'''
This script loads the initial file of a run and plots the distribution of initial binary orbits
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

subprocess.Popen(["mkdir","../{}_analysis".format(run)]) # make directory to store files for run

N_planet=4
bin_mass=0.5e18
n_bins=10

run_path='{}/{}'.format(path,run)
files=next(os.walk(run_path))[2]
files = [ fi for fi in files if fi.endswith('.dump')]
files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers

# Load the first file with the initial binaries
filename='{}/{}'.format(run_path,files[0])
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

# Calculate all binary pairs and store a, e, I
a_bin=[]
e_bin=[]
I_bin=[]

i=N_planet+1
# print numpy.linalg.norm(pos[i,:]),numpy.linalg.norm(vel[i,:]),numpy.sqrt(pf.G*pf.M_sun/numpy.linalg.norm(pos[i,:]))
primary_list=[]
single_list=[]
while i < len(pos)-1:
    if i not in single_list:
        j=i+1

        # use rebound to calculate the orbit, better handling for certain cases, e.g. circular orbits
        sim = rebound.Simulation()
        sim.G=pf.G
        sim.add(x=pos[i,0],y=pos[i,1],z=pos[i,2],
        vx=vel[i,0],vy=vel[i,1],vz=vel[i,2],
        m=bin_mass)
        sim.add(x=pos[j,0],y=pos[j,1],z=pos[j,2],
        vx=vel[j,0],vy=vel[j,1],vz=vel[j,2],
        m=bin_mass)
        orbit = sim.particles[1].calculate_orbit(sim.particles[0])

        a_bin.append(orbit.a)
        e_bin.append(orbit.e)
        I_bin.append(orbit.inc)

    i+=2

fig = pyplot.figure()
gs = gridspec.GridSpec(3, 1)
ax1 = pyplot.subplot(gs[0,0])
ax2 = pyplot.subplot(gs[1,0])
ax3 = pyplot.subplot(gs[2,0])
ax1.set_xlabel('a (m)')
ax2.set_xlabel('e')
ax3.set_xlabel('I')
ax1.set_ylabel('N')
ax2.set_ylabel('N')
ax3.set_ylabel('N')

fig.suptitle('run: {}'.format(run))

ax1.hist(a_bin,bins=n_bins)
ax2.hist(e_bin,bins=n_bins)
ax3.hist(I_bin,bins=n_bins)

pyplot.tight_layout()
save_path = '../{}_analysis/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
print save_path
pyplot.savefig(save_path, bbox_inches='tight')
pyplot.show()
# pyplot.close()

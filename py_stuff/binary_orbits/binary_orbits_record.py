'''
This script calculate the binary orbits as a function of simualtion time.
At each timestep we initialise empty lists to track:

primary_list : a list of the primary particle indices, for objects that are bound
single_list : a list of the primary particle indices, for objects that have become unbound

These are saved at the last timestep, making lists of binary objects that did, or did not,
survive the whole simulation.

We also store the binary object as a function of time, for as long as it is bound.
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

# path='/Users/jrobinson/bin_migration/'
# dirs=["prone_122","prone_123","prone_124","prone_125","prone_126",
# "prone_127","prone_128","prone_129","prone_130","prone_131","prone_132","prone_133"]
# dirs=['prone_120']

path='/Volumes/Mirkwood/bin_migration'
dirs=next(os.walk("{}/data".format(path)))[1]
dirs = [d for d in dirs if ('prone' in d) and ('test' not in d)]
# dirs=['prone_1000_timesteptest0.25','prone_1000_timesteptest0.75']
# dirs=['prone_1040_timesteptest0.25','prone_1335_timesteptest0.25']
dirs=['prone_1000_timesteptest0.25','prone_1000_timesteptest0.75',
'prone_1040_timesteptest0.25','prone_1040_timesteptest0.75',
'prone_1335_timesteptest0.25','prone_1335_timesteptest0.75']

print dirs

short=1 # if short==1 we only fidn the first and last timestep and save locally
# otherwise we find the entire binary orbit, but save to the external drive

# Variables for calculating orbits
N_planet=4 # the first N_planet lines in the data file are planet coordinates, after it is primary, secondary...
bin_mass=0.5e18 # The mass of a binary particle (m_primary = m_secondary)
# multi_search=0

# for run in dirs[::-1]: # go backwards through list
for run in dirs: # go forwards through list

    # only perform orbit calculation for certain run numbers
    # if (int(run.split('_')[-1])<1435):# or (int(run.split('_')[-1])>1380):
    #     continue

    subprocess.Popen(["mkdir","../analysis_dirs/{}_analysis".format(run)]) # make directory to store files for run

    # Load the files
    run_path='{}/data/{}'.format(path,run)
    files=next(os.walk(run_path))[2]
    files = [ fi for fi in files if fi.endswith('.dump')]
    files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers

    #---------------------------------------------------------------------------
    # Uncomment this line to find orbits for only the first and last timesteps
    if short==1:
        files=[files[0],files[-1]] # calculate only first and last
    #---------------------------------------------------------------------------

    primary_fname='../analysis_dirs/{}_analysis/{}_primary_list.txt'.format(run,run)
    single_fname='../analysis_dirs/{}_analysis/{}_single_list.txt'.format(run,run)

    # time1 = time.time()

    print 'calculate pair lists'
    for _i in range(len(files)):
        # Load data file
        filename='{}/{}'.format(run_path,files[_i])
        print filename
        with open(filename) as f:
            dat = [line.rstrip() for line in f]
        dat = [line.split() for line in dat]
        print dat[0] # header line
        t=float(dat[0][3]) # timestamp
        dat = dat[1:] # strip header line
        dat = numpy.array(dat,dtype=float) # convert data list into array

        # copy first and last data files for future analysis
        if _i==0 or _i==(len(files)-1):
            print "copy {}".format(filename)
            os.system ("cp -v {} ../analysis_dirs/{}_analysis/{}".format(filename, run, files[_i]))

        # data format:
        # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0

        pos=dat[:,7:10]*pf.AU
        vel=dat[:,10:13]*pf.AU/((24.0*60*60*365.25)/(2.0*numpy.pi)) #convert units!

        # _time1 = time.time()

        # search binary pairs
        i=N_planet+1
        primary_list=[] # initialise empty lists for this timestep to store bound
        single_list=[] # and unbound particle indices
        while i < len(pos)-1:
            if i not in single_list:
                j=i+1

                if short==1:
                    binary_record_file="../analysis_dirs/{}_analysis/binary_{:03d}_{:03d}_short.txt".format(run,i,j)
                else:
                    binary_record_file="{}/analysis/{}_analysis/binary_{:03d}_{:03d}.txt".format(path,run,i,j)

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
                # check stability
                if orbit.e>0.0 and orbit.e<1.0: # add a/R_hill check too
                    # print 'bound'
                    primary_list.append(i)

                    # print (t,
                    # orbit.a,orbit.e,orbit.inc,orbit.omega,orbit.Omega,orbit.f,
                    # dat[i,1],dat[i,2],dat[i,3],dat[i,4],dat[i,5],dat[i,6],
                    # dat[j,1],dat[j,2],dat[j,3],dat[j,4],dat[j,5],dat[j,6])

                    # Record the binary orbit
                    if _i==0: # wipe file to avoid appending errors
                        # print "wipe {}".format(binary_record_file)
                        with open(binary_record_file, 'w') as file:
                            file.write("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n".format(t,
                            orbit.a,orbit.e,orbit.inc,orbit.omega,orbit.Omega,orbit.f,
                            dat[i,1],dat[i,2],dat[i,3],dat[i,4],dat[i,5],dat[i,6],
                            dat[j,1],dat[j,2],dat[j,3],dat[j,4],dat[j,5],dat[j,6]))
                    else:
                        with open(binary_record_file, 'a') as file:
                            file.write("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n".format(t,
                            orbit.a,orbit.e,orbit.inc,orbit.omega,orbit.Omega,orbit.f,
                            dat[i,1],dat[i,2],dat[i,3],dat[i,4],dat[i,5],dat[i,6],
                            dat[j,1],dat[j,2],dat[j,3],dat[j,4],dat[j,5],dat[j,6]))

                else:
                    single_list.append(i)
            i+=2

        # _time2 = time.time()
        # print '%0.3f ms' % ((_time2-_time1)*1000.0)

    # save the lists
    numpy.savetxt(primary_fname,primary_list)
    numpy.savetxt(single_fname,single_list)

    # time2 = time.time()
    # print '%0.3f ms' % ((time2-time1)*1000.0)

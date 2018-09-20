'''
This script searches first of all calculates the binary orbits at each timestep.
If e>1 we assume the binary is unbound and the primary particle index is recorded.
For binaries that remain bound we then recalculate the orbital elements as a function of time
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
run='prone_133'

# path='/Users/jrobinson/bin_migration/migration'
# run='binSim'
# # run='binSim_20_0'

subprocess.Popen(["mkdir","../{}_analysis".format(run)]) # make directory to store files for run

N_planet=4
bin_mass=0.5e18

run_path='{}/{}'.format(path,run)
files=next(os.walk(run_path))[2]
files = [ fi for fi in files if fi.endswith('.dump')]
files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers

primary_fname='../{}_analysis/{}_primary_list.txt'.format(run,run)
single_fname='../{}_analysis/{}_single_list.txt'.format(run,run)

time1 = time.time()

# attempt to load the primary and single lists, if not calculate them
try:
    primary_list=numpy.loadtxt(primary_fname).astype(int)
    single_list=numpy.loadtxt(single_fname).astype(int)
    print 'pair lists loaded'
    # print primary_list
except:
    print 'calculate pair lists'
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

        _time1 = time.time()

        # search binary pairs
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
                # check stability
                if orbit.e>0.0 and orbit.e<1.0: # add a/R_hill check too
                    # print 'bound'
                    primary_list.append(i)
                else:
                    single_list.append(i)
            i+=2
        # print primary_list,len(primary_list)
        # print single_list,len(single_list)

        _time2 = time.time()
        print '%0.3f ms' % ((_time2-_time1)*1000.0)

        # break
    # save the lists
    numpy.savetxt(primary_fname,primary_list)
    numpy.savetxt(single_fname,single_list)

time2 = time.time()
print '%0.3f ms' % ((time2-time1)*1000.0)

fig = pyplot.figure()
gs = gridspec.GridSpec(4, 1)
ax1 = pyplot.subplot(gs[0,0])
ax2 = pyplot.subplot(gs[1,0],sharex=ax1)
ax3 = pyplot.subplot(gs[2,0],sharex=ax1)
ax4 = pyplot.subplot(gs[3,0],sharex=ax1)

ax4.set_xlabel('t (Myr)')
ax1.set_ylabel('a (km)')
# ax1.set_ylabel('a/r_hill_mutual')
ax2.set_ylabel('e')
# ax3.set_ylabel('I')
ax3.set_ylabel('I (degrees)')
ax4.set_ylabel('a_hel/q_hel (primary)')
fig.suptitle('run: {}'.format(run))

# pick a binary that survived and re calculate orbits to plot
# for i in primary_list:
i=primary_list[0]
j=i+1

a_bin=[]
e_bin=[]
I_bin=[]
time=[]
r_hill=[]

a_hel=[]
e_hel=[]
a_hel_N=[]
e_hel_N=[]

# for k in range(len(files)):
for k in range(0,len(files),1000):

    filename='{}/{}'.format(run_path,files[k])
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
    print orbit
    # # check stability
    # if orbit.e>0.0 and orbit.e<1.0: # add a/R_hill check too
    #     # print 'bound'
    #     primary_list.append(i)
    # else:
    #     single_list.append(i)

    a_bin.append(orbit.a)
    e_bin.append(orbit.e)
    I_bin.append(orbit.inc)
    time.append(t)

    # record the helio orbit of the primary
    a_hel.append(dat[i,1]*pf.AU)
    e_hel.append(dat[i,2])

    # record the helio orbit of Neptune
    a_hel_N.append(dat[N_planet,1]*pf.AU)
    e_hel_N .append(dat[N_planet,2])

    # # add mutual hill radius scale
    # # r_hill=pow(planet_mass[k-1]/(3.0*pf.M_sun),1.0/3.0) ?
    # mi=sim.particles[0].m
    # mj=sim.particles[1].m
    # ai=dat[i,1]*pf.AU
    # aj=dat[j,1]*pf.AU
    # r_hill_mutual=pow(mi+mj/(3.0*pf.M_sun),1.0/3.0)*(((ai*mi)+(aj*mj))/(mi+mj))
    # r_hill.append(r_hill_mutual)

I_bin=numpy.array(I_bin)/(numpy.pi/180.0)

ax1.plot(time,numpy.array(a_bin)/1e3)
# ax1.plot(time,numpy.array(a_bin)/numpy.array(r_hill))
ax2.plot(time,e_bin)
ax3.plot(time,I_bin)

# # plot the helio centric a/q for binary, neptune and the main neptune resonances
# q_hel=numpy.array(a_hel)*(1.0-numpy.array(e_hel))
# ax4.plot(time,numpy.array(a_hel)/q_hel)
# q_hel_N=numpy.array(a_hel_N)*(1.0-numpy.array(e_hel_N))
# ax4.plot(time,numpy.array(a_hel_N)/q_hel_N,label="1:1",c='k',linestyle=':',alpha=0.5)
# ax4.plot(time,(4.0/3.0)*numpy.array(a_hel_N)/q_hel_N,label="4:3",c='k',linestyle=':',alpha=0.5)
# ax4.plot(time,(3.0/2.0)*numpy.array(a_hel_N)/q_hel_N,label="3:2",c='k',linestyle='-',alpha=0.5)
# ax4.plot(time,2.0*numpy.array(a_hel_N)/q_hel_N,label="2:1",c='k',linestyle='--',alpha=0.5)

# plot only the heliocentric a for binary, neptune and resonances
ax4.plot(time,numpy.array(a_hel)/pf.AU)
ax4.plot(time,numpy.array(a_hel_N)/pf.AU,label="1:1",c='k',linestyle=':',alpha=0.5)
ax4.plot(time,(4.0/3.0)*numpy.array(a_hel_N)/pf.AU,label="4:3",c='k',linestyle=':',alpha=0.5)
ax4.plot(time,(3.0/2.0)*numpy.array(a_hel_N)/pf.AU,label="3:2",c='k',linestyle='--',alpha=0.5)
ax4.plot(time,2.0*numpy.array(a_hel_N)/pf.AU,label="2:1",c='k',linestyle='-',alpha=0.5)
ax4.set_ylabel('a_hel(AU)')

ax4.legend()

save_path = '../{}_analysis/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
print save_path
pyplot.savefig(save_path, bbox_inches='tight')
pyplot.show()
# pyplot.close()

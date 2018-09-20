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

# path='/Users/jrobinson/bin_migration/'
# run='prone_121'

path='/Volumes/Mirkwood/bin_migration/'
dirs=next(os.walk(path))[1]
dirs = [d for d in dirs if 'prone' in d]

N_planet=4
bin_mass=0.5e18
no_survive=0 # If there are no survivors, equals 1

for run in dirs:

    if int(run.split('_')[1])<1049:
        continue

    print run

    fig = pyplot.figure()
    fig.set_size_inches(15, 10)
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
    # ax4.set_ylabel('a_hel/q_hel (primary)')

    primary_fname='../analysis_dirs/{}_analysis/{}_primary_list.txt'.format(run,run)
    single_fname='../analysis_dirs/{}_analysis/{}_single_list.txt'.format(run,run)

    primary_list=numpy.loadtxt(primary_fname).astype(int)
    single_list=numpy.loadtxt(single_fname).astype(int)
    print 'pair lists loaded'
    if len(primary_list)==0:
        no_survive=1
        print "no survivors!"
        continue

    i=primary_list[0]
    j=i+1

    print i,j
    fig.suptitle('run: {}, binary pair: {}, {}'.format(run,i,j))

    # Load the binary data
    binary_record_file="../analysis_dirs/{}_analysis/binary_{:03d}_{:03d}.txt".format(run,i,j)
    filename='{}'.format(binary_record_file)
    print filename
    with open(filename) as f:
        dat = [line.rstrip() for line in f]
    dat = [line.split() for line in dat]
    dat = numpy.array(dat,dtype=float) # convert data list into array
    # print numpy.shape(dat)
    # print dat
    # data format:
    # t,a,e,inc,omega,Omega,f
    time=dat[:,0]
    a_bin=dat[:,1]
    e_bin=dat[:,2]
    I_bin=dat[:,3]
    # print dat[:,0][-1],dat[:,8][-1]
    a_hel=dat[:,7]*pf.AU
    e_hel=dat[:,8]

    # Load the planet heliocentric data
    planet_record_file="../analysis_dirs/{}_analysis/planet_{:03d}.txt".format(run,N_planet)
    filename='{}'.format(planet_record_file)
    print filename
    with open(filename) as f:
        dat_p = [line.rstrip() for line in f]
    dat_p = [line.split() for line in dat_p]
    dat_p = numpy.array(dat_p,dtype=float) # convert data list into array
    a_hel_N=(dat_p[:,1]*pf.AU)
    e_hel_N=(dat_p[:,2])

    I_bin=numpy.array(I_bin)/(numpy.pi/180.0)

    ax1.plot(time,numpy.array(a_bin)/1e3)
    ax2.plot(time,e_bin)
    ax3.plot(time,I_bin)

    print a_hel[-1],a_hel_N[-1]
    print a_hel[-1]/pf.AU,a_hel_N[-1]/pf.AU

    # plot only the heliocentric a for binary, neptune and resonances
    ax4.plot(time,numpy.array(a_hel)/pf.AU)

    res_11=[1,1]
    res_43=[4,3]
    res_32=[3,2]
    res_21=[2,1]
    ress=[res_11,res_43,res_32,res_21]
    res_ls=[':','--','-.','-']
    for _i in range(len(ress)):
        _j=ress[_i][0]
        _k=ress[_i][1]
        ax4.plot(time,((numpy.array(a_hel_N))*((float(j)/_k)**(2.0/3.0)))/pf.AU,label="{}:{}".format(_j,_k),c='k',linestyle=res_ls[_i],alpha=0.5)

    # ax4.plot(time,(4.0/3.0)*numpy.array(a_hel_N)/pf.AU,label="4:3",c='k',linestyle='--',alpha=0.5)
    # ax4.plot(time,(3.0/2.0)*numpy.array(a_hel_N)/pf.AU,label="3:2",c='k',linestyle='-.',alpha=0.5)
    # ax4.plot(time,2.0*numpy.array(a_hel_N)/pf.AU,label="2:1",c='k',linestyle='-',alpha=0.5)
    ax4.set_ylabel('a_hel (AU)')

    ax4.legend()

    # add mutual hill radius scale
    mi=bin_mass
    mj=mi
    ai=a_hel
    aj=dat[:,13]*pf.AU
    r_hill_mutual=pow((mi+mj)/(3.0*pf.M_sun),1.0/3.0)*(((ai*mi)+(aj*mj))/(mi+mj))

    print a_bin
    print "r_hill:",r_hill_mutual
    print numpy.array(a_bin)/r_hill_mutual
    # # Add additional y scale to ax1
    # km_to_hill_radii=1e3*
    # ax2 = ax1.twinx()
    # mn, mx = ax2.get_ylim()
    # ax2.set_ylim(mn*1000, mx*1000)
    # ax2.set_ylabel('Sv')
    _ax1 = ax1.twinx()
    colours=pyplot.rcParams['axes.prop_cycle'].by_key()['color']
    # _ax1.plot(time,numpy.array(a_bin)/r_hill_mutual, color=colours[1])
    _ax1.plot(time,numpy.array(a_bin)/r_hill_mutual, color='k',alpha=0.5)
    _ax1.set_ylabel('a_bin/r_hill_mutual')

    save_path = '../analysis_dirs/{}_analysis/{}_{}_{:03d}_{:03d}.png'.format(run,os.path.basename(__file__).split('.')[0],run,i,j)
    print save_path
    pyplot.savefig(save_path, bbox_inches='tight')
    # pyplot.show()
    pyplot.close()

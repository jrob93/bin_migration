import numpy
import matplotlib.pyplot as pyplot
import os
from natsort import natsorted, ns
import subprocess
import matplotlib.gridspec as gridspec
# import sys
# sys.path.insert(0, '../../../grav_cloud/python_stuff/')
import py_func as pf

# path='/Users/jrobinson/bin_migration/migration'
# run='binSim'
path='/Volumes/Mirkwood/bin_migration/data'
run='prone_120'
local_path='../{}_analysis'.format(run)

plot_option=2
use_local_files=1
lim=100.0
n=100
N_planet=4
colours=pf.pyplot_colours

# subprocess.Popen(["mkdir",run]) # make directory to store images

if use_local_files==1:
    run_path=local_path

else:
    run_path='{}/{}'.format(path,run)


files=next(os.walk(run_path))[2]
files = [ fi for fi in files if fi.endswith('.dump')]
files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers


last_file=files[-1]
n_step=10000
files=files[::n_step]
files.append(last_file)
print files
# exit()

# Make a separate plot for a range of timesteps
if plot_option==0:
    for i in range(len(files)):
        filename='{}/{}'.format(run_path,files[i])
        print filename
        # continue
        with open(filename) as f:
            dat = [line.rstrip() for line in f]
        dat = [line.split() for line in dat]
        print dat[0] # header line
        t=float(dat[0][3]) # timestamp
        dat = dat[1:] # strip header line
        dat = numpy.array(dat,dtype=float) # convert data list into array

        # data format:
        # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0

        pos_xy=dat[:,7:9]

        fig = pyplot.figure()
        gs = gridspec.GridSpec(1, 1)
        ax1 = pyplot.subplot(gs[0,0])
        ax1.set_aspect("equal")
        ax1.set_xlabel('x(AU)')
        ax1.set_ylabel('y(AU)')
        ax1.set_xlim(-lim,lim)
        ax1.set_ylim(-lim,lim)
        fig.suptitle('run: {}, time: {:.2e} Myr'.format(run,t))
        ax1.scatter(pos_xy[0,0],pos_xy[0,1],color='k',marker='x') # plot sun

        # plot planets and add planet orbits
        for j in range(1,N_planet+1):
            # ax1.scatter(pos_xy[1:N_planet+1,0],pos_xy[1:N_planet+1,1],color=colours[j])
            ax1.scatter(pos_xy[j,0],pos_xy[j,1],color=colours[j])
            orb=dat[j,1:7]
            # print orb
            orb_pos=pf.planet_orbit(orb,n)
            ax1.plot(orb_pos[:,0],orb_pos[:,1],color=colours[j])
        ax1.scatter(pos_xy[N_planet+2:,0],pos_xy[N_planet+2:,1]) # plot small bodies

        save_path = '../{}_analysis/{}.png'.format(run,files[i])
        print save_path
        pyplot.savefig(save_path, bbox_inches='tight')
        # pyplot.show()
        pyplot.close()

        #break

# Make a single figure showing the first and last timestep
if plot_option==1:

    fig = pyplot.figure()
    gs = gridspec.GridSpec(1, 2)
    ax1 = pyplot.subplot(gs[0,0])
    ax1.set_aspect("equal")
    ax1.set_xlabel('x(AU)')
    ax1.set_ylabel('y(AU)')
    ax1.set_xlim(-lim,lim)
    ax1.set_ylim(-lim,lim)

    ax2 = pyplot.subplot(gs[0,1])
    ax2.set_aspect("equal")
    ax2.set_xlabel('x(AU)')
    ax2.set_ylabel('y(AU)')
    ax2.set_xlim(-lim,lim)
    ax2.set_ylim(-lim,lim)

    fig.suptitle('run: {}'.format(run))

    a=[ax1,ax2]

    files=[files[0],files[-1]]

    for i in range(len(files)):
        filename='{}/{}'.format(run_path,files[i])
        print filename
        # continue
        with open(filename) as f:
            dat = [line.rstrip() for line in f]
        dat = [line.split() for line in dat]
        print dat[0] # header line
        t=float(dat[0][3]) # timestamp
        dat = dat[1:] # strip header line
        dat = numpy.array(dat,dtype=float) # convert data list into array

        # data format:
        # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0

        pos_xy=dat[:,7:9]

        a[i].set_title('time: {:.2e} Myr'.format(t))
        a[i].scatter(pos_xy[0,0],pos_xy[0,1],color='k',marker='x') # plot sun

        # plot planets and add planet orbits
        for j in range(1,N_planet+1):
            # ax1.scatter(pos_xy[1:N_planet+1,0],pos_xy[1:N_planet+1,1],color=colours[j])
            a[i].scatter(pos_xy[j,0],pos_xy[j,1],color=colours[j])
            orb=dat[j,1:7]
            # print orb
            orb_pos=pf.planet_orbit(orb,n)
            a[i].plot(orb_pos[:,0],orb_pos[:,1],color=colours[j])
        # a[i].scatter(pos_xy[N_planet+2:,0],pos_xy[N_planet+2:,1]) # plot small bodies
        a[i].scatter(pos_xy[N_planet+1:,0],pos_xy[N_planet+1:,1]) # plot small bodies

    save_path = '../{}_analysis/{}_summary.png'.format(run,os.path.basename(__file__).split('.')[0])
    print save_path
    # pyplot.savefig(save_path, bbox_inches='tight')
    pyplot.show()
    # pyplot.close()

# Make a single figure showing the first and last timestep, add binary and single planetesimals
if plot_option==2:

    fig = pyplot.figure()
    gs = gridspec.GridSpec(1, 2)
    ax1 = pyplot.subplot(gs[0,0])
    ax1.set_aspect("equal")
    ax1.set_xlabel('x(AU)')
    ax1.set_ylabel('y(AU)')
    ax1.set_xlim(-lim,lim)
    ax1.set_ylim(-lim,lim)

    ax2 = pyplot.subplot(gs[0,1])
    ax2.set_aspect("equal")
    ax2.set_xlabel('x(AU)')
    # ax2.set_ylabel('y(AU)')
    ax2.set_xlim(-lim,lim)
    ax2.set_ylim(-lim,lim)
    ax2.set_yticklabels([])

    fig.suptitle('run: {}'.format(run))

    a=[ax1,ax2]

    files=[files[0],files[-1]]

    for i in range(len(files)):
        filename='{}/{}'.format(run_path,files[i])
        # filename='{}/{}'.format(local_path,files[i])

        print filename
        # continue
        with open(filename) as f:
            dat = [line.rstrip() for line in f]
        dat = [line.split() for line in dat]
        print dat[0] # header line
        t=float(dat[0][3]) # timestamp
        dat = dat[1:] # strip header line
        dat = numpy.array(dat,dtype=float) # convert data list into array

        # data format:
        # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0

        pos_xy=dat[:,7:9]

        a[i].set_title('time: {:.2e} Myr'.format(t))
        a[i].scatter(pos_xy[0,0],pos_xy[0,1],color='k',marker='x') # plot sun

        # plot planets and add planet orbits
        for j in range(1,N_planet+1):
            # ax1.scatter(pos_xy[1:N_planet+1,0],pos_xy[1:N_planet+1,1],color=colours[j])
            a[i].scatter(pos_xy[j,0],pos_xy[j,1],color=colours[j+2],zorder=2)
            orb=dat[j,1:7]
            # print orb
            orb_pos=pf.planet_orbit(orb,n)
            a[i].plot(orb_pos[:,0],orb_pos[:,1],color=colours[j+2],zorder=2)

        secondary_size=15

        if i==0: # all bodies are binary
            a[i].scatter(pos_xy[N_planet+1::2,0],pos_xy[N_planet+1::2,1],color=colours[0],zorder=0)
            a[i].scatter(pos_xy[N_planet+1+1::2,0],pos_xy[N_planet+1+1::2,1],color=colours[0],zorder=1,s=secondary_size,edgecolors=['k']*len(pos_xy[N_planet+1+1::2,0]))

        else: # load lists to plot binaries and singles separately
            primary_fname='../{}_analysis/{}_primary_list.txt'.format(run,run)
            # single_fname='../{}_analysis/{}_single_list.txt'.format(run,run)
            primary_list=numpy.loadtxt(primary_fname).astype(int)
            # single_list=numpy.loadtxt(single_fname).astype(int)
            print len(dat)
            for k in range(N_planet+1,len(dat),2):
                if k in primary_list: # plot binaries
                    a[i].scatter(pos_xy[k,0],pos_xy[k,1],color=colours[1],zorder=0)
                    a[i].scatter(pos_xy[k+1,0],pos_xy[k+1,1],color=colours[1],zorder=1,s=secondary_size,edgecolors=['k'])
                else: #plot singles
                    a[i].scatter(pos_xy[k,0],pos_xy[k,1],color='k',alpha=0.5,zorder=0,marker='^',edgecolors='none')
                    a[i].scatter(pos_xy[k+1,0],pos_xy[k+1,1],color='k',alpha=0.5,zorder=0,marker='^',edgecolors='none')

    save_path = '../{}_analysis/{}_summary.png'.format(run,os.path.basename(__file__).split('.')[0])
    print save_path
    pyplot.savefig(save_path, bbox_inches='tight')
    pyplot.show()
    # pyplot.close()

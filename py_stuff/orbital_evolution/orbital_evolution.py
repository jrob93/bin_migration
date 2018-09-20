'''
TO DO: ADD A HISTOGRAM THAT SHOWS FINAL BINARY INCLINATIONS
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
import seaborn as sns
import scipy
import matplotlib.patches as patches
from matplotlib.patches import ConnectionPatch
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

# path='/Users/jrobinson/bin_migration/'
# # run='200_dump'
# run='prone_121'

# path='/Users/jrobinson/bin_migration/migration'
# run='binSim'
# # run='binSim_20_0'

path='/Users/jrobinson/bin_migration/py_stuff'
dirs=next(os.walk(path))[1]
dirs = [d for d in dirs if 'prone' in d]

dirs=['prone_120_analysis']

# subprocess.Popen(["mkdir","../{}_analysis".format(run)]) # make directory to store files for run

N_planet=4
bin_mass=0.5e18
n_bins=10 # bins for inclination histogram
wide_cut=0.05 # a/r_hill, the normalised binary separation above which we consider the binary to be 'wide'

for run in dirs:

    run_path='{}/{}'.format(path,run)
    files=next(os.walk(run_path))[2]
    files = [ fi for fi in files if fi.endswith('.dump')]
    files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers

    # calculate only the first and last state
    files=[files[0],files[-1]]

    # primary_fname='../{}_analysis/{}_primary_list.txt'.format(run,run)
    # single_fname='../{}_analysis/{}_single_list.txt'.format(run,run)

    run_name="_".join(run.split("_")[:2])
    primary_fname='../{}/{}_primary_list.txt'.format(run,run_name)
    single_fname='../{}/{}_single_list.txt'.format(run,run_name)

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

        # save the lists
        numpy.savetxt(primary_fname,primary_list)
        numpy.savetxt(single_fname,single_list)

    # now plot the orbital elements of the binaries at first and last timestep

    fig = pyplot.figure()
    fig.set_size_inches(15, 10)

    # gs = gridspec.GridSpec(2, 2, width_ratios=[0.75,0.25])
    # ax1 = pyplot.subplot(gs[0,0])
    # ax2 = pyplot.subplot(gs[1,0],sharex=ax1)
    # ax3 = pyplot.subplot(gs[1,1],sharey=ax2)
    gs = gridspec.GridSpec(2, 1)
    ax1 = pyplot.subplot(gs[0,0])
    ax2 = pyplot.subplot(gs[1,0],sharex=ax1)

    ax1.set_ylabel('e')
    ax1.set_xlabel('a (m)')
    # ax2.set_ylabel('I (rad)')
    ax2.set_ylabel('I (degrees)')
    ax2.set_xlabel('a (m)')
    # ax3.set_xlabel('I (degrees)')
    # ax3.set_ylabel('N')

    fig.suptitle('run: {}'.format(run))

    files=[files[0],files[-1]]

    for k in range(len(files)):

        a_bin=[]
        e_bin=[]
        I_bin=[]

        if k==0:
            a_bin_survivor=[]
            e_bin_survivor=[]
            I_bin_survivor=[]
            a_bin_survivor=[]
            e_bin_survivor=[]
            I_bin_survivor=[]

        time=[]

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

        # search binary pairs
        i=N_planet+1

        while i < len(pos)-1:
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
                # primary_list.append(i)
                a_bin.append(orbit.a)
                e_bin.append(orbit.e)
                I_bin.append(orbit.inc)
                time.append(t)

            # record the orbits that go on to survive
            if k==0:
                if i in primary_list:
                    a_bin_survivor.append(orbit.a)
                    e_bin_survivor.append(orbit.e)
                    I_bin_survivor.append(orbit.inc)

            i+=2

        I_bin=numpy.array(I_bin)/(numpy.pi/180.0) # convert to degrees

        ax1.scatter(a_bin,e_bin,label="t={}Myr".format(t),zorder=2)
        ax2.scatter(a_bin,I_bin,label="N_binaries={}".format(len(a_bin)),zorder=2)

        colours=pyplot.rcParams['axes.prop_cycle'].by_key()['color']

        # _ax2 = ax2.twiny()
        # if k==0:
        #     # Add inclination histogram for initial conditions
        #     center, hist, width=pf.hist_dist(I_bin,n_bins)
        #     _ax2.invert_xaxis()
        #     _ax2.barh(y=center, height=width, width=hist, align='center', color=colours[0], edgecolor='k',zorder=0,alpha=0.5,label='t={}Myr'.format(t))
        #     _ax2.tick_params('x', colors=colours[0])
        # else:
        #     # Add inclination histogram for final states
        #     center, hist, width=pf.hist_dist(I_bin,n_bins)
        #     _ax2.invert_xaxis()
        #     _ax2.barh(y=center, height=width, width=hist, align='center', color=colours[1], edgecolor='k',zorder=0,alpha=0.5,label='t={}Myr'.format(t))
        #     _ax2.tick_params('x', colors=colours[1])

        if k==0:

            n_order=100
            hist_scale_e = numpy.amax(a_bin)/n_order# 0.25e6
            hist_scale_I = numpy.amax(a_bin)/n_order

            # Add eccentricity histogram for initial conditions
            center, hist, width=pf.hist_dist(e_bin,n_bins)
            hist=hist*hist_scale_e
            ax1.barh(y=center, height=width, width=hist, align='center', color=colours[0], edgecolor='k',zorder=0,alpha=0.5)#,label='t={}Myr'.format(t))

            # add an inset for initial eccentricity
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            axins = inset_axes(ax1,width="20%",  height="30%",loc='upper left', bbox_to_anchor=(0.05, 0.25, 0.5, 0.75), bbox_transform=ax1.transAxes) # bbox_to_anchor xpos, _ , xlen
            axins.barh(y=center, height=width, width=hist, align='center', color=colours[0], edgecolor='k',zorder=0,alpha=0.5)#,label='t={}Myr'.format(t))

            # ip = InsetPosition(axins, [1.0, 1.0, 1.0, 1.0]) #posx, posy, width, height
            # axins.set_axes_locator(ip)
            #
            x1, x2, y1, y2 = 0, 2.5e6, numpy.amin(e_bin), numpy.amax(e_bin) # specify the limits
            axins.set_xlim(x1, x2) # apply the x-limits
            axins.set_ylim(y1, y2) # apply the y-limits
            axins.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
            pyplot.xticks(visible=False)

            #add Rectangle and lines from ax3 to ax6
            ax1.add_patch(patches.Rectangle((x1, y1),x2-x1,y2-y1,fill=False,color='r'))
            xy1 = (x1,y1)
            xy2 = (x1,y2)
            xy3 = (x2,y1)
            xy4 = (x2,y2)
            con = ConnectionPatch(xyA=xy3, xyB=xy3, coordsA="data", coordsB="data",axesA=axins, axesB=ax1, color="red",zorder=0)
            axins.add_artist(con)
            con = ConnectionPatch(xyA=xy1, xyB=xy1, coordsA="data", coordsB="data",axesA=axins, axesB=ax1, color="red",zorder=0)
            axins.add_artist(con)

            # from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
            # axins = zoomed_inset_axes(ax1, 2.5, loc=2) # zoom-factor: 2.5, location: upper-left
            # axins.barh(y=center, height=width, width=hist, align='center', color=colours[0], edgecolor='k',zorder=0,alpha=0.5)#,label='t={}Myr'.format(t))
            # x1, x2, y1, y2 = 0, 2.5e6, 0.0002, 0.00023 # specify the limits
            # axins.set_xlim(x1, x2) # apply the x-limits
            # axins.set_ylim(y1, y2) # apply the y-limits
            # pyplot.yticks(visible=False)
            # pyplot.xticks(visible=False)
            # from mpl_toolkits.axes_grid1.inset_locator import mark_inset
            # mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")

            # Add inclination histogram for initial conditions
            center, hist, width=pf.hist_dist(I_bin,n_bins)
            hist=hist*hist_scale_I
            ax2.barh(y=center, height=width, width=hist, align='center', color=colours[0], edgecolor='k',zorder=0,alpha=0.5)#,label='t={}Myr'.format(t))

        else:
            # Add eccentricity histogram for final states
            center, hist, width=pf.hist_dist(e_bin,n_bins)
            hist=hist*hist_scale_e
            # plot bars on the left
            ax1.barh(y=center, height=width, width=hist, align='center', color=colours[1], edgecolor='k',zorder=0,alpha=0.5)#,label='t={}Myr'.format(t))
            # plot bars on the right
            # ax1.barh(y=center, height=width, width=-hist, align='center', color=colours[1], edgecolor='k',zorder=0,alpha=0.5,left=numpy.amax(ax1.get_xlim()))

            # Add inclination histogram for final states
            center, hist, width=pf.hist_dist(I_bin,n_bins)
            hist=hist*hist_scale_I
            # plot bars on the left
            ax2.barh(y=center, height=width, width=hist, align='center', color=colours[1], edgecolor='k',zorder=0,alpha=0.5)#,label='t={}Myr'.format(t))
            # plot bars on the right]
            # ax2.barh(y=center, height=width, width=-hist, align='center', color=colours[1], edgecolor='k',zorder=0,alpha=0.5,left=numpy.amax(ax2.get_xlim()))

    # ax3.set_ylabel('I (degrees)')
    # ax3.set_xlabel('density')

    I_bin_survivor=numpy.array(I_bin_survivor)/(numpy.pi/180.0) # convert to degrees

    print len(a_bin),len(a_bin_survivor)
    # plot the initial conditions of the binaries that survive
    ax1.scatter(a_bin_survivor,e_bin_survivor,edgecolor='k',facecolors='none',zorder=2)
    ax2.scatter(a_bin_survivor,I_bin_survivor,edgecolor='k',facecolors='none',zorder=2)

    # plot a line to show evolution from initial to final positions
    if k==1:
        for j in range(len(a_bin)):
            ax1.plot([a_bin_survivor[j],a_bin[j]],[e_bin_survivor[j],e_bin[j]],color='k',alpha=0.2,zorder=2)
            ax2.plot([a_bin_survivor[j],a_bin[j]],[I_bin_survivor[j],I_bin[j]],color='k',alpha=0.2,zorder=2)

    ax1.legend()
    ax2.legend()
    # ax3.legend()

    # pyplot.tight_layout()

    # save_path = '../{}_analysis/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
    save_path = '../{}/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)

    print save_path
    pyplot.savefig(save_path, bbox_inches='tight')
    pyplot.show()
    # pyplot.close()
    break

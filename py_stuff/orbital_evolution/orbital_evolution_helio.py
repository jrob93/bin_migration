'''
This script plots the a, e, I heliocentric orbital elements of binaries in the simulation.
Initial and final heliocentric orbits are plotted for the binaires.
The initial conditions of binaries that survive are highlighted,
and they are connected to their final positions by a line.
Binaries that evolve to be within the cold classical region are also highlighted.
In addition we mark the heliocentric orbits of single bodies at the end of the simualtion.
The width of the binaries is indicated by the markersize.

In all cases we only consider objects on bound heliocentric orbits:
a_hel>0.0
e_hel<1.0
a_hel<a_hel_lim
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
import pandas as pd
import pickle
import io
import matplotlib.patches as patches
from matplotlib.patches import ConnectionPatch

# path='/Users/jrobinson/bin_migration/'
# # run='200_dump'
# run='prone_121'

path='/Users/jrobinson/bin_migration/py_stuff/analysis_dirs'
# run='prone_1239'
dirs=next(os.walk(path))[1]
dirs = [d for d in dirs if 'prone' in d]

# dirs=['prone_1002_analysis']
dirs=['prone_1040_timesteptest0.75_analysis']

# subprocess.Popen(["mkdir","../{}_analysis".format(run)]) # make directory to store files for run

#-------------------------------------------------------------------------------
# set up initial variables
N_planet=4 # number of planets, used to skip straight to binary data
bin_mass=0.5e18 # kg, mass of a binary component
mi=bin_mass
mj=mi
a_hel_lim=200*pf.AU # m, maximum semimajor axis of interest, beyond this we cut any surviving binaries
wide_cut=0.05 # a/r_hill, the normalised binary separation above which we consider the binary to be 'wide'
#-------------------------------------------------------------------------------
# plotting and analysis variables
CC_size=75 # marker size for crosses to denote the point is in the Cold Classical region
survivor_lw=2 # linewidth to highlight survivors
binary_load=0 # do we load binary data from files?
s1=25 # marker size for tight binaries
s2=90 # marker size for wide binaries
no_survive=0 # If there are no survivors, equals 1

for run in dirs:

    # if (int(run.split('_')[1])>1000) or (int(run.split('_')[1])<130):
    #     continue

    #-------------------------------------------------------------------------------
    # define the path to data, create file list, and define analysis files to try load
    # run_path='{}/{}_analysis'.format(path,run)
    run_path='{}/{}'.format(path,run)
    files=next(os.walk(run_path))[2]
    files = [ fi for fi in files if fi.endswith('.dump')]
    files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers
    # primary_fname='../{}_analysis/{}_primary_list.txt'.format(run,run)
    # single_fname='../{}_analysis/{}_single_list.txt'.format(run,run)
    # fname_df_orb_survivor='../{}_analysis/{}_df_orb_survivor.txt'.format(run,run)
    # fname_df_orb_initial='../{}_analysis/{}_df_orb_initial.txt'.format(run,run)
    # fname_df_orb_final='../{}_analysis/{}_df_orb_final.txt'.format(run,run)

    run_name="_".join(run.split("_")[:2])
    primary_fname='../{}/{}_primary_list.txt'.format(run,run_name)
    single_fname='../{}/{}_single_list.txt'.format(run,run_name)
    fname_df_orb_survivor='../{}/{}_df_orb_survivor.txt'.format(run,run_name)
    fname_df_orb_initial='../{}/{}_df_orb_initial.txt'.format(run,run_name)
    fname_df_orb_final='../{}/{}_df_orb_final.txt'.format(run,run_name)

    #-------------------------------------------------------------------------------
    # attempt to load binary data frames
    try:
        df_orb_initial=pd.read_csv(fname_df_orb_initial,sep="\t",index_col=0) #load the dataframe, ignoring the index column
        df_orb_final=pd.read_csv(fname_df_orb_final,sep="\t",index_col=0)
        df_orb_survivor=pd.read_csv(fname_df_orb_survivor,sep="\t",index_col=0)
        print "binary dataframes loaded"
    except:
        binary_load=1
        print "will load individual binary orbits"
    #-------------------------------------------------------------------------------
    # attempt to load the primary and single lists, if not calculate them
    try:
        primary_list=numpy.loadtxt(primary_fname).astype(int)
        single_list=numpy.loadtxt(single_fname).astype(int)
        print 'pair lists loaded'
        # print primary_list
        if len(primary_list)==0:
            no_survive=1
            print "no survivors!"
    except:
        print 'calculate pair lists'

        # calculate only the first and last state
        files=[files[0],files[-1]]

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

            # _time1 = time.time()

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
    #-------------------------------------------------------------------------------
    # now plot the heliocentric orbital elements of the primary binary particles at first and last timestep

    fig = pyplot.figure()
    fig.set_size_inches(15, 10)
    gs = gridspec.GridSpec(2, 2, width_ratios=[0.75,0.25])
    ax1 = pyplot.subplot(gs[0,0])
    ax2 = pyplot.subplot(gs[1,0],sharex=ax1)
    ax3 = pyplot.subplot(gs[0,1])
    ax4 = pyplot.subplot(gs[1,1],sharex=ax3)

    ax1.set_ylabel('e')
    ax1.set_xlabel('a (AU)')
    ax2.set_ylabel('I (degress)')
    ax2.set_xlabel('a (AU)')
    ax3.set_ylabel('e')
    ax3.set_xlabel('a (AU)')
    ax4.set_ylabel('I (degress)')
    ax4.set_xlabel('a (AU)')
    fig.suptitle('run: {}'.format(run))

    # only consider the initial and final file
    files=[files[0],files[-1]]

    # record the number of objects ejected
    ejection_binary=0
    ejection_single=0

    colours=pyplot.rcParams['axes.prop_cycle'].by_key()['color']

    for k in range(len(files)):

        a_hel=[]
        e_hel=[]
        I_hel=[]
        r_hill=[]
        a_bin=[]

        if k==0:
            a_hel_survivor=[]
            e_hel_survivor=[]
            I_hel_survivor=[]
            r_hill_survivor=[]
            a_bin_survivor=[]
            a_hel_single=[]
            e_hel_single=[]
            I_hel_single=[]

        time=[]

        # Load the .dump files to obtain the heliocentric orbital elements
        # Not that we assume binary helio elements = primary helio elements
        filename='{}/{}'.format(run_path,files[k])
        print filename
        with open(filename) as f:
            dat = [line.rstrip() for line in f]
        dat = [line.split() for line in dat]
        # print dat[0] # header line
        t=float(dat[0][3]) # timestamp
        dat = dat[1:] # strip header line
        dat = numpy.array(dat,dtype=float) # convert data list into array

        # data format:
        # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0

        a=dat[:,1]*pf.AU # m
        e=dat[:,2]
        I=dat[:,3] # radians
        pos=dat[:,7:10]*pf.AU # m
        vel=dat[:,10:13]*pf.AU/((24.0*60*60*365.25)/(2.0*numpy.pi)) # m/s

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
                if a[i]>0.0 and a[i]<a_hel_lim and e[i]<1.0: # consider if bound heliocentrically
                    a_hel.append(a[i])
                    e_hel.append(e[i])
                    I_hel.append(I[i])
                    time.append(t)
                    r_hill.append(pow((mi+mj)/(3.0*pf.M_sun),1.0/3.0)*(((a[i]*mi)+(a[j]*mj))/(mi+mj)))
                else: # binary is heliocentrically unbound
                    ejection_binary+=2

            # record the single objects (which are still bound to the solar system)
            else:
                if a[i]>0.0 and a[i]<a_hel_lim and e[i]<1.0:
                    a_hel_single.append(a[i])
                    e_hel_single.append(e[i])
                    I_hel_single.append(I[i])
                else:
                    ejection_single+=1
                if a[j]>0.0 and a[j]<a_hel_lim and e[j]<1.0:
                    a_hel_single.append(a[j])
                    e_hel_single.append(e[j])
                    I_hel_single.append(I[j])
                else:
                    ejection_single+=1

            # record the initial orbits that go on to survive
            if k==0: # at the first timestep
                if i in primary_list: # keep only objects that survive
                    if a[i]>0.0 and e[i]<1.0: # record only heliocentrically bound binaries
                        a_hel_survivor.append(a[i])
                        e_hel_survivor.append(e[i])
                        I_hel_survivor.append(I[i])
                        r_hill_survivor.append(pow((mi+mj)/(3.0*pf.M_sun),1.0/3.0)*(((a[i]*mi)+(a[j]*mj))/(mi+mj)))

            i+=2

        if binary_load==1:
            # Load the binary data in order to calculate binary separation a_bin/r_hill
            if k==0:
                _i=N_planet+1
                while _i < len(pos)-1:
                    _j=_i+1
                    # binary_record_file="../{}_analysis/binary_{:03d}_{:03d}.txt".format(run,_i,_j)
                    binary_record_file="../{}/binary_{:03d}_{:03d}.txt".format(run,_i,_j)
                    filename='{}'.format(binary_record_file)
                    print filename
                    with open(filename) as f:
                        dat = [line.rstrip() for line in f]
                    dat = [line.split() for line in dat]
                    dat = numpy.array(dat,dtype=float) # convert data list into array
                    # data format:
                    # t,a,e,inc,omega,Omega,f
                    a_bin.append(dat[0,1])

                    if no_survive==1:
                        continue

                    if _i in primary_list: # keep only objects that survive
                        if a[_i]>0.0 and e[_i]<1.0: # record only heliocentrically bound binaries
                            a_bin_survivor.append(dat[0,1])

                    _i+=2
            if k!=0:

                if no_survive==1:
                    continue

                for _i in primary_list:
                    _j=_i+1
                    if a[_i]>0.0 and a[_i]<a_hel_lim and e[_i]<1.0: # check helio orbits
                        # binary_record_file="../{}_analysis/binary_{:03d}_{:03d}.txt".format(run,_i,_j)
                        binary_record_file="../{}/binary_{:03d}_{:03d}.txt".format(run,_i,_j)
                        filename='{}'.format(binary_record_file)
                        print filename
                        with open(filename) as f:
                            dat = [line.rstrip() for line in f]
                        dat = [line.split() for line in dat]
                        dat = numpy.array(dat,dtype=float) # convert data list into array
                        # data format:
                        # t,a,e,inc,omega,Omega,f
                        a_bin.append(dat[-1,1])

            # create dataframes
            df_orb=pd.DataFrame()
            df_orb['a_hel(m)']=a_hel
            df_orb['e_hel']=e_hel
            df_orb['I_hel(rad)']=I_hel
            df_orb['r_hill(m)']=r_hill
            df_orb['a_bin(m)']=a_bin
            df_orb['a/r_hill']=df_orb['a_bin(m)']/df_orb['r_hill(m)']

            if k==0:
                # create dataframes
                df_orb_survivor=pd.DataFrame()
                df_orb_survivor['a_hel(m)']=a_hel_survivor
                df_orb_survivor['e_hel']=e_hel_survivor
                df_orb_survivor['I_hel(rad)']=I_hel_survivor
                df_orb_survivor['r_hill(m)']=r_hill_survivor
                df_orb_survivor['a_bin(m)']=a_bin_survivor
                df_orb_survivor['a/r_hill']=df_orb_survivor['a_bin(m)']/df_orb_survivor['r_hill(m)']

            # save dataframes
            if k==0:
                df_orb_survivor.to_csv(fname_df_orb_survivor,na_rep='NaN',sep="\t")
                df_orb.to_csv(fname_df_orb_initial,na_rep='NaN',sep="\t")
                print "save {}, {}".format(fname_df_orb_survivor,fname_df_orb_initial)
            else:
                df_orb.to_csv(fname_df_orb_final,na_rep='NaN',sep="\t")
                print "save {}".format(fname_df_orb_final)

        else: # binary data had already been loaded and saved as a dataframe, load this
            if k==0:
                df_orb=df_orb_initial
            else:
                df_orb=df_orb_final

        # Here we plot the initial (k=0) and final (k!=0) binary aei conditions
        # Small points for a/r_hill < wide_cut, large points for a/r_hill > wide_cut
        if k==0:
            col=colours[0]
        else:
            col=colours[1]

        print df_orb

        # tight binaries
        df_orb1=df_orb[df_orb['a/r_hill']<wide_cut]
        print "number of tight binaries = {}".format(len(df_orb1))
        ax1.scatter(df_orb1['a_hel(m)']/pf.AU,df_orb1['e_hel'],label="t={}Myr".format(t),s=s1,color=col)
        ax2.scatter(df_orb1['a_hel(m)']/pf.AU,df_orb1['I_hel(rad)']/(numpy.pi/180.0),label="N_binaries={}".format(len(a_hel)),s=s1,color=col)
        ax3.scatter(df_orb1['a_hel(m)']/pf.AU,df_orb1['e_hel'],label="t={}Myr".format(t),s=s1,color=col)
        ax4.scatter(df_orb1['a_hel(m)']/pf.AU,df_orb1['I_hel(rad)']/(numpy.pi/180.0),label="N_binaries={}".format(len(a_hel)),s=s1,color=col)
        # wide binaries
        df_orb2=df_orb[df_orb['a/r_hill']>=wide_cut]
        print "number of wide binaries = {}".format(len(df_orb2))
        ax1.scatter(df_orb2['a_hel(m)']/pf.AU,df_orb2['e_hel'],s=s2,color=col,label=None)
        ax2.scatter(df_orb2['a_hel(m)']/pf.AU,df_orb2['I_hel(rad)']/(numpy.pi/180.0),s=s2,color=col,label=None)
        ax3.scatter(df_orb2['a_hel(m)']/pf.AU,df_orb2['e_hel'],s=s2,color=col,label=None)
        ax4.scatter(df_orb2['a_hel(m)']/pf.AU,df_orb2['I_hel(rad)']/(numpy.pi/180.0),s=s2,color=col,label=None)

        # check numbers of bodies. Note, a binary object consists of two bodies
        print "number of bound binary objects, on closed helio orbits:",len(a_hel)*2
        print "number of single objects, on closed helio orbits:",len(a_hel_single)
        print "number of bound binary objects, on open helio orbits:",ejection_binary
        print "number of single objects, on open helio orbits:",ejection_single
        print "total number of bodies = {}".format(len(a_hel)*2+len(a_hel_single)+ejection_binary+ejection_single)

    # plot the initial conditions of the binaries that survive, for wide and tight binaries
    df_orb_survivor1=df_orb_survivor[df_orb_survivor['a/r_hill']<wide_cut]
    print len(df_orb_survivor1)
    ax1.scatter(df_orb_survivor1['a_hel(m)']/pf.AU,df_orb_survivor1['e_hel'],edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s1)
    ax2.scatter(df_orb_survivor1['a_hel(m)']/pf.AU,df_orb_survivor1['I_hel(rad)']/(numpy.pi/180.0),edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s1)
    ax3.scatter(df_orb_survivor1['a_hel(m)']/pf.AU,df_orb_survivor1['e_hel'],edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s1)
    ax4.scatter(df_orb_survivor1['a_hel(m)']/pf.AU,df_orb_survivor1['I_hel(rad)']/(numpy.pi/180.0),edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s1)
    df_orb_survivor2=df_orb_survivor[df_orb_survivor['a/r_hill']>=wide_cut]
    print len(df_orb_survivor2)
    ax1.scatter(df_orb_survivor2['a_hel(m)']/pf.AU,df_orb_survivor2['e_hel'],edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s2)
    ax2.scatter(df_orb_survivor2['a_hel(m)']/pf.AU,df_orb_survivor2['I_hel(rad)']/(numpy.pi/180.0),edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s2)
    ax3.scatter(df_orb_survivor2['a_hel(m)']/pf.AU,df_orb_survivor2['e_hel'],edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s2)
    ax4.scatter(df_orb_survivor2['a_hel(m)']/pf.AU,df_orb_survivor2['I_hel(rad)']/(numpy.pi/180.0),edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s2)

    # plot all single objects
    ax1.scatter(numpy.array(a_hel_single)/pf.AU,e_hel_single,c='k',alpha=0.5,edgecolors='none',marker='^')
    ax2.scatter(numpy.array(a_hel_single)/pf.AU,numpy.array(I_hel_single)/(numpy.pi/180.0),c='k',alpha=0.5,edgecolors='none',label="N_singles={}".format(len(a_hel_single)),marker='^')
    ax3.scatter(numpy.array(a_hel_single)/pf.AU,e_hel_single,c='k',alpha=0.5,edgecolors='none',marker='^')
    ax4.scatter(numpy.array(a_hel_single)/pf.AU,numpy.array(I_hel_single)/(numpy.pi/180.0),c='k',alpha=0.5,edgecolors='none',label="N_singles={}".format(len(a_hel_single)),marker='^')

    # plot a line to show evolution from initial to final positions
    for j in range(len(a_hel)):
        ax1.plot(numpy.array([a_hel_survivor[j],a_hel[j]])/pf.AU,[e_hel_survivor[j],e_hel[j]],color='k',alpha=0.2)
        ax2.plot(numpy.array([a_hel_survivor[j],a_hel[j]])/pf.AU,numpy.array([I_hel_survivor[j],I_hel[j]])/(numpy.pi/180.0),color='k',alpha=0.2)
        ax3.plot(numpy.array([a_hel_survivor[j],a_hel[j]])/pf.AU,[e_hel_survivor[j],e_hel[j]],color='k',alpha=0.2)
        ax4.plot(numpy.array([a_hel_survivor[j],a_hel[j]])/pf.AU,numpy.array([I_hel_survivor[j],I_hel[j]])/(numpy.pi/180.0),color='k',alpha=0.2)

    # draw lines approximating the cold classical belt
    a_min=40.0#*pf.AU
    a_max=46.0#*pf.AU
    q=36.0
    e_min=(1.0-(q/a_min))
    e_max=(1.0-(q/a_max))
    print "e_min = {}, e_max = {}".format(e_min,e_max)
    I_min=0.0
    I_max=5.0#*numpy.pi/180.0
    ax1.plot([a_min,a_min],[0.0,e_min],'k',zorder=0)
    ax1.plot([a_max,a_max],[0.0,e_max],'k',zorder=0)
    ax2.plot([a_min,a_min],[I_min,I_max],'k',zorder=0)
    ax2.plot([a_max,a_max],[I_min,I_max],'k',zorder=0)
    ax2.plot([a_min,a_max],[I_max,I_max],'k',zorder=0)
    ax3.plot([a_min,a_min],[0.0,e_min],'k',zorder=0)
    ax3.plot([a_max,a_max],[0.0,e_max],'k',zorder=0)
    ax4.plot([a_min,a_min],[I_min,I_max],'k',zorder=0)
    ax4.plot([a_max,a_max],[I_min,I_max],'k',zorder=0)
    ax4.plot([a_min,a_max],[I_max,I_max],'k',zorder=0)
    # plot the e limit, which is slightly curved
    n_plot=100
    a_lim=numpy.linspace(a_min,a_max,n_plot)
    e_lim=(1.0-(q/a_lim))
    #print a_lim
    #print e_lim
    ax1.plot(a_lim,e_lim,'k',zorder=0)
    ax3.plot(a_lim,e_lim,'k',zorder=0)

    # highlight the surviving binaries that lie in the cold classical region of interest
    df_orb['a_hel(AU)']=df_orb['a_hel(m)']/pf.AU
    df_orb['I_hel(d)']=df_orb['I_hel(rad)']/(numpy.pi/180.0)
    print "pre cuts\n",df_orb
    df_orb=df_orb[(df_orb['a_hel(AU)']>a_min) & (df_orb['a_hel(AU)']<a_max)]
    df_orb=df_orb[(df_orb['I_hel(d)']>I_min) & (df_orb['I_hel(d)']<I_max)]
    df_orb=df_orb[(df_orb['e_hel'])<(1.0-(q/df_orb['a_hel(AU)']))]
    print "post cuts\n",df_orb
    # df_orb['a_hel(AU)']=df_orb['a_hel(m)']/pf.AU
    # df_orb['I_hel(d)']=df_orb['I_hel(rad)']/(numpy.pi/180.0)
    # df_orb['q_hel(AU)']=df_orb['a_hel(AU)']*(1.0-df_orb['e_hel'])
    # print df_orb
    # df_orb=df_orb[(df_orb['I_hel(d)']>I_min) & (df_orb['I_hel(d)']<I_max)]
    # df_orb=df_orb[(df_orb['q_hel(AU)']<q)]
    # print df_orb

    # print "test e:"
    # print df_orb['e_hel']
    # print (1.0-(q/df_orb['a_hel(AU)']))
    # a_test=numpy.linspace(0,100)
    # ax1.plot(a_test,(1.0-(q/a_test)))

    if len(df_orb)==0:
        print 'There are no cold classicals binaries formed'
    else:
        print 'There are {} cold classical binaries formed'.format(len(df_orb))
        ax1.scatter(df_orb['a_hel(AU)'],df_orb['e_hel'],c='r',marker='x',label=None,s=CC_size,zorder=2)
        ax2.scatter(df_orb['a_hel(AU)'],df_orb['I_hel(d)'],c='r',marker='x',label=None,s=CC_size,zorder=2)
        ax3.scatter(df_orb['a_hel(AU)'],df_orb['e_hel'],c='r',marker='x',label=None,s=CC_size,zorder=2)
        ax4.scatter(df_orb['a_hel(AU)'],df_orb['I_hel(d)'],c='r',marker='x',label=None,s=CC_size,zorder=2)

    # repeat for single objects
    df_hel_single=pd.DataFrame(numpy.array([a_hel_single,e_hel_single,I_hel_single]).T,columns=['a_hel(m)','e_hel','I_hel(rad)'])
    df_hel_single['a_hel(AU)']=df_hel_single['a_hel(m)']/pf.AU
    df_hel_single['I_hel(d)']=df_hel_single['I_hel(rad)']/(numpy.pi/180.0)
    df_hel_single=df_hel_single[(df_hel_single['a_hel(AU)']>a_min) & (df_hel_single['a_hel(AU)']<a_max)]
    df_hel_single=df_hel_single[(df_hel_single['I_hel(d)']>I_min) & (df_hel_single['I_hel(d)']<I_max)]
    df_hel_single=df_hel_single[(df_hel_single['e_hel']<(1.0-(q/df_hel_single['a_hel(AU)'])))]

    # df_hel_single['a_hel(AU)']=df_hel_single['a_hel(m)']/pf.AU
    # df_hel_single['I_hel(d)']=df_hel_single['I_hel(rad)']/(numpy.pi/180.0)
    # df_hel_single['q_hel(AU)']=df_hel_single['a_hel(AU)']*(1.0-df_hel_single['e_hel'])
    # df_hel_single=df_hel_single[(df_hel_single['I_hel(d)']>I_min) & (df_hel_single['I_hel(d)']<I_max) & (df_hel_single['q_hel(AU)']>q)]

    if len(df_hel_single)==0:
        print 'There are no cold classicals singles formed'
    else:
        print 'There are {} cold classical singles formed'.format(len(df_hel_single))
        ax1.scatter(df_hel_single['a_hel(AU)'],df_hel_single['e_hel'],c='r',marker='x',label=None,s=CC_size)
        ax2.scatter(df_hel_single['a_hel(AU)'],df_hel_single['I_hel(d)'],c='r',marker='x',label=None,s=CC_size)
        ax3.scatter(df_hel_single['a_hel(AU)'],df_hel_single['e_hel'],c='r',marker='x',label=None,s=CC_size)
        ax4.scatter(df_hel_single['a_hel(AU)'],df_hel_single['I_hel(d)'],c='r',marker='x',label=None,s=CC_size)

    # # draw lines approximating the cold classical belt
    # a_min=42.0#*pf.AU
    # a_max=47.5#*pf.AU
    # e_min=0.14
    # e_max=0.24
    # I_min=0.0
    # I_max=6.0#*numpy.pi/180.0
    # ax1.plot([a_min,a_min],[0.0,e_min],'k',zorder=0)
    # ax1.plot([a_max,a_max],[0.0,e_max],'k',zorder=0)
    # ax1.plot([a_min,a_max],[e_min,e_max],'k',zorder=0)
    # ax2.plot([a_min,a_min],[I_min,I_max],'k',zorder=0)
    # ax2.plot([a_max,a_max],[I_min,I_max],'k',zorder=0)
    # ax2.plot([a_min,a_max],[I_max,I_max],'k',zorder=0)
    # ax3.plot([a_min,a_min],[0.0,e_min],'k',zorder=0)
    # ax3.plot([a_max,a_max],[0.0,e_max],'k',zorder=0)
    # ax3.plot([a_min,a_max],[e_min,e_max],'k',zorder=0)
    # ax4.plot([a_min,a_min],[I_min,I_max],'k',zorder=0)
    # ax4.plot([a_max,a_max],[I_min,I_max],'k',zorder=0)
    # ax4.plot([a_min,a_max],[I_max,I_max],'k',zorder=0)
    #
    # # highlight the surviving binaries that lie in the cold classical region of interest
    # df_orb['a_hel(AU)']=df_orb['a_hel(m)']/pf.AU
    # df_orb['I_hel(d)']=df_orb['I_hel(rad)']/(numpy.pi/180.0)
    # df_orb=df_orb[(df_orb['a_hel(AU)']>a_min) & (df_orb['a_hel(AU)']<a_max)]
    # df_orb=df_orb[(df_orb['I_hel(d)']>I_min) & (df_orb['I_hel(d)']<I_max)]
    # m=(e_max-e_min)/(a_max-a_min) # y = mx + c
    # c=e_min-(m*a_min)
    # df_orb=df_orb[(df_orb['e_hel']<((m*df_orb['a_hel(AU)'])+c))]
    # if len(df_orb)==0:
    #     print 'There are no cold classicals binaries formed'
    # else:
    #     print 'There are {} cold classical binaries formed'.format(len(df_orb))
    #     ax1.scatter(df_orb['a_hel(AU)'],df_orb['e_hel'],c='r',marker='x',label=None,s=CC_size)
    #     ax2.scatter(df_orb['a_hel(AU)'],df_orb['I_hel(d)'],c='r',marker='x',label=None,s=CC_size)
    #     ax3.scatter(df_orb['a_hel(AU)'],df_orb['e_hel'],c='r',marker='x',label=None,s=CC_size)
    #     ax4.scatter(df_orb['a_hel(AU)'],df_orb['I_hel(d)'],c='r',marker='x',label=None,s=CC_size)
    #
    # # repeat for single objects
    # df_hel_single=pd.DataFrame(numpy.array([a_hel_single,e_hel_single,I_hel_single]).T,columns=['a_hel','e_hel','I_hel'])
    # df_hel_single['a_hel']=df_hel_single['a_hel']/pf.AU
    # df_hel_single['I_hel']=df_hel_single['I_hel']/(numpy.pi/180.0)
    # df_hel_single=df_hel_single[(df_hel_single['a_hel']>a_min) & (df_hel_single['a_hel']<a_max)]
    # df_hel_single=df_hel_single[(df_hel_single['I_hel']>I_min) & (df_hel_single['I_hel']<I_max)]
    # m=(e_max-e_min)/(a_max-a_min)# y = mx + c
    # c=e_min-(m*a_min)
    # df_hel_single=df_hel_single[(df_hel_single['e_hel']<((m*df_hel_single['a_hel'])+c))]
    # if len(df_hel_single)==0:
    #     print 'There are no cold classicals singles formed'
    # else:
    #     print 'There are {} cold classical singles formed'.format(len(df_hel_single))
    #     ax1.scatter(df_hel_single['a_hel'],df_hel_single['e_hel'],c='r',marker='x',label=None,s=CC_size)
    #     ax2.scatter(df_hel_single['a_hel'],df_hel_single['I_hel'],c='r',marker='x',label=None,s=CC_size)
    #     ax3.scatter(df_hel_single['a_hel'],df_hel_single['e_hel'],c='r',marker='x',label=None,s=CC_size)
    #     ax4.scatter(df_hel_single['a_hel'],df_hel_single['I_hel'],c='r',marker='x',label=None,s=CC_size)

    # Plot an axes that shows a zoom in of the CC region
    # Set zoom axis limits
    a_nudge=1.0
    e_nudge=0.05
    I_nudge=1.0
    x1=a_min-a_nudge
    x2=a_max+a_nudge
    y1=0.0
    y2=e_max+e_nudge
    ax3.set_xlim(x1,x2)
    ax3.set_ylim(y1,y2)

    _x1=a_min-a_nudge
    _x2=a_max+a_nudge
    _y1=0.0
    _y2=I_max+I_nudge
    ax4.set_xlim(_x1,_x2)
    ax4.set_ylim(_y1,_y2)

    #add Rectangle and lines from ax3 to ax6
    ax1.add_patch(patches.Rectangle((x1, y1),x2-x1,y2-y1,fill=False,color='r'))
    xy1 = (x1,y1)
    xy2 = (x1,y2)
    xy3 = (x2,y1)
    xy4 = (x2,y2)
    con = ConnectionPatch(xyA=xy2, xyB=xy2, coordsA="data", coordsB="data",axesA=ax3, axesB=ax1, color="red")
    ax3.add_artist(con)
    con = ConnectionPatch(xyA=xy1, xyB=xy1, coordsA="data", coordsB="data",axesA=ax3, axesB=ax1, color="red")
    ax3.add_artist(con)

    ax2.add_patch(patches.Rectangle((_x1, _y1),_x2-_x1,_y2-_y1,fill=False,color='r'))
    xy1 = (_x1,_y1)
    xy2 = (_x1,_y2)
    xy3 = (_x2,_y1)
    xy4 = (_x2,_y2)
    con = ConnectionPatch(xyA=xy2, xyB=xy2, coordsA="data", coordsB="data",axesA=ax4, axesB=ax2, color="red")
    ax4.add_artist(con)
    con = ConnectionPatch(xyA=xy1, xyB=xy1, coordsA="data", coordsB="data",axesA=ax4, axesB=ax2, color="red")
    ax4.add_artist(con)

    ax1.legend()
    ax2.legend()

    # save_path = '../{}_analysis/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
    save_path = '../{}/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)

    print save_path
    pyplot.savefig(save_path, bbox_inches='tight')
    # pyplot.show()
    pyplot.close()
    # break

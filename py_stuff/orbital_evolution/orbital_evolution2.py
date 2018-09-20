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
import pandas as pd

# path='/Users/jrobinson/bin_migration/'
# # run='200_dump'
# run='prone_121'

# path='/Users/jrobinson/bin_migration/migration'
# run='binSim'
# # run='binSim_20_0'

path='/Users/jrobinson/bin_migration/py_stuff'
dirs=next(os.walk(path))[1]
dirs = [d for d in dirs if 'prone' in d]

#dirs=['prone_1335_timesteptest0.75_analysis']
dirs=['prone_1000_timesteptest0.25','prone_1000_timesteptest0.75']

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
n_bins=10 # bins for inclination histogram
#-------------------------------------------------------------------------------

for run in dirs:

    # if int(run.split('_')[1])<1240:
    #     continue

    run_path='{}/{}'.format(path,run)
    files=next(os.walk(run_path))[2]
    files = [ fi for fi in files if fi.endswith('.dump')]
    files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers
    print run_path

    # calculate only the first and last state
    files=[files[0],files[-1]]

    # primary_fname='../{}_analysis/{}_primary_list.txt'.format(run,run)
    # single_fname='../{}_analysis/{}_single_list.txt'.format(run,run)

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
        e_bin=[]
        I_bin=[]

        if k==0:
            a_hel_survivor=[]
            e_hel_survivor=[]
            I_hel_survivor=[]
            r_hill_survivor=[]
            a_bin_survivor=[]
            e_bin_survivor=[]
            I_bin_survivor=[]
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
                    binary_record_file="../{}/binary_{:03d}_{:03d}_short.txt".format(run,_i,_j)
                    filename='{}'.format(binary_record_file)
                    print filename
                    with open(filename) as f:
                        dat = [line.rstrip() for line in f]
                    dat = [line.split() for line in dat]
                    dat = numpy.array(dat,dtype=float) # convert data list into array
                    # data format:
                    # t,a,e,inc,omega,Omega,f
                    a_bin.append(dat[0,1])
                    e_bin.append(dat[0,2])
                    I_bin.append(dat[0,3])

                    if no_survive==1:
                        _i+=2
                        continue

                    if _i in primary_list: # keep only objects that survive
                        if a[_i]>0.0 and e[_i]<1.0: # record only heliocentrically bound binaries
                            a_bin_survivor.append(dat[0,1])
                            e_bin_survivor.append(dat[0,2])
                            I_bin_survivor.append(dat[0,3])

                    _i+=2

            if k!=0:

                if no_survive==1:
                    continue

                for _i in primary_list:
                    _j=_i+1
                    if a[_i]>0.0 and a[_i]<a_hel_lim and e[_i]<1.0: # check helio orbits
                        # binary_record_file="../{}_analysis/binary_{:03d}_{:03d}.txt".format(run,_i,_j)
                        binary_record_file="../{}/binary_{:03d}_{:03d}_short.txt".format(run,_i,_j)
                        filename='{}'.format(binary_record_file)
                        print filename
                        with open(filename) as f:
                            dat = [line.rstrip() for line in f]
                        dat = [line.split() for line in dat]
                        dat = numpy.array(dat,dtype=float) # convert data list into array
                        # data format:
                        # t,a,e,inc,omega,Omega,f
                        a_bin.append(dat[-1,1])
                        e_bin.append(dat[-1,2])
                        I_bin.append(dat[-1,3])

            # create dataframes
            df_orb=pd.DataFrame()
            df_orb['a_hel(m)']=a_hel
            df_orb['e_hel']=e_hel
            df_orb['I_hel(rad)']=I_hel
            df_orb['r_hill(m)']=r_hill
            df_orb['a_bin(m)']=a_bin
            df_orb['e_bin']=e_bin
            df_orb['I_bin(rad)']=I_bin
            df_orb['a/r_hill']=df_orb['a_bin(m)']/df_orb['r_hill(m)']

            if k==0:
                # create dataframes
                df_orb_survivor=pd.DataFrame()
                df_orb_survivor['a_hel(m)']=a_hel_survivor
                df_orb_survivor['e_hel']=e_hel_survivor
                df_orb_survivor['I_hel(rad)']=I_hel_survivor
                df_orb_survivor['r_hill(m)']=r_hill_survivor
                df_orb_survivor['a_bin(m)']=a_bin_survivor
                df_orb_survivor['e_bin']=e_bin_survivor
                df_orb_survivor['I_bin(rad)']=I_bin_survivor
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
                a_bin=df_orb_initial['a_bin(m)']
                e_bin=df_orb_initial['e_bin']
                I_bin=df_orb_initial['I_bin(rad)']
            else:
                df_orb=df_orb_final
                a_bin=df_orb_final['a_bin(m)']
                e_bin=df_orb_final['e_bin']
                I_bin=df_orb_final['I_bin(rad)']
                a_bin_survivor=df_orb_survivor['a_bin(m)']
                e_bin_survivor=df_orb_survivor['e_bin']
                I_bin_survivor=df_orb_survivor['I_bin(rad)']


        # Here we plot the initial (k=0) and final (k!=0) binary aei conditions
        # Small points for a/r_hill < wide_cut, large points for a/r_hill > wide_cut
        if k==0:
            col=colours[0]
        else:
            col=colours[1]

        # print df_orb

        # check numbers of bodies. Note, a binary object consists of two bodies
        print "number of bound binary objects, on closed helio orbits:",len(a_hel)*2
        print "number of single objects, on closed helio orbits:",len(a_hel_single)
        print "number of bound binary objects, on open helio orbits:",ejection_binary
        print "number of single objects, on open helio orbits:",ejection_single
        print "total number of bodies = {}".format(len(a_hel)*2+len(a_hel_single)+ejection_binary+ejection_single)

        # tight binaries
        df_orb1=df_orb[df_orb['a/r_hill']<wide_cut]
        print "number of tight binaries = {}".format(len(df_orb1))
        ax1.scatter(df_orb1['a_bin(m)'],df_orb1['e_bin'],label="t={}Myr".format(t),s=s1,color=col)
        ax2.scatter(df_orb1['a_bin(m)'],df_orb1['I_bin(rad)']/(numpy.pi/180.0),label="N_binaries={}".format(len(a_bin)),s=s1,color=col)

        # wide binaries
        df_orb2=df_orb[df_orb['a/r_hill']>=wide_cut]
        print "number of wide binaries = {}".format(len(df_orb2))
        ax1.scatter(df_orb2['a_bin(m)'],df_orb2['e_bin'],label=None,s=s2,color=col)
        ax2.scatter(df_orb2['a_bin(m)'],df_orb2['I_bin(rad)']/(numpy.pi/180.0),label=None,s=s2,color=col)

        I_bin=numpy.array(I_bin)/(numpy.pi/180.0) # convert to degrees

        # ax1.scatter(a_bin,e_bin,label="t={}Myr".format(t),zorder=2)
        # ax2.scatter(a_bin,I_bin,label="N_binaries={}".format(len(a_bin)),zorder=2)

        # colours=pyplot.rcParams['axes.prop_cycle'].by_key()['color']

        if k==0:

            print len(a_bin)
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

            # Add inclination histogram for final states
            center, hist, width=pf.hist_dist(I_bin,n_bins)
            hist=hist*hist_scale_I
            # plot bars on the left
            ax2.barh(y=center, height=width, width=hist, align='center', color=colours[1], edgecolor='k',zorder=0,alpha=0.5)#,label='t={}Myr'.format(t))


    I_bin_survivor=numpy.array(I_bin_survivor)/(numpy.pi/180.0) # convert to degrees

    print len(a_bin),len(a_bin_survivor)
    # # plot the initial conditions of the binaries that survive
    # ax1.scatter(a_bin_survivor,e_bin_survivor,edgecolor='k',facecolors='none',zorder=2)
    # ax2.scatter(a_bin_survivor,I_bin_survivor,edgecolor='k',facecolors='none',zorder=2)

    # plot the initial conditions of the binaries that survive, for wide and tight binaries
    df_orb_survivor1=df_orb_survivor[df_orb_survivor['a/r_hill']<wide_cut]
    print "len df_orb_survivor1={}".format(len(df_orb_survivor1))
    ax1.scatter(df_orb_survivor1['a_bin(m)'],df_orb_survivor1['e_bin'],edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s1)
    ax2.scatter(df_orb_survivor1['a_bin(m)'],df_orb_survivor1['I_bin(rad)']/(numpy.pi/180.0),edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s1)
    df_orb_survivor2=df_orb_survivor[df_orb_survivor['a/r_hill']>=wide_cut]
    print "len df_orb_survivor2={}".format(len(df_orb_survivor2))
    ax1.scatter(df_orb_survivor2['a_bin(m)'],df_orb_survivor2['e_bin'],edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s2)
    ax2.scatter(df_orb_survivor2['a_bin(m)'],df_orb_survivor2['I_bin(rad)']/(numpy.pi/180.0),edgecolor='k',facecolors='none',lw=survivor_lw,label=None,s=s2)

    # plot a line to show evolution from initial to final positions
    if k==1:
        for j in range(len(a_bin)):
            ax1.plot([a_bin_survivor[j],a_bin[j]],[e_bin_survivor[j],e_bin[j]],color='k',alpha=0.2,zorder=2)
            ax2.plot([a_bin_survivor[j],a_bin[j]],[I_bin_survivor[j],I_bin[j]],color='k',alpha=0.2,zorder=2)


    ax1.legend()
    ax2.legend()

    # save_path = '../{}_analysis/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
    save_path = '../{}/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)

    print save_path
    pyplot.savefig(save_path, bbox_inches='tight')
    # pyplot.show()
    pyplot.close()
    # break

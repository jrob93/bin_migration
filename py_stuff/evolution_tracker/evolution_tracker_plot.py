'''
In this script we bring together elements of other scripts.
Using a series of plots we show the evolution of the heliocentric elements of simulation particles,
as binaries and singles, and we also show how the mutual binary orbits evolve, as a function of simulation time.

Add a horizontal line for
'''
import numpy
import matplotlib.pyplot as pyplot
import os
from natsort import natsorted, ns
import subprocess
import matplotlib.gridspec as gridspec
import py_func as pf
import pandas as pd
import rebound

# set up initial variables
N_planet=4 # number of planets, used to skip straight to binary data
bin_mass=0.5e18 # kg, mass of a binary component
mi=bin_mass
mj=mi
a_hel_lim=1000*pf.AU # m, maximum semimajor axis of interest, beyond this we cut any heliocentrically bound particles
file_sample=50 # number of images to sample the range of dump files
a_nep=27.0 # neptune semi major axis of interest, AU

# load planet data for modern values of planet positions
df_planet=pd.read_csv("../gen_planet_data.txt",sep="\t",index_col=0)
print df_planet

# establish path for the runs
path='/Users/jrobinson/bin_migration/py_stuff'
# dirs=next(os.walk(path))[1]
# dirs = [d for d in dirs if 'prone' in d]
dirs=['prone_1200_analysis']
print dirs

for run in dirs:

    subprocess.Popen(["mkdir","../{}/{}_figs".format(run,os.path.basename(__file__).split('.')[0])]) # make directory to store plots for that run

    #---------------------------------------------------------------------------
    # Load the planet orbital elements data
    for i in range(1,N_planet+1):
        _i=i-1
        print _i
        planet_record_file="../{}/planet_{:03d}.txt".format(run,i)
        filename='{}'.format(planet_record_file)
        print filename
        try:
            with open(filename) as f:
                dat = [line.rstrip() for line in f]
            dat = [line.split() for line in dat]
            dat = numpy.array(dat,dtype=float) # convert data list into array
        except:
            print 'error, skip'
            break

        if _i==0: # initialise data array

            with open("../{}/0.dump".format(run)) as f:
                _dat = [line.rstrip() for line in f]
            header = _dat[0] # header line
            print header
            dat_len=len(dat[:,0])
            a_planet=numpy.zeros((dat_len,N_planet))
            e_planet=numpy.zeros((dat_len,N_planet))
            I_planet=numpy.zeros((dat_len,N_planet))

        # data format:
        # t,a,e,inc,omega,Omega,f

        a_planet[:,_i]=dat[:,1]
        e_planet[:,_i]=dat[:,2]
        I_planet[:,_i]=dat[:,3]

        if _i==0:
            time=dat[:,0]
    #---------------------------------------------------------------------------
    # Set up the figure
    fig = pyplot.figure()
    fig.set_size_inches(15, 10)
    gs = gridspec.GridSpec(3, 2, height_ratios=[1,1,0.5])
    ax1 = pyplot.subplot(gs[0,0])
    ax2 = pyplot.subplot(gs[1,0],sharex=ax1)
    ax3 = pyplot.subplot(gs[0,1])
    ax4 = pyplot.subplot(gs[1,1],sharex=ax3)
    ax5 = pyplot.subplot(gs[2,:])
    ax1.set_ylabel('e_hel')
    ax1.set_xlabel('a_hel (AU)')
    ax2.set_ylabel('I_hel (degress)')
    ax2.set_xlabel('a_hel (AU)')
    ax3.set_ylabel('e_bin')
    ax3.set_xlabel('a_bin (m)')
    ax4.set_ylabel('I_bin (degress)')
    ax4.set_xlabel('a_bin (m)')
    ax5.set_xlabel('t (Myr)')
    ax5.set_ylabel('a_hel (AU)')

    print run
    params=numpy.array(header.split()[6:]).astype(float)
    print params
    fig.suptitle('{}: N_a_i={}AU, N_a_f={}AU, tau_a={}Myr, tau_e={}Myr, tau_i={}Myr, N_e_o={}, N_e_i={}deg, tmax={}Myr, dt={}d'.format(
    run,params[0],params[1],params[2]/1e6,params[3]/1e6,params[4]/1e6,params[5],params[6],params[7]/1e6,params[8]))

    colours=pf.pyplot_colours
    #---------------------------------------------------------------------------
    # plot the planet orbital elements as a function of time
    planet_name=["Jupiter","Saturn","Uranus","Neptune"]
    for j in range(N_planet):
        ax5.plot(time,a_planet[:,j],label=planet_name[j],zorder=0)
    # also plot q and Q
    for j in range(N_planet):
        q=a_planet[:,j]*(1.0-e_planet[:,j])
        Q=a_planet[:,j]*(1.0+e_planet[:,j])
        ax5.plot(time,q,linestyle=':',color=colours[j],alpha=0.3,zorder=0)
        ax5.plot(time,Q,linestyle=':',color=colours[j],alpha=0.3,zorder=0)
    # plot current planet data as reference points
    for k in range(N_planet):
        k=k+5
        ax5.scatter(time[-1],df_planet['a(m)'].iloc[k]/pf.AU,zorder=2,edgecolor='k')
    # ax1.legend()
    ax5.axhline(a_nep,time[0],time[-1],color='k',alpha=0.5)
    #---------------------------------------------------------------------------
    # draw lines approximating the cold classical belt
    a_min=40.0#*pf.AU
    a_max=46.0#*pf.AU
    q=36.0
    e_min=(1.0-(q/a_min))
    e_max=(1.0-(q/a_max))
    print "e_min = {}, e_max = {}".format(e_min,e_max)
    I_min=0.0
    I_max=5.0#*numpy.pi/180.0
    ax1.plot([a_min,a_min],[0.0,e_min],'k',zorder=2)
    ax1.plot([a_max,a_max],[0.0,e_max],'k',zorder=2)
    ax2.plot([a_min,a_min],[I_min,I_max],'k',zorder=2)
    ax2.plot([a_max,a_max],[I_min,I_max],'k',zorder=2)
    ax2.plot([a_min,a_max],[I_max,I_max],'k',zorder=2)
    # plot the e limit, which is slightly curved
    n_plot=100
    a_lim=numpy.linspace(a_min,a_max,n_plot)
    e_lim=(1.0-(q/a_lim))
    ax1.plot(a_lim,e_lim,'k',zorder=2)
    #---------------------------------------------------------------------------
    # load the binary data, file by file
    load_path="/Volumes/Mirkwood/bin_migration/data/{}".format("_".join(run.split('_')[:2]))
    print load_path
    files=next(os.walk(load_path))[2]
    files = [ fi for fi in files if fi.endswith('.dump')]
    files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers
    files=files[::len(files)/file_sample] # only sample some of the dump files
    print files

    # record the number of objects ejected
    ejection_particle=[]

    for fi in files:
        print fi
        with open("{}/{}".format(load_path,fi)) as f:
            dat = [line.rstrip() for line in f]
        dat = [line.split() for line in dat]
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

        # search binary pairs and record binaries and singles
        a_hel=[]
        e_hel=[]
        I_hel=[]
        a_hel_single=[]
        e_hel_single=[]
        I_hel_single=[]
        a_bin=[]
        e_bin=[]
        I_bin=[]

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
                if a[i]>0.0 and e[i]<1.0 and a[i]<a_hel_lim: # consider if bound heliocentrically
                    a_hel.append(a[i])
                    e_hel.append(e[i])
                    I_hel.append(I[i])
                    a_bin.append(orbit.a)
                    e_bin.append(orbit.e)
                    I_bin.append(orbit.inc)
                else: # binary is heliocentrically unbound
                    ejection_particle.append(i)
                    ejection_particle.append(j)

            # record the single objects (which are still bound to the solar system)
            else:
                if a[i]>0.0 and e[i]<1.0 and a[i]<a_hel_lim:
                    a_hel_single.append(a[i])
                    e_hel_single.append(e[i])
                    I_hel_single.append(I[i])
                else:
                    ejection_particle.append(i)
                if a[j]>0.0 and e[j]<1.0 and a[j]<a_hel_lim:
                    a_hel_single.append(a[j])
                    e_hel_single.append(e[j])
                    I_hel_single.append(I[j])
                else:
                    ejection_particle.append(j)

            i+=2
        #---------------------------------------------------------------------------
        # Check numbers of particles: either binary, single, or ejected
        print "number of binaries = {}".format(len(a_bin))
        print "number of objects in binary pairs = {}".format(len(a_bin)*2)
        print "number of singles = {}".format(len(a_hel_single))
        print "ejected particles = {}".format(len(numpy.unique(ejection_particle)))
        print "total particles = {}".format(len(a_hel_single)+(len(a_bin)*2)+len(numpy.unique(ejection_particle)))
        #---------------------------------------------------------------------------
        # Unit conversions for plotting
        a_hel=numpy.array(a_hel)/pf.AU
        I_hel=numpy.array(I_hel)*(180.0/numpy.pi)
        a_hel_single=numpy.array(a_hel_single)/pf.AU
        I_hel_single=numpy.array(I_hel_single)*(180.0/numpy.pi)
        I_bin=numpy.array(I_bin)*(180.0/numpy.pi)
        #---------------------------------------------------------------------------
        # Plot stuff for this timestep
        l1=ax5.axvline(t,0,30,color='k') # line to track simulation time
        s1_0=ax1.scatter(a_hel,e_hel,color=colours[0])
        s2_0=ax2.scatter(a_hel,I_hel,color=colours[0])
        s1_1=ax1.scatter(a_hel_single,e_hel_single,color='k',alpha=0.5,marker='^',edgecolors='none')
        s2_1=ax2.scatter(a_hel_single,I_hel_single,color='k',alpha=0.5,marker='^',edgecolors='none')
        s3=ax3.scatter(a_bin,e_bin,color=colours[0])
        s4=ax4.scatter(a_bin,I_bin,color=colours[0])
        #---------------------------------------------------------------------------
        # Find objects that lie in the cold classical region:
        # binary objects
        df_hel=pd.DataFrame(numpy.array([a_hel,e_hel,I_hel]).T,columns=['a_hel(AU)','e_hel','I_hel(d)'])
        df_hel=df_hel[(df_hel['a_hel(AU)']>a_min) & (df_hel['a_hel(AU)']<a_max)]
        df_hel=df_hel[(df_hel['I_hel(d)']>I_min) & (df_hel['I_hel(d)']<I_max)]
        df_hel=df_hel[(df_hel['e_hel']<(1.0-(q/df_hel['a_hel(AU)'])))]
        if len(df_hel)>0:
            print 'There are {} cold classical binaries formed'.format(len(df_hel))
            s1_2=ax1.scatter(df_hel['a_hel(AU)'],df_hel['e_hel'],c='r',marker='x')
            s2_2=ax2.scatter(df_hel['a_hel(AU)'],df_hel['I_hel(d)'],c='r',marker='x')
        # repeat for single objects
        df_hel_single=pd.DataFrame(numpy.array([a_hel_single,e_hel_single,I_hel_single]).T,columns=['a_hel(AU)','e_hel','I_hel(d)'])
        df_hel_single=df_hel_single[(df_hel_single['a_hel(AU)']>a_min) & (df_hel_single['a_hel(AU)']<a_max)]
        df_hel_single=df_hel_single[(df_hel_single['I_hel(d)']>I_min) & (df_hel_single['I_hel(d)']<I_max)]
        df_hel_single=df_hel_single[(df_hel_single['e_hel']<(1.0-(q/df_hel_single['a_hel(AU)'])))]
        if len(df_hel_single)>0:
            print 'There are {} cold classical singles formed'.format(len(df_hel_single))
            s1_3=ax1.scatter(df_hel_single['a_hel(AU)'],df_hel_single['e_hel'],c='r',marker='x')
            s2_3=ax2.scatter(df_hel_single['a_hel(AU)'],df_hel_single['I_hel(d)'],c='r',marker='x')
        #---------------------------------------------------------------------------
        # Find semimajor axis of Neptune
        t_mask=numpy.where(time==t)
        print "sim time = {} Myr, a_neptune = {} AU".format(time[t_mask],a_planet[:,N_planet-1][t_mask])
        #---------------------------------------------------------------------------
        # save and clear the figure for the next timestep
        save_path = '../{}/{}_figs/{}_{}_{:07d}.png'.format(run,os.path.basename(__file__).split('.')[0],os.path.basename(__file__).split('.')[0],run,int(fi.split('/')[-1].split('.')[0]))
        print save_path
        pyplot.savefig(save_path, bbox_inches='tight')
        # pyplot.show()
        # pyplot.close()
        # break
        # Clear the plots
        for ax_obj in [l1,s1_0,s2_0,s1_1,s2_1,s3,s4]:
            ax_obj.remove()
        if len(df_hel)>0:
            s1_2.remove()
            s2_2.remove()
        if len(df_hel_single)>0:
            s1_3.remove()
            s2_3.remove()
        #---------------------------------------------------------------------------

        # break

# pyplot.show()
pyplot.close()

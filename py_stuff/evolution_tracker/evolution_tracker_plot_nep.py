'''
In this script we bring together elements of other scripts.
Using a series of plots we show the evolution of the heliocentric elements of simulation particles,
as binaries and singles, and we also show how the mutual binary orbits evolve, as a function of simulation time.

Here we make a single figure to show the activity around a=27 for neptune

Note, we have assumed that the heliocentric orbit of the bound binary object is equal to the orbit of the primary
'''
import numpy
import matplotlib.pyplot as pyplot
import os
from natsort import natsorted, ns
import subprocess
import matplotlib
import matplotlib.gridspec as gridspec
import py_func as pf
import pandas as pd
import rebound
from matplotlib.patches import Ellipse

def CC_bounds(df_hel):
    df_hel=df_hel[(df_hel['a_hel(AU)']>a_min) & (df_hel['a_hel(AU)']<a_max)]
    df_hel=df_hel[(df_hel['I_hel(d)']>I_min) & (df_hel['I_hel(d)']<I_max)]
    df_hel=df_hel[(df_hel['e_hel']<(1.0-(q/df_hel['a_hel(AU)'])))]
    return df_hel

# set up initial variables
N_planet=4 # number of planets, used to skip straight to binary data
bin_mass=0.5e18 # kg, mass of a binary component
mi=bin_mass
mj=mi
a_hel_lim=1000*pf.AU # m, maximum semimajor axis of interest, beyond this we cut any heliocentrically bound particles
# file_sample=5 # number of images to sample the range of dump files
a_nep=27.5 # neptune semi major axis of interest, AU
wide_cut=0.05 # a/r_hill, the normalised binary separation above which we consider the binary to be 'wide'
CC_size=75 # marker size for crosses to denote the point is in the Cold Classical region
s1=25 # marker size for tight binaries
s2=90 # marker size for wide binaries
e_bins=[0.0,0.1,0.3,1.0] # define the bins for eccentricty

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

    subprocess.Popen(["mkdir","../analysis_dirs/{}/{}_figs".format(run,os.path.basename(__file__).split('.')[0])]) # make directory to store plots for that run

    #---------------------------------------------------------------------------
    # Load the planet orbital elements data
    for i in range(1,N_planet+1):
        _i=i-1
        print _i
        planet_record_file="../analysis_dirs/{}/planet_{:03d}.txt".format(run,i)
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

            with open("../analysis_dirs/{}/0.dump".format(run)) as f:
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
    # ax2 = pyplot.subplot(gs[1,0],sharex=ax1)
    ax2 = pyplot.subplot(gs[1,0])
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
    # files=files[::len(files)/file_sample] # only sample some of the dump files
    # print files

    #---------------------------------------------------------------------------
    # find file when a_neptune is closest to 27 AU
    a_nep_time=[]
    a_neptune=a_planet[:,N_planet-1]
    nep_i=numpy.abs(a_neptune-a_nep).argmin()
    # for _nep_i in [nep_i-1,nep_i,nep_i+1]:
    for _nep_i in [nep_i]:
        print time[_nep_i],a_neptune[_nep_i]
        a_nep_time.append(time[_nep_i])
    a_nep_ind=(numpy.array(a_nep_time)*1000).astype(int) # find dump file index
    a_nep_ind=numpy.insert(a_nep_ind,0,0)
    #---------------------------------------------------------------------------

    # record the number of objects ejected
    ejection_particle=[]

    col_i=0
    for fi in files:
        # print fi
        if int(fi.split('.')[0]) not in a_nep_ind:
            continue
        with open("{}/{}".format(load_path,fi)) as f:
            dat = [line.rstrip() for line in f]
        dat = [line.split() for line in dat]
        t=float(dat[0][3]) # timestamp
        # if (t!=0) and (t not in a_nep_time):
        #     continue
        print fi,t

        dat = dat[1:] # strip header line
        dat = numpy.array(dat,dtype=float) # convert data list into array

        # data format:
        # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0

        a=dat[:,1]*pf.AU # m
        e=dat[:,2]
        I=dat[:,3] # radians
        pos=dat[:,7:10]*pf.AU # m
        vel=dat[:,10:13]*pf.AU/((24.0*60*60*365.25)/(2.0*numpy.pi)) # m/s

        #---------------------------------------------------------------------------

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
        ind_bin=[]
        r_hill=[]

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
                    ind_bin.append(i)
                    r_hill.append(pow((mi+mj)/(3.0*pf.M_sun),1.0/3.0)*(((a[i]*mi)+(a[j]*mj))/(mi+mj)))

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
        # Make dataframes for the orbits
        df_orb_all=pd.DataFrame()
        df_orb_all['a_hel(m)']=a_hel
        df_orb_all['e_hel']=e_hel
        df_orb_all['I_hel(rad)']=I_hel
        df_orb_all['r_hill(m)']=r_hill
        df_orb_all['a_bin(m)']=a_bin
        df_orb_all['e_bin']=e_bin
        df_orb_all['I_bin(rad)']=I_bin
        df_orb_all['a/r_hill']=df_orb_all['a_bin(m)']/df_orb_all['r_hill(m)']
        # print df_orb_all
        #---------------------------------------------------------------------------
        # define a dataframe to store the orbits for connecting start and end points of orbits
        df_orb=pd.DataFrame(numpy.array([numpy.array(a)/pf.AU,e,numpy.array(I)*(180.0/numpy.pi)]).T,columns=['a_hel(AU)','e_hel','I_hel(d)'])
        df_bin=pd.DataFrame(numpy.array([ind_bin,a_bin,e_bin,numpy.array(I_bin)*(180.0/numpy.pi)]).T,columns=['i','a_bin(m)','e_bin','I_bin(d)'])
        df_bin['i']=df_bin['i'].astype(int)
        if col_i==0:
            df_orb_0=df_orb
            df_bin_0=df_bin
        #---------------------------------------------------------------------------
        # add connecting lines
        if col_i>0:
            # connect heliocentric orbits
            for l in range(N_planet+1,len(df_orb_0),2):
                a_orb_i=df_orb.iloc[l]['a_hel(AU)']
                e_orb_i=df_orb.iloc[l]['e_hel']
                I_orb_i=df_orb.iloc[l]['I_hel(d)']
                a_orb_j=df_orb.iloc[l+1]['a_hel(AU)']
                e_orb_j=df_orb.iloc[l+1]['e_hel']
                I_orb_j=df_orb.iloc[l+1]['I_hel(d)']
                # we always connect the i'th particles
                ax1.plot([df_orb_0.iloc[l]['a_hel(AU)'],a_orb_i],[df_orb_0.iloc[l]['e_hel'],e_orb_i],color='k',alpha=0.5)
                ax2.plot([df_orb_0.iloc[l]['a_hel(AU)'],a_orb_i],[df_orb_0.iloc[l]['I_hel(d)'],I_orb_i],color='k',alpha=0.5)
                # only connect the i to j particles for single objects
                if l not in ind_bin:
                    ax1.plot([df_orb_0.iloc[l+1]['a_hel(AU)'],a_orb_j],[df_orb_0.iloc[l+1]['e_hel'],e_orb_j],color='k',alpha=0.5)
                    ax2.plot([df_orb_0.iloc[l+1]['a_hel(AU)'],a_orb_j],[df_orb_0.iloc[l+1]['I_hel(d)'],I_orb_j],color='k',alpha=0.5)
            # connect binary orbits
            inds=numpy.array(df_bin['i'])
            inds_0=numpy.array(df_bin_0['i'])
            for l in inds:
                if l in inds_0:
                    bin_0=df_bin_0[df_bin_0['i']==l]
                    a_bin_0=float(bin_0['a_bin(m)'])
                    e_bin_0=float(bin_0['e_bin'])
                    I_bin_0=float(bin_0['I_bin(d)'])
                    bin_orb=df_bin[df_bin['i']==l]
                    a_bin_orb=float(bin_orb['a_bin(m)'])
                    e_bin_orb=float(bin_orb['e_bin'])
                    I_bin_orb=float(bin_orb['I_bin(d)'])
                    ax3.plot([a_bin_0,a_bin_orb],[e_bin_0,e_bin_orb],color='k',alpha=0.5)
                    ax4.plot([a_bin_0,a_bin_orb],[I_bin_0,I_bin_orb],color='k',alpha=0.5)
        #---------------------------------------------------------------------------
        # Check numbers of particles: either binary, single, or ejected
        print "number of binaries = {}".format(len(a_bin))
        print "number of objects in binary pairs = {}".format(len(a_bin)*2)
        print "number of singles = {}".format(len(a_hel_single))
        print "ejected particles = {}".format(len(numpy.unique(ejection_particle)))
        print "total particles = {}".format(len(a_hel_single)+(len(a_bin)*2)+len(numpy.unique(ejection_particle)))
        # #---------------------------------------------------------------------------
        # # Unit conversions for plotting
        # a_hel=numpy.array(a_hel)/pf.AU
        # I_hel=numpy.array(I_hel)*(180.0/numpy.pi)
        # a_hel_single=numpy.array(a_hel_single)/pf.AU
        # I_hel_single=numpy.array(I_hel_single)*(180.0/numpy.pi)
        # I_bin=numpy.array(I_bin)*(180.0/numpy.pi)
        # #---------------------------------------------------------------------------
        # # Plot stuff for this timestep
        # l1=ax5.axvline(t,0,30,color='k') # line to track simulation time
        # s1_0=ax1.scatter(a_hel,e_hel,color=colours[col_i])#,alpha=0.5)
        # s2_0=ax2.scatter(a_hel,I_hel,color=colours[col_i])#,alpha=0.5)
        # s1_1=ax1.scatter(a_hel_single,e_hel_single,color=colours[col_i],alpha=0.5,marker='^',edgecolors='none')
        # s2_1=ax2.scatter(a_hel_single,I_hel_single,color=colours[col_i],alpha=0.5,marker='^',edgecolors='none')
        # s3=ax3.scatter(a_bin,e_bin,color=colours[col_i],label="t={}Myr, a_neptune={}AU".format(t,float(a_neptune[numpy.where(time==t)])))#,alpha=0.5)
        # s4=ax4.scatter(a_bin,I_bin,color=colours[col_i])#,alpha=0.5)
        # #---------------------------------------------------------------------------
        # # Find objects that lie in the cold classical region:
        # # binary objects
        # df_hel=pd.DataFrame(numpy.array([a_hel,e_hel,I_hel]).T,columns=['a_hel(AU)','e_hel','I_hel(d)'])
        # df_hel=df_hel[(df_hel['a_hel(AU)']>a_min) & (df_hel['a_hel(AU)']<a_max)]
        # df_hel=df_hel[(df_hel['I_hel(d)']>I_min) & (df_hel['I_hel(d)']<I_max)]
        # df_hel=df_hel[(df_hel['e_hel']<(1.0-(q/df_hel['a_hel(AU)'])))]
        # if len(df_hel)>0:
        #     print 'There are {} cold classical binaries formed'.format(len(df_hel))
        #     s1_2=ax1.scatter(df_hel['a_hel(AU)'],df_hel['e_hel'],c='r',marker='x')
        #     s2_2=ax2.scatter(df_hel['a_hel(AU)'],df_hel['I_hel(d)'],c='r',marker='x')
        # # repeat for single objects
        # df_hel_single=pd.DataFrame(numpy.array([a_hel_single,e_hel_single,I_hel_single]).T,columns=['a_hel(AU)','e_hel','I_hel(d)'])
        # df_hel_single=df_hel_single[(df_hel_single['a_hel(AU)']>a_min) & (df_hel_single['a_hel(AU)']<a_max)]
        # df_hel_single=df_hel_single[(df_hel_single['I_hel(d)']>I_min) & (df_hel_single['I_hel(d)']<I_max)]
        # df_hel_single=df_hel_single[(df_hel_single['e_hel']<(1.0-(q/df_hel_single['a_hel(AU)'])))]
        # if len(df_hel_single)>0:
        #     print 'There are {} cold classical singles formed'.format(len(df_hel_single))
        #     s1_3=ax1.scatter(df_hel_single['a_hel(AU)'],df_hel_single['e_hel'],c='r',marker='x')
        #     s2_3=ax2.scatter(df_hel_single['a_hel(AU)'],df_hel_single['I_hel(d)'],c='r',marker='x')
        #---------------------------------------------------------------------------
        # Unit conversions for plotting
        df_orb_all['a_hel(AU)']=numpy.array(df_orb_all['a_hel(m)'])/pf.AU
        df_orb_all['I_hel(d)']=numpy.array(df_orb_all['I_hel(rad)'])*(180.0/numpy.pi)
        a_hel_single=numpy.array(a_hel_single)/pf.AU
        I_hel_single=numpy.array(I_hel_single)*(180.0/numpy.pi)
        df_orb_all['I_bin(d)']=numpy.array(df_orb_all['I_bin(rad)'])*(180.0/numpy.pi)
        #---------------------------------------------------------------------------
        # Plot stuff for this timestep
        l1=ax5.axvline(t,0,30,color='k',alpha=0.5) # line to track simulation time

        df_single=pd.DataFrame(numpy.array([a_hel_single,e_hel_single,I_hel_single]).T,columns=['a_hel(AU)','e_hel','I_hel(d)'])
        print len(df_single)
        df_single_CC=CC_bounds(df_single)
        df_single=pd.concat([df_single,df_single_CC]).drop_duplicates(keep=False)
        print len(df_single),len(df_single_CC)

        s1_1=ax1.scatter(df_single['a_hel(AU)'],df_single['e_hel'],marker='^',facecolor=matplotlib.colors.colorConverter.to_rgba('black', alpha=0.5),edgecolor='None')
        s2_1=ax2.scatter(df_single['a_hel(AU)'],df_single['I_hel(d)'],marker='^',facecolor=matplotlib.colors.colorConverter.to_rgba('black', alpha=0.5),edgecolor='None')
        print 'There are {} cold classical singles formed'.format(len(df_single_CC))
        s1_3=ax1.scatter(df_single_CC['a_hel(AU)'],df_single_CC['e_hel'],marker='^',facecolor=matplotlib.colors.colorConverter.to_rgba('black', alpha=0.5),edgecolor='k')
        s2_3=ax2.scatter(df_single_CC['a_hel(AU)'],df_single_CC['I_hel(d)'],marker='^',facecolor=matplotlib.colors.colorConverter.to_rgba('black', alpha=0.5),edgecolor='k')

        # find the bin numbers for eccentricities
        e_bin_bin_num=pf.find_bin_num(numpy.array(df_orb_all['e_bin']),e_bins)
        df_orb_all['e_bin_bin_num']=e_bin_bin_num

        # tight binaries
        df_o=df_orb_all[df_orb_all['a/r_hill']<wide_cut]
        print len(df_o)
        df_o_CC=CC_bounds(df_o)
        df_o=pd.concat([df_o,df_o_CC]).drop_duplicates(keep=False)
        print len(df_o),len(df_o_CC)
        for j in range(len(df_o)):
            a_hel=float(df_o.iloc[j]['a_hel(AU)'])
            e_hel=float(df_o.iloc[j]['e_hel'])
            I_hel=float(df_o.iloc[j]['I_hel(d)'])
            a_bin=float(df_o.iloc[j]['a_bin(m)'])
            e_bin=float(df_o.iloc[j]['e_bin'])
            I_bin=float(df_o.iloc[j]['I_bin(d)'])
            p_num=int(df_o.iloc[j]['e_bin_bin_num'])+3

            s1_0=ax1.scatter(a_hel,e_hel,color=colours[col_i],s=s1,marker=(p_num,1,0))#,alpha=0.5)
            s2_0=ax2.scatter(a_hel,I_hel,color=colours[col_i],s=s1,marker=(p_num,1,0))#,alpha=0.5)
            s3=ax3.scatter(a_bin,e_bin,color=colours[col_i],s=s1,marker=(p_num,1,0))#,label="t={}Myr, a_neptune={}AU".format(t,float(a_neptune[numpy.where(time==t)])))#,alpha=0.5)
            s4=ax4.scatter(a_bin,I_bin,color=colours[col_i],s=s1,marker=(p_num,1,0))#,alpha=0.5)

        for j in range(len(df_o_CC)):
            a_hel=float(df_o_CC.iloc[j]['a_hel(AU)'])
            e_hel=float(df_o_CC.iloc[j]['e_hel'])
            I_hel=float(df_o_CC.iloc[j]['I_hel(d)'])
            a_bin=float(df_o_CC.iloc[j]['a_bin(m)'])
            e_bin=float(df_o_CC.iloc[j]['e_bin'])
            I_bin=float(df_o_CC.iloc[j]['I_bin(d)'])
            p_num=int(df_o_CC.iloc[j]['e_bin_bin_num'])+3

            # tight binaries in cold classical region
            s1_0=ax1.scatter(a_hel,e_hel,color=colours[col_i],s=s1,marker=(p_num,1,0),edgecolor='k')#,alpha=0.5)
            s2_0=ax2.scatter(a_hel,I_hel,color=colours[col_i],s=s1,marker=(p_num,1,0),edgecolor='k')#,alpha=0.5)
            s3=ax3.scatter(a_bin,e_bin,color=colours[col_i],s=s1,marker=(p_num,1,0),edgecolor='k')#,label="t={}Myr, a_neptune={}AU".format(t,float(a_neptune[numpy.where(time==t)])))#,alpha=0.5)
            s4=ax4.scatter(a_bin,I_bin,color=colours[col_i],s=s1,marker=(p_num,1,0),edgecolor='k')#,alpha=0.5)

        # wide binaries
        df_o=df_orb_all[df_orb_all['a/r_hill']>=wide_cut]
        print len(df_o)
        df_o_CC=CC_bounds(df_o)
        df_o=pd.concat([df_o,df_o_CC]).drop_duplicates(keep=False)
        print len(df_o),len(df_o_CC)
        for j in range(len(df_o)):
            a_hel=float(df_o.iloc[j]['a_hel(AU)'])
            e_hel=float(df_o.iloc[j]['e_hel'])
            I_hel=float(df_o.iloc[j]['I_hel(d)'])
            a_bin=float(df_o.iloc[j]['a_bin(m)'])
            e_bin=float(df_o.iloc[j]['e_bin'])
            I_bin=float(df_o.iloc[j]['I_bin(d)'])
            p_num=int(df_o.iloc[j]['e_bin_bin_num'])+3

            s1_0=ax1.scatter(a_hel,e_hel,color=colours[col_i],s=s2,marker=(p_num,1,0))#,alpha=0.5)
            s2_0=ax2.scatter(a_hel,I_hel,color=colours[col_i],s=s2,marker=(p_num,1,0))#,alpha=0.5)
            s3=ax3.scatter(a_bin,e_bin,color=colours[col_i],s=s2,marker=(p_num,1,0))#,label="t={}Myr, a_neptune={}AU".format(t,float(a_neptune[numpy.where(time==t)])))#,alpha=0.5)
            s4=ax4.scatter(a_bin,I_bin,color=colours[col_i],s=s2,marker=(p_num,1,0))#,alpha=0.5)

        for j in range(len(df_o_CC)):
            a_hel=float(df_o_CC.iloc[j]['a_hel(AU)'])
            e_hel=float(df_o_CC.iloc[j]['e_hel'])
            I_hel=float(df_o_CC.iloc[j]['I_hel(d)'])
            a_bin=float(df_o_CC.iloc[j]['a_bin(m)'])
            e_bin=float(df_o_CC.iloc[j]['e_bin'])
            I_bin=float(df_o_CC.iloc[j]['I_bin(d)'])
            p_num=int(df_o_CC.iloc[j]['e_bin_bin_num'])+3

            # wide binaries in cold classical region
            s1_0=ax1.scatter(a_hel,e_hel,color=colours[col_i],s=s2,marker=(p_num,1,0),edgecolor='k')#,alpha=0.5)
            s2_0=ax2.scatter(a_hel,I_hel,color=colours[col_i],s=s2,marker=(p_num,1,0),edgecolor='k')#,alpha=0.5)
            s3=ax3.scatter(a_bin,e_bin,color=colours[col_i],s=s2,marker=(p_num,1,0),edgecolor='k')#,label="t={}Myr, a_neptune={}AU".format(t,float(a_neptune[numpy.where(time==t)])))#,alpha=0.5)
            s4=ax4.scatter(a_bin,I_bin,color=colours[col_i],s=s2,marker=(p_num,1,0),edgecolor='k')#,alpha=0.5)

        #---------------------------------------------------------------------------
        # Find semimajor axis of Neptune
        t_mask=numpy.where(time==t)
        print "sim time = {} Myr, a_neptune = {} AU".format(time[t_mask],a_planet[:,N_planet-1][t_mask])
        # ---------------------------------------------------------------------------


        col_i+=1

ax3.legend()

# # ELLIPSES
# ylims=ax1.get_ylim()
# xlims=ax1.get_xlim()
# print xlims,ylims
#
# def ax_scale(val,lims):
#     # return ((lims[1]-lims[0])*val)+lims[0]
#     return (val-lims[0])/(lims[1]-lims[0])
# def ax_scale_inv(val,lims):
#     return ((lims[1]-lims[0])*val)+lims[0]
#     # return (val-lims[0])/(lims[1]-lims[0])
#
# print ax_scale(numpy.array(xlims),numpy.array(xlims))
# print ax_scale(numpy.array(ylims),numpy.array(ylims))
#
# print ax1.transData.transform((xlims[0], ylims[0]))
# print ax1.transData.transform((xlims[1], ylims[1]))
# ax1.set_xlim(xlims[0],xlims[1])
# ax1.set_ylim(ylims[0],ylims[1])
# size=10
# # inv = ax1.transData.inverted()
# zero=ax1.transData.transform((xlims[0], ylims[0]))
# print "zero={}".format(zero)
# zero=ax1.transData.inverted().transform(zero)
# ax1.scatter(zero[0],zero[1],marker='+',color='k')
#
# # use the scale of display to set x and y axis scaling
# print (xlims[0], ylims[0]),(xlims[1], ylims[1])
# zero_0=ax1.transData.transform((xlims[0], ylims[0]))
# zero_1=ax1.transData.transform((xlims[1], ylims[1]))
# print "zeros",zero_0,zero_1
# x_disp_scale=zero_1[0]-zero_0[0]
# y_disp_scale=zero_1[1]-zero_0[1]
# sizex_disp=20
# for j in range(len(df_orb_all)):
#     a=float(df_orb_all.iloc[j]['a_hel(AU)'])
#     e=float(df_orb_all.iloc[j]['e_hel'])
#     ae=(a,e)
#
#     sizey_disp=sizex_disp*numpy.sqrt(1.0-float(df_orb_all.iloc[j]['e_bin'])**2.0)
#     # sizey_disp=sizex_disp*0.5
#     # print float(df_orb_all.iloc[j]['e_bin']),sizex_disp,sizey_disp
#     sizex=((xlims[1]-xlims[0])/(x_disp_scale))*sizex_disp
#     sizey=((ylims[1]-ylims[0])/(y_disp_scale))*sizey_disp
#
#     wh=[sizex,sizey]
#     print ae,wh
#     ax1.scatter(ae[0],ae[1],marker='x',color='k')
#     a=45#numpy.random.rand()*360
#     # These work except when setting an angle
#     # el=Ellipse(xy=ae, width=wh[0]*numpy.cos(a*(numpy.pi/180.0)), height=wh[1]*numpy.sin(a*(numpy.pi/180.0)), angle=a, facecolor="none", edgecolor='r')#,transform=ax1.transData)
#     el=Ellipse(xy=ae, width=wh[0], height=wh[1], facecolor="none", edgecolor='r')#, angle=a)
#
#     # ae_disp=ax1.transData.transform(ae)
#     # wh_disp=(sizex_disp,sizey_disp)
#     # print ae_disp,wh_disp
#     # el=Ellipse(xy=ae_disp, width=sizex_disp, height=sizey_disp,facecolor="none", edgecolor='r',transform=None)
#
#     # t2 = matplotlib.transforms.Affine2D().rotate_deg(-45) + ax1.transData
#     # el.set_transform(el)
#
#     # # Try rotate the points
#     # ax1.scatter(ae[0]+wh[0]/2.0,ae[1],marker='x',color='k')
#     # ax1.scatter(ae[0]-wh[0]/2.0,ae[1],marker='x',color='k')
#     #
#     # a=a*(numpy.pi/180.0)
#     # wh_disp=ax1.transData.transform(wh)
#     # wh_disp[0]=wh_disp[0]*numpy.cos(a)
#     # wh_disp[1]=wh_disp[1]*numpy.sin(a)
#     # wh=ax1.transData.inverted().transform(wh_disp)
#     #
#     # ax1.scatter(ae[0]+(wh[0]/2.0),ae[1]+(wh[1]/2.0),marker='x',color='k')
#     # ax1.scatter(ae[0]-(wh[0]/2.0),ae[1]-(wh[1]/2.0),marker='x',color='k')
#
#     # ax1.scatter(ae[0]+(wh[0]/2.0),ae[1]+(wh[0]/2.0)*numpy.sin(a),marker='x',color='k')
#     # ax1.scatter(ae[0]-(wh[0]/2.0)*numpy.cos(a),ae[1]-(wh[0]/2.0)*numpy.sin(a),marker='x',color='k')
#
#     ax1.add_patch(el)
#
# # # Ellipse pos okay, can't get shape?
# # for j in range(len(df_orb_all)):
# #     a=float(df_orb_all.iloc[j]['a_hel(AU)'])
# #     e=float(df_orb_all.iloc[j]['e_hel'])
# #     ae=ax1.transData.transform((a, e))
# #     wh=[size,size*numpy.sqrt(1-float(df_orb_all.iloc[j]['e_bin'])**2.0)]
# #     print ae,wh
# #     ae=ax1.transData.inverted().transform(ae)
# #     ax1.scatter(ae[0],ae[1],marker='x',color='k')
# #     el=Ellipse(xy=ae, width=1, height=0.05, color='r')
# #     ax1.add_artist(el)
#
# # # Ellipse shape correct but offset?
# # for j in range(len(df_orb_all)):
# #     a=float(df_orb_all.iloc[j]['a_hel(AU)'])
# #     e=float(df_orb_all.iloc[j]['e_hel'])
# #     ae=ax1.transData.transform((a, e))
# #     wh=[size,size*numpy.sqrt(1-float(df_orb_all.iloc[j]['e_bin'])**2.0)]
# #     print ae,wh
# #     el=Ellipse(xy=(ae[0],ae[1]), width=wh[0], height=wh[1], angle=numpy.random.rand()*360, color='r',transform=None)
# #     ae=ax1.transData.inverted().transform(ae)
# #     ax1.scatter(ae[0],ae[1],marker='x',color='k')
# #     ax1.add_artist(el)

# save and clear the figure for the next timestep
save_path = '../analysis_dirs/{}/{}_figs/{}_{}_{:07d}.png'.format(run,os.path.basename(__file__).split('.')[0],os.path.basename(__file__).split('.')[0],run,int(fi.split('/')[-1].split('.')[0]))
print save_path
pyplot.savefig(save_path, bbox_inches='tight')
pyplot.show()
# pyplot.close()

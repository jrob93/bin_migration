import numpy
import matplotlib.pyplot as pyplot
import os
from natsort import natsorted, ns
import subprocess
import matplotlib.gridspec as gridspec
import sys
sys.path.insert(0, '../../../grav_cloud/python_stuff/')
import py_func as pf
import pandas as pd

N_planet=4

# load planet data
df_planet=pd.read_csv("../gen_planet_data.txt",sep="\t",index_col=0)
print df_planet

path='/Users/jrobinson/bin_migration/py_stuff'

df_fname='../df_run_params.txt'
df_run_params=pd.read_csv(df_fname,sep="\t")

# Find all Unique values
df_tau=df_run_params.drop_duplicates(['tau_a(Myr)','tau_e(Myr)','tau_i(Myr)','tmax(Myr)'])
print df_tau
print len(df_tau)

# Find all runs with the same parameters
i=0
for i in range(len(df_tau)):
    vals=numpy.array(df_tau[['tau_a(Myr)','tau_e(Myr)','tau_i(Myr)','tmax(Myr)']].iloc[i])
    dirs=[]
    for j in range(len(df_run_params)):
        if (df_run_params['tau_a(Myr)'].iloc[j]==vals[0]) and (df_run_params['tau_e(Myr)'].iloc[j]==vals[1]) and (df_run_params['tau_i(Myr)'].iloc[j]==vals[2]) and (df_run_params['tmax(Myr)'].iloc[j]==vals[3]):
            dirs.append(df_run_params.iloc[j]['run'])

    # dirs=['prone_120_analysis', 'prone_121_analysis']
    print vals
    print dirs
    # exit()

    # for run in dirs:
    run=dirs[0]
    # subprocess.Popen(["mkdir","../{}_analysis".format(run)]) # make directory to store images

    # Load the planet data
    for i in range(1,N_planet+1):
        _i=i-1
        # planet_record_file="../{}_analysis/planet_{:03d}.txt".format(run,i)
        planet_record_file="../{}/planet_{:03d}.txt".format(run,i)
        filename='{}'.format(planet_record_file)
        print filename
        try:
            with open(filename) as f:
                dat = [line.rstrip() for line in f]
            dat = [line.split() for line in dat]
            dat = numpy.array(dat,dtype=float) # convert data list into array
        except:
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

    print numpy.shape(a_planet)
    print numpy.shape(time)

    fig = pyplot.figure()
    fig.set_size_inches(15, 10)
    gs = gridspec.GridSpec(3, 1)
    ax1 = pyplot.subplot(gs[0,0])
    ax2 = pyplot.subplot(gs[1,0],sharex=ax1)
    ax3 = pyplot.subplot(gs[2,0],sharex=ax1)

    ax3.set_xlabel('t (Myr)')
    ax1.set_ylabel('a (AU)')
    ax2.set_ylabel('e')
    # ax3.set_ylabel('I (rad)')
    ax3.set_ylabel('I (degrees)')

    # /* The options we understand. */
    # static struct argp_option options[] = {
    #   {"nep_a_i",  'a', "nepai",      0,  "Initial Neptune Semi (AU)"},
    #   {"nep_a_f",    'A', "nepaf",      0,  "Final Neptune Semi (AU)" },
    #   {"tau_a",   't', "taua", 0,      "Migration Timescale (yr)" },
    #   {"tau_e",   'E', "taue", 0,      "Eccentricity damping timescale (yr)" },
    #   {"tau_i",   'I', "taui", 0,      "Inclination damping timescale (yr)" },
    #   {"nepeo", 'e', "nepeo", 0, "Neptune's original eccentricity."},
    #   {"nepio", 'i', "nepio", 0, "Neptune's original inclincation (deg)."},
    #   {"disk", 'd', "IncDisk", 0, "Include a binary particle small disk? Provided integer is the number of binaries in the disk."},
    #   {"tmax", 'm', "tmax", 0, "Maximum simulation time (yr)"},
    #   {"jump", 'j', "jump", 0, "Run Jump simulation? Provided integer is number of single particles to read. # of binaries comes from -d integer."},
    #   {"timeStep", 's', "timeStep", 0.5, "The integration timestep in days."},
    #   { 0 }
    # };
    print run
    params=numpy.array(header.split()[6:]).astype(float)
    print params
    # fig.suptitle('run: {}, {}'.format(run," ".join(header.split()[4:])))
    fig.suptitle('N_a_i={}AU, N_a_f={}AU, tau_a={}Myr, tau_e={}Myr, tau_i={}Myr, N_e_o={}, N_e_i={}deg, tmax={}Myr, dt={}d'.format(
    params[0],params[1],params[2]/1e6,params[3]/1e6,params[4]/1e6,params[5],params[6],params[7]/1e6,params[8]))

    # break

    planet_name=["Jupiter","Saturn","Uranus","Neptune"]
    for j in range(N_planet):
        ax1.plot(time,a_planet[:,j],label=planet_name[j],zorder=0)
        ax2.plot(time,e_planet[:,j],zorder=0)
        # ax3.plot(time,I_planet[:,j])
        ax3.plot(time,numpy.array(I_planet[:,j])/(numpy.pi/180.0),zorder=0)

    ax1.legend()

    colours=pyplot.rcParams['axes.prop_cycle'].by_key()['color']

    for j in range(N_planet):
        # also plot q and Q
        q=a_planet[:,j]*(1.0-e_planet[:,j])
        Q=a_planet[:,j]*(1.0+e_planet[:,j])
        ax1.plot(time,q,linestyle=':',color=colours[j],alpha=0.3,zorder=0)
        ax1.plot(time,Q,linestyle=':',color=colours[j],alpha=0.3,zorder=0)

    # plot planet data
    for k in range(N_planet):
        k=k+5
        ax1.scatter(time[-1],df_planet['a(m)'].iloc[k]/pf.AU,zorder=2,edgecolor='k')
        ax2.scatter(time[-1],df_planet['e'].iloc[k],zorder=2,edgecolor='k')
        ax3.scatter(time[-1],df_planet['I(rad)'].iloc[k]/(numpy.pi/180.0),zorder=2,edgecolor='k')

    # save_path = '../{}_analysis/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
    save_path = './{}_{:d}_{:d}_{:d}_{:d}.png'.format(os.path.basename(__file__).split('.')[0],
    int(vals[0]),int(vals[1]),int(vals[2]),int(vals[3]))
    print save_path
    pyplot.savefig(save_path, bbox_inches='tight')
    pyplot.show()
    # pyplot.close()
    # break

import numpy
import matplotlib.pyplot as pyplot
import os
from natsort import natsorted, ns
import subprocess
import matplotlib.gridspec as gridspec
import sys
sys.path.insert(0, '../../../grav_cloud/python_stuff/')
# import py_func as pf
import pandas as pd

path='/Users/jrobinson/bin_migration/py_stuff'
dirs=next(os.walk(path))[1]
dirs = [d for d in dirs if 'prone' in d]
dirs=natsorted(dirs, key=lambda y: y.lower()) # sort unpadded numbers

print dirs

df_fname='df_run_params.txt'

df=open(df_fname,'w')
df.write('run\tN_a_i(AU)\tN_a_f(AU)\ttau_a(Myr)\ttau_e(Myr)\ttau_i(Myr)\tN_e_o\tN_i_o(deg)\ttmax(Myr)\tdt(d)\n')

for run in dirs:

    with open("./{}/0.dump".format(run)) as f:
        _dat = [line.rstrip() for line in f]
    header = _dat[0] # header line
    # print header

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

    # print run
    params=numpy.array(header.split()[6:]).astype(float)
    # print params
    df.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(run,params[0],params[1],params[2]/1e6,params[3]/1e6,params[4]/1e6,params[5],params[6],params[7]/1e6,params[8]))

df.close()

df_run_params=pd.read_csv(df_fname,sep="\t")
print df_run_params
print len(df_run_params)

# Find all Unique values
# df_tau=df_run_params[['tau_a(Myr)','tau_e(Myr)','tau_i(Myr)']]
df_tau=df_run_params.drop_duplicates(['tau_a(Myr)','tau_e(Myr)','tau_i(Myr)','tmax(Myr)'])
print df_tau#[['run','tau_a(Myr)','tau_e(Myr)','tau_i(Myr)','tmax(Myr)']]

# Find all runs with the same parameters
for i in range(len(df_tau)):
    vals=numpy.array(df_tau[['tau_a(Myr)','tau_e(Myr)','tau_i(Myr)','tmax(Myr)']].iloc[i])
    print vals
    runs=[]
    for j in range(len(df_run_params)):
        if (df_run_params['tau_a(Myr)'].iloc[j]==vals[0]) and (df_run_params['tau_e(Myr)'].iloc[j]==vals[1]) and (df_run_params['tau_i(Myr)'].iloc[j]==vals[2]) and (df_run_params['tmax(Myr)'].iloc[j]==vals[3]):
            runs.append(df_run_params.iloc[j]['run'])
    print runs
    print len(runs)

import numpy
import matplotlib.pyplot as pyplot
import os
from natsort import natsorted, ns
import subprocess
import matplotlib.gridspec as gridspec
import sys
sys.path.insert(0, '../../../grav_cloud/python_stuff/')
import py_func as pf

path='/Users/jrobinson/bin_migration/'
run='200_dump'

# path='/Users/jrobinson/bin_migration/migration'
# run='binSim'
# # run='binSim_20_0'

lim=50.0
n=100
N_planet=4

subprocess.Popen(["mkdir","../{}_analysis".format(run)]) # make directory to store images

run_path='{}/{}'.format(path,run)
files=next(os.walk(run_path))[2]
files = [ fi for fi in files if fi.endswith('.dump')]
files=natsorted(files, key=lambda y: y.lower()) # sort unpadded numbers

# time=[]
# a_planet=numpy.zeros((len(files),N_planet))
# e_planet=numpy.zeros((len(files),N_planet))
# I_planet=numpy.zeros((len(files),N_planet))
# for i in range(len(files)):
#     filename='{}/{}'.format(run_path,files[i])
#     print filename
#     with open(filename) as f:
#         dat = [line.rstrip() for line in f]
#     dat = [line.split() for line in dat]
#     print dat[0] # header line
#     t=float(dat[0][3]) # timestamp
#     dat = dat[1:] # strip header line
#     dat = numpy.array(dat,dtype=float) # convert data list into array
#     time.append(t)
#     # data format:
#     # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0
#     a_planet[i,:]=dat[1:N_planet+1,1]
#     e_planet[i,:]=dat[1:N_planet+1,2]
#     I_planet[i,:]=dat[1:N_planet+1,3]

# LOWER RESOLUTION
time=[]
file_nums=range(0,len(files),1000)
print file_nums
a_planet=numpy.zeros((len(file_nums),N_planet))
e_planet=numpy.zeros((len(file_nums),N_planet))
I_planet=numpy.zeros((len(file_nums),N_planet))
j=0
for i in file_nums:
    filename='{}/{}'.format(run_path,files[i])
    print filename
    with open(filename) as f:
        dat = [line.rstrip() for line in f]
    dat = [line.split() for line in dat]
    print dat[0] # header line
    t=float(dat[0][3]) # timestamp
    dat = dat[1:] # strip header line
    dat = numpy.array(dat,dtype=float) # convert data list into array
    time.append(t)
    # data format:
    # ii,a,e,inc,omega,Omega,f,x-x0,y-y0,z-z0,vx-vx0,vy-vy0,vz-vz0
    a_planet[j,:]=dat[1:N_planet+1,1]
    e_planet[j,:]=dat[1:N_planet+1,2]
    I_planet[j,:]=dat[1:N_planet+1,3]
    j+=1
fig = pyplot.figure()
gs = gridspec.GridSpec(3, 1)
ax1 = pyplot.subplot(gs[0,0])
ax2 = pyplot.subplot(gs[1,0],sharex=ax1)
ax3 = pyplot.subplot(gs[2,0],sharex=ax1)

ax3.set_xlabel('t (Myr)')
ax1.set_ylabel('a (AU)')
ax2.set_ylabel('e')
ax3.set_ylabel('I (rad)')
fig.suptitle('run: {}'.format(run))


for j in range(N_planet):
    ax1.plot(time,a_planet[:,j])
    ax2.plot(time,e_planet[:,j])
    ax3.plot(time,I_planet[:,j])

colours=pyplot.rcParams['axes.prop_cycle'].by_key()['color']

for j in range(N_planet):

    # also plot q and Q
    q=a_planet[:,j]*(1.0-e_planet[:,j])
    Q=a_planet[:,j]*(1.0+e_planet[:,j])
    ax1.plot(time,q,linestyle=':',color=colours[j],alpha=0.3)
    ax1.plot(time,Q,linestyle=':',color=colours[j],alpha=0.3)

save_path = '../{}_analysis/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
print save_path
pyplot.savefig(save_path, bbox_inches='tight')
pyplot.show()
# pyplot.close()

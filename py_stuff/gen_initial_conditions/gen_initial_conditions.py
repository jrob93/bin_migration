import numpy
import sys
sys.path.insert(0, '../../../grav_cloud/python_stuff/')
import py_func as pf
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import rebound

pi=numpy.pi

N_bin=500

# binary orbit, min and max values
mi=0.5e18
mj=mi
a_bin=numpy.array([5e3,5e4])
e_bin=numpy.array([0.0,1.0])
I_bin=numpy.array([-90.0,90.0])*(pi/180.0)
o_bin=numpy.array([0.0,360.0])*(pi/180.0)
O_bin=numpy.array([0.0,360.0])*(pi/180.0)
f_bin=numpy.array([0.0,360.0])*(pi/180.0)

# helio orbit, min and max values
a_hel=numpy.array([30.0,35.0])*pf.AU
e_hel=numpy.array([0.0,0.05])
I_hel=numpy.array([-5.0,5.0])*(pi/180.0)
o_hel=numpy.array([0.0,360.0])*(pi/180.0)
O_hel=numpy.array([0.0,360.0])*(pi/180.0)
f_hel=numpy.array([0.0,360.0])*(pi/180.0)

pos=[]
vel=[]
m=[]

of=open('binary_initial.txt','w')

for i in range(N_bin):
    # place primary
    mu_hel=pf.G*pf.M_sun
    orb_hel=numpy.zeros(6)
    orb_hel[0]=numpy.random.uniform(a_hel[0],a_hel[1])
    orb_hel[1]=numpy.random.uniform(e_hel[0],e_hel[1])
    orb_hel[2]=numpy.random.uniform(I_hel[0],I_hel[1])
    orb_hel[3]=numpy.random.uniform(o_hel[0],o_hel[1])
    orb_hel[4]=numpy.random.uniform(O_hel[0],O_hel[1])
    orb_hel[5]=numpy.random.uniform(f_hel[0],f_hel[1])

    pi,vi=pf.pos_vel_from_orbit(orb_hel,mu_hel)
    pos.append(pi)
    vel.append(vi)
    m.append(mi)

    # place secondary
    mu_bin=pf.G*(mi+mj)
    orb_bin=numpy.zeros(6)
    orb_bin[0]=numpy.random.uniform(a_bin[0],a_bin[1])
    orb_bin[1]=numpy.random.uniform(e_bin[0],e_bin[1])
    orb_bin[2]=numpy.random.uniform(I_bin[0],I_bin[1])
    orb_bin[3]=numpy.random.uniform(o_bin[0],o_bin[1])
    orb_bin[4]=numpy.random.uniform(O_bin[0],O_bin[1])
    orb_bin[5]=numpy.random.uniform(f_bin[0],f_bin[1])

    print "orbit ",orb_bin

    pj,vj=pf.pos_vel_from_orbit(orb_bin,mu_bin)
    pj=pi+pj
    vj=vi+vj
    pos.append(pj)
    vel.append(vj)
    m.append(mj)

    of.write("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n".format(pi[0],pi[1],pi[2],vi[0],vi[1],vi[2],mi))
    of.write("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n".format(pj[0],pj[1],pj[2],vj[0],vj[1],vj[2],mj))

    print pi,vi
    print pj,vj

of.close()

pos=numpy.array(pos)
vel=numpy.array(vel)
m=numpy.array(m)

print pos.shape
# print pos

# fig = pyplot.figure()
# gs = gridspec.GridSpec(1, 1)
# ax1 = pyplot.subplot(gs[0,0])
# ax1.set_aspect("equal")
# ax1.set_xlabel('x')
# ax1.set_ylabel('y')
# # ax1.set_xlim(-lim,lim)
# # ax1.set_ylim(-lim,lim)
#
# ax1.scatter(pos[:,0],pos[:,1],alpha=0.2) # plot small bodies
#
# # save_path = '{}/{}.png'.format(run,files[i])
# # print save_path
# # pyplot.savefig(save_path, bbox_inches='tight')
# pyplot.show()
# # pyplot.close()

# Calculate all binary pairs and store a, e, I
a_bin=[]
e_bin=[]
I_bin=[]

i=0
while i < len(pos)-1:
    j=i+1

    print pos[i,:],vel[i,:]
    print pos[j,:],vel[j,:]

    # use rebound to calculate the orbit, better handling for certain cases, e.g. circular orbits
    sim = rebound.Simulation()
    sim.G=pf.G
    sim.add(x=pos[i,0],y=pos[i,1],z=pos[i,2],
    vx=vel[i,0],vy=vel[i,1],vz=vel[i,2],
    m=m[i])
    sim.add(x=pos[j,0],y=pos[j,1],z=pos[j,2],
    vx=vel[j,0],vy=vel[j,1],vz=vel[j,2],
    m=m[j])
    orbit = sim.particles[1].calculate_orbit(sim.particles[0])
    print orbit

    a_bin.append(orbit.a)
    e_bin.append(orbit.e)
    I_bin.append(orbit.inc)

    i+=2

fig = pyplot.figure()
gs = gridspec.GridSpec(3, 1)
ax1 = pyplot.subplot(gs[0,0])
ax2 = pyplot.subplot(gs[1,0])
ax3 = pyplot.subplot(gs[2,0])
ax1.set_xlabel('a (m)')
ax2.set_xlabel('e')
ax3.set_xlabel('I')
ax1.set_ylabel('N')
ax2.set_ylabel('N')
ax3.set_ylabel('N')

# fig.suptitle('run: {}'.format(run))

n_bins=100
ax1.hist(a_bin,bins=n_bins)
ax2.hist(e_bin,bins=n_bins)
ax3.hist(I_bin,bins=n_bins)

pyplot.tight_layout()

# save_path = '{}/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
# print save_path
# pyplot.savefig(save_path, bbox_inches='tight')
pyplot.show()
# pyplot.close()

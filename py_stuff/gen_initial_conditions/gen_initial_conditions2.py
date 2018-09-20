import numpy
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import rebound

pi=numpy.pi
cos = numpy.cos
sin = numpy.sin
G=6.67428e-11 # m3 kg-1 s-2
M_sun=1.98855e30 # kg
AU=1.496e11# m

#-------------------------------------------------------------------------------
def r_elliptical(a,e,theta):
    '''
    Function to find distance, r, at true anomaly, theta, around a keplerian orbit

    Parameters
    ----------
    a
        semi major axis (m)
    e
        orbital eccentricity
    theta
        true anomaly
    '''
    r = (a*(1-e**2))/(1+e*cos(theta))
    return r
#-------------------------------------------------------------------------------
def find_e_p(OMEGA,omega,I):
    '''Function to find the normalised Lenz vector along the apsidal line, uses PSS eq11.39a'''
    e_p=numpy.zeros(3)
    e_p[0]=(cos(OMEGA)*cos(omega))-(cos(I)*sin(OMEGA)*sin(omega)) #x
    e_p[1]=(sin(OMEGA)*cos(omega))+(cos(I)*cos(OMEGA)*sin(omega)) #y
    e_p[2]=sin(I)*sin(omega) #z
    return e_p
#-------------------------------------------------------------------------------
def find_e_Q(OMEGA,omega,I):
    '''Function to find the normalised e_Q vector which is perp to h and e_a and
    lies in the orbital plane, uses PSS eq11.39b'''
    e_Q=numpy.zeros(3)
    e_Q[0]=-(cos(OMEGA)*sin(omega))-(cos(I)*sin(OMEGA)*cos(omega)) #x
    e_Q[1]=-(sin(OMEGA)*sin(omega))+(cos(I)*cos(OMEGA)*cos(omega)) #y
    e_Q[2]=sin(I)*cos(omega) #z
    return e_Q
#-------------------------------------------------------------------------------
def pos_vel_from_orbit(orb,mu):
    '''function to find xyz and vxvyvz from the particle orbit

    Parameters
    ----------
    orb
        list of orbital elements [a(m),e,I(rad),OMEGA(rad),omega(rad),f_true(rad)]
    mu
        mu=G(M+m)
    '''
    a=orb[0]
    e=orb[1]
    I=orb[2]
    OMEGA=orb[3]
    omega=orb[4]
    f_true = orb[5]
    n=numpy.sqrt(mu/(a**3))
    eta=numpy.sqrt(1.0-(e**2))
    r=r_elliptical(a,e,f_true)
    #-------------------------------------------------------------------------------
    #determine orbital basis vectors
    e_p=find_e_p(OMEGA,omega,I)
    e_Q=find_e_Q(OMEGA,omega,I)
    #-------------------------------------------------------------------------------
    #define the r(x,y,z) position array
    pos=numpy.zeros(3)
    pos=r*((cos(f_true)*e_p)+(sin(f_true)*e_Q))#PSS eq 11.36a
    vel=numpy.zeros(3)
    vel=(n*a/eta)*((-sin(f_true)*e_p)+((e+cos(f_true))*e_Q))
    return pos,vel
#-------------------------------------------------------------------------------

N_bin=500

# binary orbit, min and max values
mi=0.5e18
mj=mi
a_bin=numpy.array([5e6,5e7])
e_bin=numpy.array([0.0,1.0])
#I_bin=numpy.array([-90.0,90.0])*(pi/180.0)
I_bin=numpy.array([0.0,180.0])*(pi/180.0)
o_bin=numpy.array([0.0,360.0])*(pi/180.0)
O_bin=numpy.array([0.0,360.0])*(pi/180.0)
f_bin=numpy.array([0.0,360.0])*(pi/180.0)

# helio orbit, min and max values
a_hel=numpy.array([30.0,35.0])*AU
e_hel=numpy.array([0.0,0.05])
I_hel=numpy.array([0.0,5.0])*(pi/180.0)
o_hel=numpy.array([0.0,360.0])*(pi/180.0)
O_hel=numpy.array([0.0,360.0])*(pi/180.0)
f_hel=numpy.array([0.0,360.0])*(pi/180.0)

pos=[]
vel=[]
m=[]

of=open('binary_initial.txt','w')

for i in range(N_bin):
    # place primary
    mu_hel=G*M_sun
    orb_hel=numpy.zeros(6)
    orb_hel[0]=numpy.random.uniform(a_hel[0],a_hel[1])
    orb_hel[1]=numpy.random.uniform(e_hel[0],e_hel[1])
    orb_hel[2]=numpy.random.uniform(I_hel[0],I_hel[1])
    orb_hel[3]=numpy.random.uniform(o_hel[0],o_hel[1])
    orb_hel[4]=numpy.random.uniform(O_hel[0],O_hel[1])
    orb_hel[5]=numpy.random.uniform(f_hel[0],f_hel[1])

    pi,vi=pos_vel_from_orbit(orb_hel,mu_hel)
    pos.append(pi)
    vel.append(vi)
    m.append(mi)

    # place secondary
    mu_bin=G*(mi+mj)
    orb_bin=numpy.zeros(6)
    orb_bin[0]=numpy.random.uniform(a_bin[0],a_bin[1])
    orb_bin[1]=numpy.random.uniform(e_bin[0],e_bin[1])
    orb_bin[2]=numpy.random.uniform(I_bin[0],I_bin[1])
    orb_bin[3]=numpy.random.uniform(o_bin[0],o_bin[1])
    orb_bin[4]=numpy.random.uniform(O_bin[0],O_bin[1])
    orb_bin[5]=numpy.random.uniform(f_bin[0],f_bin[1])

    print "orbit ",orb_bin
    print "binary inclination: ",i,orb_bin[2]/(numpy.pi/180.0)

    pj,vj=pos_vel_from_orbit(orb_bin,mu_bin)
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
    sim.G=G
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

    print "binary inclination: ",orbit.inc/(numpy.pi/180.0)

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
ax3.hist(numpy.array(I_bin)/(numpy.pi/180.0),bins=n_bins)

pyplot.tight_layout()

# save_path = '{}/{}_{}.png'.format(run,os.path.basename(__file__).split('.')[0],run)
# print save_path
# pyplot.savefig(save_path, bbox_inches='tight')
pyplot.show()
# pyplot.close()

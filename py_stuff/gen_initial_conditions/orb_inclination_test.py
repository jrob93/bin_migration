import numpy
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import rebound
from mpl_toolkits.mplot3d import Axes3D

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

N=100

mu=G*M_sun

I1=135.0*(pi/180.0)
I2=-45.0*(pi/180.0)
I3=45.0*(pi/180.0)

f_true=numpy.linspace(0.0,2.0*pi,N)

# define empty arrays for position and velocity
pos1=numpy.zeros((N,3))
pos2=numpy.zeros((N,3))
pos3=numpy.zeros((N,3))
vel1=numpy.zeros((N,3))
vel2=numpy.zeros((N,3))
vel3=numpy.zeros((N,3))

# find pos and vel for all values of f_true
for i in range(N):
    orb1=[1.0*AU,0.4,I1,0.0,0.0,f_true[i]]
    orb2=[1.0*AU,0.4,I2,0.0,0.0,f_true[i]]
    orb3=[1.0*AU,0.4,I3,0.0,0.0,f_true[i]]

    pos1[i,:],vel1[i,:]=pos_vel_from_orbit(orb1,mu)
    pos2[i,:],vel2[i,:]=pos_vel_from_orbit(orb2,mu)
    pos3[i,:],vel3[i,:]=pos_vel_from_orbit(orb3,mu)

print pos1[0,:],vel1[0,:]
print pos2[0,:],vel2[0,:]
# print pos3[0,:],vel3[0,:]

# Now double check the orbit, calculate the orbit from position and velocity at f_true=0
sim = rebound.Simulation()
sim.G=G
sim.add(x=0.0,y=0.0,z=0.0,
vx=0.0,vy=0.0,vz=0.0,
m=M_sun)
sim.add(x=pos1[0,0],y=pos1[0,1],z=pos1[0,2],
vx=vel1[0,0],vy=vel1[0,1],vz=vel1[0,2],
m=1.0)
sim.add(x=pos2[0,0],y=pos2[0,1],z=pos2[0,2],
vx=vel2[0,0],vy=vel2[0,1],vz=vel2[0,2],
m=1.0)
sim.add(x=pos3[0,0],y=pos3[0,1],z=pos3[0,2],
vx=vel3[0,0],vy=vel3[0,1],vz=vel3[0,2],
m=1.0)

orbit1 = sim.particles[1].calculate_orbit(sim.particles[0])
print orbit1
print 'inclination = {} degrees'.format(orbit1.inc/(pi/180.0))
orbit2 = sim.particles[2].calculate_orbit(sim.particles[0])
print orbit2
print 'inclination = {} degrees'.format(orbit2.inc/(pi/180.0))
orbit3 = sim.particles[3].calculate_orbit(sim.particles[0])
print orbit3
print 'inclination = {} degrees'.format(orbit3.inc/(pi/180.0))


fig = pyplot.figure()

ax1 = fig.add_subplot(111)
ax1 = fig.add_subplot(111,projection='3d')
ax1.set_aspect("equal")
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_zlabel('z (m)')

colours=pyplot.rcParams['axes.prop_cycle'].by_key()['color']

# optional: shift along the y axis to better display orbit 2 for comparison
pos2[:,1]+=AU

# plot the orbit in 3d space
ax1.plot(pos1[:,0],pos1[:,1],pos1[:,2],c=colours[0],zorder=0)
ax1.plot(pos2[:,0],pos2[:,1],pos2[:,2],c=colours[1],zorder=0)
ax1.plot(pos3[:,0],pos3[:,1],pos3[:,2],c=colours[2],zorder=0)

# plot the position in orbit at f_true=0
ax1.scatter(pos1[0,0],pos1[0,1],pos1[0,2],c=colours[0],edgecolor='k',zorder=2)
ax1.scatter(pos2[0,0],pos2[0,1],pos2[0,2],c=colours[1],edgecolor='k',zorder=2)
ax1.scatter(pos3[0,0],pos3[0,1],pos3[0,2],c=colours[2],edgecolor='k',zorder=2)

# Plot the velocity vector at f_true=0, note that we scale the vector length to be visible on the figure
vel1=vel1*1e6
vel2=vel2*1e6
vel3=vel3*1e6

ax1.plot([pos1[0,0],pos1[0,0]+vel1[0,0]],[pos1[0,1],pos1[0,1]+vel1[0,1]],[pos1[0,2],pos1[0,2]+vel1[0,2]],c='k',zorder=1)
ax1.plot([pos2[0,0],pos2[0,0]+vel2[0,0]],[pos2[0,1],pos2[0,1]+vel2[0,1]],[pos2[0,2],pos2[0,2]+vel2[0,2]],c='k',zorder=1)
ax1.plot([pos3[0,0],pos3[0,0]+vel3[0,0]],[pos3[0,1],pos3[0,1]+vel3[0,1]],[pos3[0,2],pos3[0,2]+vel3[0,2]],c='k',zorder=1)

pyplot.show()

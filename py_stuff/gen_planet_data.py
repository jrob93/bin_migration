import numpy
import rebound
import os

# Add planets to plot current solar system positions

# planet_list=["Sun","Jupiter"]#,"Saturn","Uranus","Neptune"]
planet_list=["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"]

sim = rebound.Simulation()
sim.G=6.67428e-11
sim.units = ('m', 's', 'kg')
for name in planet_list:
    sim.add(name)
sim.status()

orbits=sim.calculate_orbits(primary=sim.particles[0])

f=open("{}.txt".format(os.path.basename(__file__).split('.')[0]),'w')
f.write("index\ta(m)\te\tI(rad)\tOMEGA(rad)\tomega(rad)\tf_true(rad)\tm(kg)\tx(m)\ty(m)\tz(m)\tvx(ms^-1)\tvy(ms^-1)\tvz(ms^-1)\n")
for i in range(len(planet_list)):
    p=sim.particles[i]
    if i==0:
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(i,'nan','nan','nan','nan','nan','nan',p.m,p.x,p.y,p.z,p.vx,p.vy,p.vz))
    else:
        o=orbits[i-1]
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(i,o.a,o.e,o.inc,o.Omega,o.omega,o.f,p.m,p.x,p.y,p.z,p.vx,p.vy,p.vz))
f.close()

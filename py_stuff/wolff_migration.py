import numpy
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec

import sys
sys.path.insert(0, '/Users/jrobinson/grav_cloud/python_stuff')
import py_func as pf

N=100

a0=1.0
af=2.0
e0=0.0
i0=0.0

ta=1.0
te=1.0
ti=1.0

t0=0.0
tf=1.0
t=numpy.linspace(t0,tf,N)

aN=af+(a0-af)*numpy.exp(-t/ta)
eN=e0*numpy.exp(-t/te)
iN=i0*numpy.exp(-t/ti)

fig = pyplot.figure() #open figure once

gs = gridspec.GridSpec(3,1)
ax1 = pyplot.subplot(gs[0,0])
ax2 = pyplot.subplot(gs[1,0],sharex=ax1)
ax3 = pyplot.subplot(gs[2,0],sharex=ax1)

ax1.plot(t,aN)
ax2.plot(t,eN)
ax3.plot(t,iN)

ax1.set_ylabel('a')
# ax1.set_xlabel('t')
ax2.set_ylabel('e')
# ax2.set_xlabel('t')
ax3.set_ylabel('i')
ax3.set_xlabel('t')

pyplot.show()

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
import matplotlib.image as mpimg

ims_1=next(os.walk('planet_elements_plotter'))[2]
ims_1 = [ i for i in ims_1 if i.endswith(".png") and i.startswith("planet_elements_plot_group") ]
ims_1 = [ '{}/{}'.format('planet_elements_plotter',i) for i in ims_1 ]
ims_1.sort()
print ims_1

ims_2=next(os.walk('orbital_evolution'))[2]
ims_2 = [ i for i in ims_2 if i.endswith(".png") and i.startswith("orbital_evolution_helio_group") ]
ims_2 = [ '{}/{}'.format('orbital_evolution',i) for i in ims_2 ]
ims_2.sort()
print ims_2

ims_3=next(os.walk('orbital_evolution'))[2]
ims_3 = [ i for i in ims_3 if i.endswith(".png") and i.startswith("orbital_evolution2_group") ]
ims_3 = [ '{}/{}'.format('orbital_evolution',i) for i in ims_3 ]
ims_3.sort()
print ims_3

print len(ims_1)

fig = pyplot.figure()
gs = gridspec.GridSpec(6, 3)

ax1 = pyplot.subplot(gs[0,0])
ax2 = pyplot.subplot(gs[0,1])
ax3 = pyplot.subplot(gs[0,2])
ax4 = pyplot.subplot(gs[1,0])
ax5 = pyplot.subplot(gs[1,1])
ax6 = pyplot.subplot(gs[1,2])
ax7 = pyplot.subplot(gs[2,0])
ax8 = pyplot.subplot(gs[2,1])
ax9 = pyplot.subplot(gs[2,2])
ax10 = pyplot.subplot(gs[3,0])
ax11 = pyplot.subplot(gs[3,1])
ax12 = pyplot.subplot(gs[3,2])
ax13 = pyplot.subplot(gs[4,0])
ax14 = pyplot.subplot(gs[4,1])
ax15 = pyplot.subplot(gs[4,2])
ax16 = pyplot.subplot(gs[5,0])
ax17 = pyplot.subplot(gs[5,1])
ax18 = pyplot.subplot(gs[5,2])

a = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15,ax16,ax17,ax18]

k=0
for j in range(6):

    img1=mpimg.imread(ims_1[j])
    img2=mpimg.imread(ims_2[j])
    img3=mpimg.imread(ims_3[j])

    a[k].imshow(img1)
    a[k+1].imshow(img2)
    a[k+2].imshow(img3)

    a[k].axes.get_xaxis().set_visible(False)
    a[k].axes.get_yaxis().set_visible(False)
    a[k+1].axes.get_xaxis().set_visible(False)
    a[k+1].axes.get_yaxis().set_visible(False)
    a[k+2].axes.get_xaxis().set_visible(False)
    a[k+2].axes.get_yaxis().set_visible(False)

    k+=3

pyplot.tight_layout()
pyplot.show()

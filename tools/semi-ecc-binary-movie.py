import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import collections as cl
import sys
import multiple

filename='tlst'
dout='semi-ecc'
tscale=1.0
rscale=1.0
mscale=1.0
vscale=1.0
amin=1e-7
amax=0.1

narg=len(sys.argv)-1
if (narg>=1):
    filename = sys.argv[1]
    if (filename=='-h'): 
        print "arguments: snapshot list (suffix[time]) [tlst]; output file name [semi-ecc]; time scale [1]; r scale [1]; v scale [1]; m scale [1]; semi-min [1e-7]; semi-max [0.1] "
        quit()
if (narg>=2):    dout=sys.argv[2]
if (narg>=3):    tscale=float(sys.argv[3])
if (narg>=4):    rscale=float(sys.argv[4])
if (narg>=5):    vscale=float(sys.argv[5])
if (narg>=6):    mscale=float(sys.argv[6])
if (narg>=7):    amin=float(sys.argv[7])
if (narg>=8):    amax=float(sys.argv[8])

fl = open(filename,'r')
path = fl.read()
path = path.split()

fig,axes=plt.subplots(1,2,figsize=(8,4))

for i in path:
    dat = multiple.ParticleArray(np.loadtxt('data.'+str(i),skiprows=1))
    p1,p2,semi,ecc,r = multiple.findPair(dat)
    sel = (semi>amin) & (semi<amax)
    
    axes[0].plot(semi[sel],ecc[sel],'+',color='k')
    axes[0].set_title('T= %10.3f' % (float(i)*tscale))
    axes[0].set_xscale('log')
    axes[0].set_xlim(amin,amax)
    axes[0].set_ylim(0,1)
    axes[0].set_xlabel('semi')
    axes[0].set_ylabel('ecc')
    axes[0].legend(loc='upper left')

    r= np.sqrt((dat.pos*dat.pos).sum(axis=1))
    ekin= dat.m*(dat.vel*dat.vel).sum(axis=1)
    axes[1].plot(r,ekin,'.',color='k')
    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    axes[1].set_xlabel('r')
    axes[1].set_ylabel('ekin')
    axes[1].set_xlim(1e-7,10)
    axes[1].set_ylim(1e-1,1e4)

    fig.savefig(dout+'.'+str(i)+'.png',bbox_inches = "tight")
    axes[0].clear()
    axes[1].clear()

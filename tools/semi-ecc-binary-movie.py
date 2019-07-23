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
rmin=1e-7
rmax=10
ekmin=1e-1
ekmax=100

narg=len(sys.argv)-1
if (narg>=1):
    filename = sys.argv[1]
    if (filename=='-h'): 
        print "arguments: snapshot list (suffix[time]) [tlst]; output file name [semi-ecc]; time scale [1]; r scale [1]; v scale [1]; m scale [1]; semi-min [1e-7]; semi-max [0.1]; rmin [1e-7]; rmax [10]; ekmin [1e-1]; ekmax [100] "
        quit()
if (narg>=2):    dout=sys.argv[2]
if (narg>=3):    tscale=float(sys.argv[3])
if (narg>=4):    rscale=float(sys.argv[4])
if (narg>=5):    vscale=float(sys.argv[5])
if (narg>=6):    mscale=float(sys.argv[6])
if (narg>=7):    amin=float(sys.argv[7])
if (narg>=8):    amax=float(sys.argv[8])
if (narg>=9):    rmin=float(sys.argv[9])
if (narg>=10):   rmax=float(sys.argv[10])
if (narg>=11):   ekmin=float(sys.argv[11])
if (narg>=12):   ekmax=float(sys.argv[12])


fl = open(filename,'r')
path = fl.read()
path = path.split()

fig,axes=plt.subplots(2,2,figsize=(8,8))

for i in path:
    dat = multiple.ParticleArray(np.loadtxt('data.'+str(i),skiprows=1))
    p1,p2,semi,ecc,r = multiple.findPair(dat)
    sel = (semi>amin) & (semi<amax)
    Nbin =sel.sum()
    
    axes[0][0].plot(semi[sel],ecc[sel],'+',color='k')
    axes[0][0].set_title('T= %10.3f Nbin=%d' % (float(i)*tscale,Nbin))
    axes[0][0].set_xscale('log')
    axes[0][0].set_xlim(amin,amax)
    axes[0][0].set_ylim(0,1)
    axes[0][0].set_xlabel('semi')
    axes[0][0].set_ylabel('ecc')
    axes[0][0].legend(loc='upper left')

    axes[1][0].hist(semi[sel],bins=np.logspace(np.log10(amin),np.log10(amax),50))
    axes[1][0].set_xlabel('semi')
    axes[1][0].set_xscale('log')

    axes[0][1].hist(ecc[sel],bins=np.linspace(0,1,50))
    axes[0][1].set_xlabel('ecc')

    rcm = np.array([(dat.pos[:,k]*dat.m).sum() for k in range(3)])
    vcm = np.array([(dat.vel[:,k]*dat.m).sum() for k in range(3)])
    
    r= np.sqrt(((dat.pos-rcm)**2).sum(axis=1))
    ekin= dat.m*((dat.vel-vcm)**2).sum(axis=1)

    xbins=np.logspace(np.log10(rmin),np.log10(rmax),100)
    ybins=np.logspace(np.log10(ekmin),np.log10(ekmax),100)
    lbin=[xbins,ybins]
    counts, _, _ = np.histogram2d(r,ekin, bins=lbin)
    #counts,ybins,xbins,image = plt.hist2d(nbperi,nbecc,bins=lbin,norm=cl.LogNorm())

    axes[1][1].pcolormesh(xbins, ybins, np.log(counts.T))
    axes[1][1].set_xscale('log')
    axes[1][1].set_yscale('log')
    axes[1][1].set_xlabel('r')
    axes[1][1].set_ylabel('ekin')
    axes[1][1].set_xlim(rmin,rmax)
    axes[1][1].set_ylim(ekmin,ekmax)

    print 'T=',i,' N=',dat.m.size,' Nbin=', Nbin
    print 'Semi range:', semi.min(), semi.max()
    print 'R range:', r.min(),r.max()
    print 'Ekin range:', ekin.min(),ekin.max()

    fig.savefig(dout+'.'+str(i)+'.png',bbox_inches = "tight")
    axes[0][0].clear()
    axes[0][1].clear()
    axes[1][0].clear()
    axes[1][1].clear()

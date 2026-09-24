#!/usr/bin/env python3
# tsv110 benchmark figures + summary stats
import os, math, statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D = HERE + os.sep
FIG = os.path.join(HERE, '..', 'figs') + os.sep
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'figure.dpi':130, 'savefig.dpi':150, 'font.size':9})

def load(fn):
    rows=[]
    for line in open(D+fn):
        line=line.strip()
        if not line or line.startswith('#'): continue
        p=line.split(',')
        rows.append(dict(kernel=p[0],variant=p[1],ni=int(p[2]),nj=int(p[3]),t=float(p[4]),ips=float(p[5])))
    return rows

g=load('kernel_g.csv'); a=load('kernel_a.csv'); n=load('kernel_n.csv')
def get(rows,k,v,ni,nj):
    for r in rows:
        if r['kernel']==k and r['variant']==v and r['ni']==ni and r['nj']==nj: return r
    raise KeyError((k,v,ni,nj))
nis=[4,16,64,256,1024]; njs=[8,32,128,512,2048]
kernels=['nb','epep','mono','quad']
ktitles={'nb':'neighbour search','epep':'EP-EP force','mono':'EP-SP monopole','quad':'EP-SP quadrupole'}
colors={'nosimd':'#444444','neon4':'#1f77b4','neon1':'#ff7f0e'}

summary=[]
for k in kernels:
    sp=[]
    for ni in nis:
        for nj in njs:
            t0=get(n,'nb' if k=='nb' else k,'nosimd',ni,nj)['t'] if k=='nb' else get(n,k,'nosimd',ni,nj)['t']
            t4=get(n,k,'neon4',ni,nj)['t']; t1=get(n,k,'neon1',ni,nj)['t']
            sp.append(t0/min(t4,t1))
    gm=math.exp(sum(math.log(x) for x in sp)/len(sp))
    summary.append((k,gm,min(sp),max(sp)))
    print(f'SPEEDUP {k}: geomean={gm:.2f} min={min(sp):.2f} max={max(sp):.2f}')

# ---------- fig1: CPU peak ----------
fig,ax=plt.subplots(1,2,figsize=(9,3.2))
labels=['FP32 1T','FP32 24T','FP64 1T','FP64 24T']; vals=[20.76,479.29,9.83,234.87]
ax[0].bar(labels,vals,color=['#1f77b4','#1f77b4','#ff7f0e','#ff7f0e'],alpha=.85)
for i,v in enumerate(vals): ax[0].text(i,v*1.03,f'{v:.1f}',ha='center',fontsize=8)
ax[0].set_ylabel('GFLOP/s'); ax[0].set_title('NEON FMA peak (measured)')
ax[0].set_ylim(0,560)
ops=['vrsqrte\n(2-lane)','vrecpe\n(2-lane)','vsqrt\n(2-lane)','vdiv\n(2-lane)','vfma\n(4-lane)']
ov=[1.735,1.253,9.041,8.853,0.385]
ax[1].bar(ops,ov,color='#2ca02c',alpha=.85)
ax[1].set_ylabel('ns / vector instruction'); ax[1].set_yscale('log')
ax[1].set_title('Instruction cost (independent, 1 thread)')
for i,v in enumerate(ov): ax[1].text(i,v*1.15,f'{v:.3g}',ha='center',fontsize=8)
plt.tight_layout(); plt.savefig(FIG+'fig1_cpu_peak.png'); plt.close()

# ---------- fig2: latency + bandwidth ----------
fig,ax=plt.subplots(1,2,figsize=(9,3.2))
levels=['L1 16K','L2 256K','L3 8M','DRAM 512M']
old=[2.06,5.58,42.95,151.58]; new4k=[1.54,4.66,40.27,145.11]; thp=[1.54,4.66,28.23,82.02]
x=np.arange(4); w=0.26
ax[0].bar(x-w,old,w,label='4-byte chain, unpinned',color='#999999')
ax[0].bar(x,new4k,w,label='pointer chain, pinned, 4K pages',color='#1f77b4')
ax[0].bar(x+w,thp,w,label='pointer chain, pinned, THP(2M)',color='#2ca02c')
ax[0].set_yscale('log'); ax[0].set_xticks(x); ax[0].set_xticklabels(levels)
ax[0].set_ylabel('ns/access'); ax[0].set_title('Memory latency (pointer chase, pinned)')
ax[0].legend(fontsize=7)
bw=['DRAM copy\n1T','DRAM triad\n1T','DRAM copy\n24T','DRAM triad\n24T']
bv=[6.7,12.3,14.5,30.9]
ax[1].bar(bw,bv,color='#d62728',alpha=.85)
for i,v in enumerate(bv): ax[1].text(i,v*1.03,f'{v:.1f}',ha='center',fontsize=8)
ax[1].set_ylabel('GB/s'); ax[1].set_title('DRAM bandwidth (256 MB arrays)')
plt.tight_layout(); plt.savefig(FIG+'fig2_latency_bw.png'); plt.close()

# ---------- fig3: speedup heatmaps ----------
fig,axs=plt.subplots(2,2,figsize=(9.5,7))
vmin,vmax=0,5
for ax,k in zip(axs.ravel(),kernels):
    M=np.zeros((len(nis),len(njs)))
    for i,ni in enumerate(nis):
        for j,nj in enumerate(njs):
            t0=get(n,k,'nosimd',ni,nj)['t']; t4=get(n,k,'neon4',ni,nj)['t']; t1=get(n,k,'neon1',ni,nj)['t']
            M[i,j]=t0/min(t4,t1)
    im=ax.imshow(M,cmap='viridis',vmin=vmin,vmax=vmax,origin='lower',aspect='auto')
    for i in range(len(nis)):
        for j in range(len(njs)):
            ax.text(j,i,f'{M[i,j]:.1f}',ha='center',va='center',
                    color='white' if M[i,j]<3.5 else 'black', fontsize=8)
    ax.set_xticks(range(len(njs))); ax.set_xticklabels(njs)
    ax.set_yticks(range(len(nis))); ax.set_yticklabels(nis)
    ax.set_xlabel('n_j'); ax.set_ylabel('n_i')
    ax.set_title(f'{ktitles[k]}',fontsize=10)
fig.suptitle('NEON best / NoSimd speedup (ns per node-pair call)', fontsize=11)
fig.colorbar(im,ax=axs,shrink=.8,label='speedup')
plt.savefig(FIG+'fig3_kernel_speedup.png'); plt.close()

# ---------- fig4: per-interaction time lines (ni=1024) ----------
fig,axs=plt.subplots(2,2,figsize=(9.5,7))
for ax,k in zip(axs.ravel(),kernels):
    ni=1024
    for v in ['nosimd','neon4','neon1']:
        ys=[get(n,k,v,ni,nj)['t']/(ni*nj) for nj in njs]
        ax.plot(njs,ys,'o-',label=v,color=colors[v],ms=4)
    yg=[get(g,k,'nosimd',ni,nj)['t']/(ni*nj) for nj in njs]
    ya=[get(a,k,'nosimd',ni,nj)['t']/(ni*nj) for nj in njs]
    ax.plot(njs,yg,'s--',label='nosimd generic',color='#bbbbbb',ms=3)
    ax.plot(njs,ya,'^--',label='nosimd -mcpu=tsv110',color='#e377c2',ms=3)
    ax.set_xscale('log',base=2); ax.set_yscale('log')
    ax.set_xlabel('n_j'); ax.set_ylabel('ns / interaction')
    ax.set_title(ktitles[k],fontsize=10); ax.grid(alpha=.3,which='both')
    if k=='epep': ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig(FIG+'fig4_kernel_lines.png'); plt.close()

# ---------- fig5: e2e ----------
e2e={'base':[],'auto':[],'neon':[]}
for line in open(D+'e2e_timing.txt'):
    p=line.split()
    if len(p)>=4 and p[0]=='RESULT':
        e2e[p[1]].append(float([x for x in p if x.startswith('wall_s=')][0].split('=')[1]))
fig,ax=plt.subplots(figsize=(5.2,3.4))
names=['base','auto','neon']; med=[statistics.median(e2e[x]) for x in names]
ax.bar(names,med,color=['#999999','#e377c2','#2ca02c'],alpha=.8)
for i,x in enumerate(names):
    ax.scatter([i]*len(e2e[x]),e2e[x],color='black',zorder=3,s=18)
    ax.text(i,med[i]+1.5,f'{med[i]:.1f}s',ha='center',fontsize=9)
ax.set_ylabel('wall time (s)'); ax.set_ylim(0,95)
ax.set_title('PeTar demo, 2000 stars, t=10 Myr, 24 MPI x 1 OMP\n(3 runs; speedup %.2fx)'%(med[0]/med[2]))
plt.tight_layout(); plt.savefig(FIG+'fig5_e2e.png'); plt.close()
print(f'E2E base={med[0]} auto={med[1]} neon={med[2]} speedup={med[0]/med[2]:.2f}')

# ---------- fig6: rsqrt accuracy ----------
fig,ax=plt.subplots(figsize=(5.6,3.4))
meth=['vrsqrte','+1 NR','+2 NR','+1 cubic','f64 +2 NR']
err=[3.277e-3,1.614e-5,1.439e-7,1.663e-7,3.898e-10]
ax.bar(meth,err,color=['#d62728','#ff7f0e','#2ca02c','#1f77b4','#9467bd'],alpha=.85)
ax.set_yscale('log'); ax.set_ylabel('max relative error')
ax.axhline(7e-3,color='k',ls='--',lw=1); ax.text(3.6,9e-3,'kernel tolerance 7e-3',fontsize=7)
for i,v in enumerate(err): ax.text(i,v*1.3,f'{v:.2e}',ha='center',fontsize=7)
ax.set_title('1/sqrt estimate + Newton/cubic correction (F32/F64)')
plt.tight_layout(); plt.savefig(FIG+'fig6_rsqrt.png'); plt.close()

print('FIGS_DONE')

# ---------- fig7: N scaling ----------
scal={}
for line in open(D+'scaling.txt'):
    p=line.split()
    if len(p)>=3 and p[0]=='SCAL':
        N=int(p[1][2:]); name=p[2]
        kv={x.split('=')[0]:x.split('=')[1] for x in p[3:]}
        scal[(N,name)]=dict(wall=float(kv['wall'].rstrip('s')),total=float(kv['total_step'].rstrip('s')),
                            tf=float(kv['tree_force'].rstrip('s')),cf=float(kv['fdps_calcforce'].rstrip('s')))
Ns=sorted(set(N for N,_ in scal))
fig,ax=plt.subplots(1,2,figsize=(9,3.4))
for name,c in [('base','#999999'),('neon','#2ca02c')]:
    ys=[scal[(N,name)]['wall'] for N in Ns]
    ax[0].plot(Ns,ys,'o-',color=c,label=name)
    for N,y in zip(Ns,ys): ax[0].annotate(f'{y:.1f}s',(N,y),textcoords='offset points',xytext=(0,6),ha='center',fontsize=7)
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].minorticks_off()
ax[0].set_xticks(Ns); ax[0].set_xticklabels(Ns)
ax[0].set_xlabel('N (single stars, no binaries)'); ax[0].set_ylabel('wall time for t=2 Myr (s)')
ax[0].set_title('End-to-end scaling (24 MPI x 1 OMP)'); ax[0].legend()
sp_e2e=[scal[(N,'base')]['wall']/scal[(N,'neon')]['wall'] for N in Ns]
sp_cf=[scal[(N,'base')]['cf']/scal[(N,'neon')]['cf'] for N in Ns]
x=np.arange(len(Ns)); w=0.35
ax[1].bar(x-w/2,sp_e2e,w,label='end-to-end wall',color='#1f77b4')
ax[1].bar(x+w/2,sp_cf,w,label='FDPS calc_force (local)',color='#ff7f0e')
for i,(a,b) in enumerate(zip(sp_e2e,sp_cf)):
    ax[1].text(i-w/2,a+.03,f'{a:.2f}',ha='center',fontsize=8)
    ax[1].text(i+w/2,b+.03,f'{b:.2f}',ha='center',fontsize=8)
ax[1].set_xticks(x); ax[1].set_xticklabels(Ns); ax[1].set_xlabel('N')
ax[1].set_ylabel('speedup (NEON / base)'); ax[1].set_ylim(0,4)
ax[1].set_title('Speedup vs N'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(FIG+'fig7_scaling.png'); plt.close()
print('fig7 done')


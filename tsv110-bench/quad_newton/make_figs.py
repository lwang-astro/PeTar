#!/usr/bin/env python3
# figures for the quad Newton study
import os, math, statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
F = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figs')
os.makedirs(F, exist_ok=True)
plt.rcParams.update({'figure.dpi':130,'savefig.dpi':150,'font.size':9})

C = {'cubic':'#1f77b4','newton':'#ff7f0e','neon_cubic':'#1f77b4','neon_newton':'#ff7f0e','nosimd':'#777777'}

def parse_kv(s):
    out={}
    for tok in s.split(','):
        if '=' in tok:
            k,v=tok.split('=',1); out[k]=v
    return out

# ---------- 1) r_inv accuracy ----------
rs={}
for m in ['cubic','newton']:
    hist=[]; stat=None
    for line in open(os.path.join(D,f'rsinv_{m}.csv')):
        p=line.strip().split(',')
        if p[0]=='STAT': stat=parse_kv(line.strip())
        elif p[0]=='HIST': hist.append((float(p[2]),int(p[3])))
    hist.sort()
    x=[h[0] for h in hist]; n=[h[1] for h in hist]
    tot=sum(n)
    cs=np.cumsum(n)/tot
    rs[m]=(x,cs,stat)
    print(m, stat)

fig,ax=plt.subplots(figsize=(6,3.6))
for m,c in [('cubic','#1f77b4'),('newton','#ff7f0e')]:
    x,cs,st=rs[m]
    ax.semilogx([10**v for v in x],cs,color=c,label=f"{m} (max={float(st['max']):.2e}, p99={float(st['p99']):.2e})")
ax.axvline(7e-3,color='k',ls='--',lw=1); ax.text(8e-3,0.05,'force tolerance 7e-3',fontsize=7)
ax.set_xlim(1e-9,1e-2); ax.set_ylim(0,1.02)
ax.set_xlabel('relative error of 1/sqrt(x)'); ax.set_ylabel('CDF')
ax.set_title('r_inv accuracy: cubic correction vs +half-Newton (F32)')
ax.grid(alpha=.3); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(F,'fig1_rsinv_cdf.png')); plt.close()

# ---------- 2) per-particle force error CDF ----------
def load_dump(name):
    ea=[]; ep=[]
    for seed in [1,2,3]:
        fn=os.path.join(D,f'errdump_{name}_{seed}.csv')
        if not os.path.exists(fn): continue
        for line in open(fn):
            p=line.strip().split(',')
            if p[0]=='i': continue
            ea.append(float(p[1])); ep.append(float(p[2]))
    return np.array(ea),np.array(ep)

fig,axs=plt.subplots(1,2,figsize=(9,3.4))
for name,c in [('cubic','#1f77b4'),('newton','#ff7f0e')]:
    ea,ep=load_dump(name)
    for ax,arr,lbl in [(axs[0],ea,'relative |acc| error'),(axs[1],ep,'relative potential error')]:
        s=np.sort(arr); ax.semilogx(s, np.arange(1,len(s)+1)/len(s), color=c, label=name)
    print('dump',name,'n=',len(ea),'acc max',ea.max())
axs[0].set_title('EP-SP quadrupole: |acc| error CDF'); axs[1].set_title('potential error CDF')
for ax in axs:
    ax.axvline(7e-3,color='k',ls='--',lw=1); ax.set_xlim(1e-8,1e-1); ax.set_ylim(0,1.02)
    ax.set_xlabel('relative error'); ax.set_ylabel('CDF'); ax.grid(alpha=.3); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(F,'fig2_force_error_cdf.png')); plt.close()

# ---------- 3) error bars across dynamic-range configs ----------
rows=[]
for line in open(os.path.join(D,'errors.csv')):
    p=line.strip().split(',')
    if p[0]!='ERR': continue
    rows.append(p)
cfgs={}
for p in rows:
    variant,scale,off=p[1],float(p[4]),float(p[5])
    key=(scale,off)
    cfgs.setdefault(key,{'cubic':[],'newton':[]})
    tag='newton' if 'newton' in variant else 'cubic'
    cfgs[key][tag].append((float(p[7]),float(p[13])))  # acc_max, pot_max
fig,ax=plt.subplots(figsize=(6.5,3.4))
keys=sorted(cfgs.keys())
labels=[f"scale={int(s)}\noffset={int(o)}" for s,o in keys]
x=np.arange(len(keys)); w=0.2
for i,(tag,c) in enumerate([('cubic','#1f77b4'),('newton','#ff7f0e')]):
    acc=[max(v[0] for v in cfgs[k][tag]) for k in keys]
    pot=[max(v[1] for v in cfgs[k][tag]) for k in keys]
    ax.bar(x+(i-0.5)*2*w,acc,w,color=c,label=f'{tag} |acc| max')
    ax.bar(x+(i+0.5)*2*w,pot,w,color=c,alpha=.45,label=f'{tag} pot max')
ax.set_yscale('log'); ax.set_xticks(x); ax.set_xticklabels(labels,fontsize=8)
ax.axhline(7e-3,color='k',ls='--',lw=1)
ax.set_ylabel('max relative error'); ax.set_title('EP-SP quadrupole error vs NoSimd F64 (3 seeds, max)')
ax.legend(fontsize=7,ncol=2); ax.grid(alpha=.3,axis='y')
plt.tight_layout(); plt.savefig(os.path.join(F,'fig3_error_bars.png')); plt.close()

# ---------- 4) timing ----------
t={}
for line in open(os.path.join(D,'timing.csv')):
    p=line.strip().split(',')
    if p[0]!='TIME': continue
    v,ni,ns,tn=p[1],int(p[2]),int(p[3]),float(p[4])
    t.setdefault((v,ni,ns),[]).append(tn)
tv={k:statistics.median(v) for k,v in t.items()}
nis=[4,16,64,256,1024]; nss=[8,32,128,512,2048]
fig,axs=plt.subplots(1,2,figsize=(9.5,3.6))
ni_show=256
for v,c,lab in [('nosimd','#777777','NoSimd F64'),('neon_cubic','#1f77b4','NEON cubic'),('neon_newton','#ff7f0e','NEON +half-Newton')]:
    ys=[tv[(v,ni_show,ns)]/(ni_show*ns) for ns in nss]
    axs[0].plot(nss,ys,'o-',color=c,label=lab,ms=4)
axs[0].set_xscale('log',base=2); axs[0].set_yscale('log')
axs[0].set_xlabel('n_j (n_i=256)'); axs[0].set_ylabel('ns / interaction')
axs[0].set_title('quadrupole kernel time'); axs[0].grid(alpha=.3,which='both'); axs[0].legend(fontsize=7)
ratio=np.zeros((len(nis),len(nss)))
for i,ni in enumerate(nis):
    for j,ns in enumerate(nss):
        ratio[i,j]=tv[('neon_newton',ni,ns)]/tv[('neon_cubic',ni,ns)]
im=axs[1].imshow(ratio,cmap='coolwarm',vmin=0.97,vmax=1.03,origin='lower',aspect='auto')
for i in range(len(nis)):
    for j in range(len(nss)):
        axs[1].text(j,i,f'{ratio[i,j]:.3f}',ha='center',va='center',fontsize=7,
                    color='black' if abs(ratio[i,j]-1)<0.02 else 'white')
axs[1].set_xticks(range(len(nss))); axs[1].set_xticklabels(nss)
axs[1].set_yticks(range(len(nis))); axs[1].set_yticklabels(nis)
axs[1].set_xlabel('n_j'); axs[1].set_ylabel('n_i')
axs[1].set_title('kernel time ratio: +half-Newton / cubic')
fig.colorbar(im,ax=axs[1],shrink=.85,label='ratio')
plt.tight_layout(); plt.savefig(os.path.join(F,'fig4_timing.png')); plt.close()
gm=math.exp(sum(math.log(x) for x in ratio.ravel())/ratio.size)
print(f'TIMING geomean ratio newton/cubic = {gm:.4f}  (cells min={ratio.min():.3f} max={ratio.max():.3f})')
print('FIGS_DONE')

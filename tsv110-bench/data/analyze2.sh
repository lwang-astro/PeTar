#!/bin/bash
# extra analyses: orientation winners + GCC autovectorization evidence
cd $HOME/tsv110-bench/kernel
source $HOME/tests/PeTar-project/env.sh

python3 - <<'PYEOF'
rows=[]
for line in open('../data/kernel_n.csv'):
    line=line.strip()
    if not line or line.startswith('#'): continue
    p=line.split(','); rows.append((p[0],p[1],int(p[2]),int(p[3]),float(p[4])))
def get(k,v,ni,nj):
    for r in rows:
        if r[0]==k and r[1]==v and r[2]==ni and r[3]==nj: return r[4]
nis=[4,16,64,256,1024]; njs=[8,32,128,512,2048]
for k in ['epep','quad']:
    print(k,'winner matrix (4=neon4, 1=neon1):')
    for ni in nis:
        print(' ',ni, ['4' if get(k,'neon4',ni,nj)<get(k,'neon1',ni,nj) else '1' for nj in njs])
PYEOF

PETAR=$HOME/tests/PeTar-project/PeTar
FDPS=$HOME/tests/PeTar-project/FDPS
SDAR=$HOME/tests/PeTar-project/SDAR
g++ -O3 -mcpu=tsv110 -std=c++17 -DUSE_QUAD -I$PETAR/src -I$FDPS/src -I$SDAR/src \
    -fopt-info-vec=vec.txt -fopt-info-vec-missed=vecmiss.txt -c bench_kernels.cxx -o /tmp/bk.o 2>/dev/null
echo "== vectorized lines in soft_force.hpp: $(grep -c 'soft_force' vec.txt) =="
grep 'soft_force' vec.txt | head -5
echo "== missed loops in soft_force.hpp: $(grep -c 'soft_force' vecmiss.txt) =="
grep 'soft_force' vecmiss.txt | head -5
echo "== vector fmla count in NoSimd EP-EP function (bench_a) =="
read addr size <<< $(nm -S --demangle bench_a | grep -m1 'CalcForceEpEpWithLinearCutoffNoSimd::operator' | awk '{print $1, $2}')
echo "addr=$addr size=$size"
if [ -n "$addr" ]; then
  end=$(printf '%x' $((0x$addr + 0x$size)))
  objdump -d --start-address=0x$addr --stop-address=0x$end bench_a | grep -c 'fmla\|fmul\|fadd\|fdiv\|fsqrt' || true
  objdump -d --start-address=0x$addr --stop-address=0x$end bench_a | grep -c 'fmla v\|fmul v\|fadd v\|fdiv v\|fsqrt v' || true
fi

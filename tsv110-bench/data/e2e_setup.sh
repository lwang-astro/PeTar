#!/bin/bash
# set up two scratch PeTar builds: -mcpu autovec and NEON-integrated
set -e
ROOT=$HOME/tsv110-bench/e2e
SRC=$HOME/tests/PeTar-project/PeTar
FDPS=$HOME/tests/PeTar-project/FDPS
SDAR=$HOME/tests/PeTar-project/SDAR
mkdir -p $ROOT
cd $ROOT
if [ ! -d PeTar-auto ]; then
  cp -a $SRC PeTar-auto
  cp -a $SRC PeTar-neon
fi
rm -rf PeTar-auto/.git PeTar-neon/.git PeTar-auto/build PeTar-neon/build
cd $ROOT/PeTar-auto
./configure --prefix=$ROOT/install-auto --with-arch=tsv110 --with-fdps-prefix=$FDPS --with-sdar-prefix=$SDAR > configure.out 2>&1
cd $ROOT/PeTar-neon
cp $HOME/tsv110-bench/kernel/force_tsv110.hpp src/force_tsv110.hpp
python3 - <<'EOF'
p='src/petar.hpp'
s=open(p).read()
inc='#ifdef USE_FUGAKU\n#include "force_fugaku.hpp"\n#endif'
assert inc in s, "include anchor not found"
s=s.replace(inc, inc+'\n#ifdef USE_NEON_KERNEL\n#include "force_tsv110.hpp"\n#endif',1)
nb='#elif USE_FUGAKU\n        tree_nb.calcForceAllAndWriteBack(SearchNeighborEpEpFugaku(), system_soft, dinfo);\n'
assert nb in s, "nb anchor not found"
s=s.replace(nb, nb+'#elif defined(USE_NEON_KERNEL)\n        tree_nb.calcForceAllAndWriteBack(tsv110::SearchNeighborEpEpNeon(), system_soft, dinfo);\n',1)
anchor='#elif USE_SIMD // end use_gpu\n'
assert anchor in s, "force anchor not found"
neon_block='''#elif defined(USE_NEON_KERNEL)
        PS::F64 eps2 = EPISoft::eps*EPISoft::eps;
        PS::F64 rout2 = EPISoft::r_out*EPISoft::r_out;
        PS::F64 G= ForceSoft::grav_const;
        tree_soft.calcForceAllAndWriteBack(tsv110::CalcForceEpEpWithLinearCutoffNeon(eps2, rout2, G),
#ifdef USE_QUAD
                                           tsv110::CalcForceEpSpQuadNeon<PS::SPJQuadrupoleInAndOut>(eps2, G),
#else
                                           tsv110::CalcForceEpSpMonoNeon<PS::SPJMonopoleInAndOut>(eps2, G),
#endif
                                           system_soft,
                                           dinfo);
        
'''
s=s.replace(anchor, neon_block+anchor,1)
open(p,'w').write(s)
print("petar.hpp patched OK")
EOF
./configure --prefix=$ROOT/install-neon --with-arch=tsv110 --with-fdps-prefix=$FDPS --with-sdar-prefix=$SDAR > configure.out 2>&1
grep -q 'USE_NEON_KERNEL' Makefile || echo 'CXXFLAGS += -D USE_NEON_KERNEL' >> Makefile
echo setup_done

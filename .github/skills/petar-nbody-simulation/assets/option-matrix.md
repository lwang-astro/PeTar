# PeTar Option Matrix (Auto-generated)

This file is generated from commands discovered in PATH via <binary> -h option extraction.

Discovery rule:
- scan PATH for executable names starting with 'petar'
- classify helper tools by suffix (*.hard.debug, *.format.transfer)
- classify solver binaries by required core options (-u, -t, -o)

## Solver Binaries Included

- petar.mpi.omp.avx2.64b (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.64b
- petar.mpi.omp.avx2.agama (62 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.agama
- petar.mpi.omp.avx2.bse.agama (94 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama
- petar.mpi.omp.avx2.bse.galpy.gasdrag (107 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.bse.galpy (96 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- petar.mpi.omp.avx2.bse.gasdrag (102 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.gasdrag
- petar.mpi.omp.avx2.bse (91 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse
- petar.mpi.omp.avx2.bse.pnhermite (93 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.pnhermite
- petar.mpi.omp.avx2.dsm.galpy.gasdrag (93 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.galpy.gasdrag
- petar.mpi.omp.avx2.dsm.gasdrag (88 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.gasdrag
- petar.mpi.omp.avx2.galpy.mp (64 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.mp
- petar.mpi.omp.avx2.galpy (64 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy
- petar.mpi.omp.avx2.gasdrag (70 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gasdrag
- petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag (107 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
- petar.mpi.omp.avx2.gpu.bse.galpy (96 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy
- petar.mpi.omp.avx2.gpu.bse (91 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse
- petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag (93 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag
- petar.mpi.omp.avx2.gpu (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu
- petar.mpi.omp.avx2.kdk.g (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.kdk.g
- petar.mpi.omp.avx2.kdkdk4.g (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.kdkdk4.g
- petar.mpi.omp.avx2 (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2
- petar.mpi.omp.avx2.pnall (62 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.pnall

## Requirement-Driven Solver Filters

No fixed priority across physics scenarios. Select filters by user requirements first, then choose by performance.

Configure mapping:
- --with-interrupt -> suffix tokens: merger | bse | mobse | bseEmp | dsm
- --with-external -> suffix tokens: galpy | agama
- --with-external-hard -> suffix token: gasdrag
- --with-pn -> suffix tokens: pn*

### Interrupt Module (--with-interrupt)

- petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
  - perf-tags: gpu,avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
- petar.mpi.omp.avx2.gpu.bse.galpy
  - perf-tags: gpu,avx2,omp,mpi
  - options: 96
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy
- petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag
  - perf-tags: gpu,avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag
- petar.mpi.omp.avx2.gpu.bse
  - perf-tags: gpu,avx2,omp,mpi
  - options: 91
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse
- petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.bse.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 102
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.gasdrag
- petar.mpi.omp.avx2.bse.galpy
  - perf-tags: avx2,omp,mpi
  - options: 96
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- petar.mpi.omp.avx2.bse.agama
  - perf-tags: avx2,omp,mpi
  - options: 94
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama

### Stellar Evolution (bse/mobse/bseEmp)

- petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
  - perf-tags: gpu,avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
- petar.mpi.omp.avx2.gpu.bse.galpy
  - perf-tags: gpu,avx2,omp,mpi
  - options: 96
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy
- petar.mpi.omp.avx2.gpu.bse
  - perf-tags: gpu,avx2,omp,mpi
  - options: 91
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse
- petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.bse.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 102
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.gasdrag
- petar.mpi.omp.avx2.bse.galpy
  - perf-tags: avx2,omp,mpi
  - options: 96
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- petar.mpi.omp.avx2.bse.agama
  - perf-tags: avx2,omp,mpi
  - options: 94
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama
- petar.mpi.omp.avx2.bse.pnhermite
  - perf-tags: avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.pnhermite

### Long-timescale External: Galpy (--with-external=galpy)

- petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
  - perf-tags: gpu,avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
- petar.mpi.omp.avx2.gpu.bse.galpy
  - perf-tags: gpu,avx2,omp,mpi
  - options: 96
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy
- petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag
  - perf-tags: gpu,avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag
- petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.bse.galpy
  - perf-tags: avx2,omp,mpi
  - options: 96
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- petar.mpi.omp.avx2.dsm.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.galpy.gasdrag
- petar.mpi.omp.avx2.galpy.mp
  - perf-tags: avx2,omp,mpi
  - options: 64
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.mp
- petar.mpi.omp.avx2.galpy
  - perf-tags: avx2,omp,mpi
  - options: 64
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy

### Long-timescale External: Agama (--with-external=agama)

- petar.mpi.omp.avx2.bse.agama
  - perf-tags: avx2,omp,mpi
  - options: 94
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama
- petar.mpi.omp.avx2.agama
  - perf-tags: avx2,omp,mpi
  - options: 62
  - path: /home/lwang/bin/petar.mpi.omp.avx2.agama

### Short-timescale External Hard: Gas Drag (--with-external-hard=gasdrag)

- petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
  - perf-tags: gpu,avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
- petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag
  - perf-tags: gpu,avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag
- petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.bse.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 102
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.gasdrag
- petar.mpi.omp.avx2.dsm.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.galpy.gasdrag
- petar.mpi.omp.avx2.dsm.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 88
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.gasdrag
- petar.mpi.omp.avx2.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 70
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gasdrag

### Post-Newtonian Relativity (--with-pn)

- petar.mpi.omp.avx2.bse.pnhermite
  - perf-tags: avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.pnhermite
- petar.mpi.omp.avx2.pnall
  - perf-tags: avx2,omp,mpi
  - options: 62
  - path: /home/lwang/bin/petar.mpi.omp.avx2.pnall

### Special-purpose Features (Non-default Recommendations)

These are for special user requests and should not be selected by default.

### High-precision Tree Force (--enable-64b)

- petar.mpi.omp.avx2.64b
  - perf-tags: avx2,omp,mpi
  - options: 59
  - path: /home/lwang/bin/petar.mpi.omp.avx2.64b

### High-precision Position Representation (--enable-mpfrc)

- petar.mpi.omp.avx2.galpy.mp
  - perf-tags: avx2,omp,mpi
  - options: 64
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.mp

### Debug Build: g (--with-debug=g)

- petar.mpi.omp.avx2.kdkdk4.g
  - perf-tags: avx2,omp,mpi
  - options: 59
  - path: /home/lwang/bin/petar.mpi.omp.avx2.kdkdk4.g
- petar.mpi.omp.avx2.kdk.g
  - perf-tags: avx2,omp,mpi
  - options: 59
  - path: /home/lwang/bin/petar.mpi.omp.avx2.kdk.g

### Debug Build: assert (--with-debug=assert)

- (no matching solver detected)


## Recommended Solver Candidates (Performance-first)

Priority order: gpu > avx512 > avx2 > omp > mpi
Apply this ordering after choosing a requirement-driven filter.
- #1 petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
  - perf-tags: gpu,avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag
- #2 petar.mpi.omp.avx2.gpu.bse.galpy
  - perf-tags: gpu,avx2,omp,mpi
  - options: 96
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy
- #3 petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag
  - perf-tags: gpu,avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag
- #4 petar.mpi.omp.avx2.gpu.bse
  - perf-tags: gpu,avx2,omp,mpi
  - options: 91
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse
- #5 petar.mpi.omp.avx2.gpu
  - perf-tags: gpu,avx2,omp,mpi
  - options: 59
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu
- #6 petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 107
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- #7 petar.mpi.omp.avx2.bse.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 102
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.gasdrag
- #8 petar.mpi.omp.avx2.bse.galpy
  - perf-tags: avx2,omp,mpi
  - options: 96
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- #9 petar.mpi.omp.avx2.bse.agama
  - perf-tags: avx2,omp,mpi
  - options: 94
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama
- #10 petar.mpi.omp.avx2.dsm.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.galpy.gasdrag
- #11 petar.mpi.omp.avx2.bse.pnhermite
  - perf-tags: avx2,omp,mpi
  - options: 93
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.pnhermite
- #12 petar.mpi.omp.avx2.bse
  - perf-tags: avx2,omp,mpi
  - options: 91
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse

## Helper Tool Binaries Excluded From Simulation Commands

- petar.mpi.omp.avx2.64b.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.64b.format.transfer
- petar.mpi.omp.avx2.64b.hard.debug (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.64b.hard.debug
- petar.mpi.omp.avx2.agama.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.agama.format.transfer
- petar.mpi.omp.avx2.agama.hard.debug (62 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.agama.hard.debug
- petar.mpi.omp.avx2.bse.agama.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama.format.transfer
- petar.mpi.omp.avx2.bse.agama.hard.debug (94 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama.hard.debug
- petar.mpi.omp.avx2.bse.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.format.transfer
- petar.mpi.omp.avx2.bse.galpy.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.format.transfer
- petar.mpi.omp.avx2.bse.galpy.gasdrag.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag.format.transfer
- petar.mpi.omp.avx2.bse.galpy.gasdrag.hard.debug (107 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag.hard.debug
- petar.mpi.omp.avx2.bse.galpy.hard.debug (96 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.hard.debug
- petar.mpi.omp.avx2.bse.gasdrag.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.gasdrag.format.transfer
- petar.mpi.omp.avx2.bse.gasdrag.hard.debug (102 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.gasdrag.hard.debug
- petar.mpi.omp.avx2.bse.hard.debug (91 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.hard.debug
- petar.mpi.omp.avx2.bse.pnhermite.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.pnhermite.format.transfer
- petar.mpi.omp.avx2.bse.pnhermite.hard.debug (93 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.pnhermite.hard.debug
- petar.mpi.omp.avx2.dsm.galpy.gasdrag.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.galpy.gasdrag.format.transfer
- petar.mpi.omp.avx2.dsm.galpy.gasdrag.hard.debug (93 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.galpy.gasdrag.hard.debug
- petar.mpi.omp.avx2.dsm.gasdrag.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.gasdrag.format.transfer
- petar.mpi.omp.avx2.dsm.gasdrag.hard.debug (88 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.gasdrag.hard.debug
- petar.mpi.omp.avx2.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.format.transfer
- petar.mpi.omp.avx2.galpy.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.format.transfer
- petar.mpi.omp.avx2.galpy.hard.debug (64 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.hard.debug
- petar.mpi.omp.avx2.galpy.mp.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.mp.format.transfer
- petar.mpi.omp.avx2.galpy.mp.hard.debug (64 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.mp.hard.debug
- petar.mpi.omp.avx2.gasdrag.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gasdrag.format.transfer
- petar.mpi.omp.avx2.gasdrag.hard.debug (70 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gasdrag.hard.debug
- petar.mpi.omp.avx2.gpu.bse.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.format.transfer
- petar.mpi.omp.avx2.gpu.bse.galpy.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.format.transfer
- petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag.format.transfer
- petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag.hard.debug (107 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag.hard.debug
- petar.mpi.omp.avx2.gpu.bse.galpy.hard.debug (96 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.galpy.hard.debug
- petar.mpi.omp.avx2.gpu.bse.hard.debug (91 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.bse.hard.debug
- petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag.format.transfer
- petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag.hard.debug (93 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag.hard.debug
- petar.mpi.omp.avx2.gpu.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.format.transfer
- petar.mpi.omp.avx2.gpu.hard.debug (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.gpu.hard.debug
- petar.mpi.omp.avx2.hard.debug (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.hard.debug
- petar.mpi.omp.avx2.kdk.g.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.kdk.g.format.transfer
- petar.mpi.omp.avx2.kdk.g.hard.debug (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.kdk.g.hard.debug
- petar.mpi.omp.avx2.kdkdk4.g.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.kdkdk4.g.format.transfer
- petar.mpi.omp.avx2.kdkdk4.g.hard.debug (59 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.kdkdk4.g.hard.debug
- petar.mpi.omp.avx2.pnall.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.pnall.format.transfer
- petar.mpi.omp.avx2.pnall.hard.debug (62 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.pnall.hard.debug

## Script Tools (Non-solver Workflow Commands)

- petar.data.clear (13 options)
  - path: /home/lwang/bin/petar.data.clear
- petar.data.gether (5 options)
  - path: /home/lwang/bin/petar.data.gether
- petar.data.process (44 options)
  - path: /home/lwang/bin/petar.data.process
- petar.external.pot.movie (12 options)
  - path: /home/lwang/bin/petar.external.pot.movie
- petar.find.dt (13 options)
  - path: /home/lwang/bin/petar.find.dt
- petar.format.transfer.post (16 options)
  - path: /home/lwang/bin/petar.format.transfer.post
- petar.galev.process (20 options)
  - path: /home/lwang/bin/petar.galev.process
- petar.get.init.binary (3 options)
  - path: /home/lwang/bin/petar.get.init.binary
- petar.get.object.snap (26 options)
  - path: /home/lwang/bin/petar.get.object.snap
- petar.init (16 options)
  - path: /home/lwang/bin/petar.init
- petar.movie (87 options)
  - path: /home/lwang/bin/petar.movie
- petar.update.par (9 options)
  - path: /home/lwang/bin/petar.update.par

## Other Petar Commands In PATH

- petar.bse (44 options)
  - path: /home/lwang/bin/petar.bse
- petar.external.agama (7 options)
  - path: /home/lwang/bin/petar.external.agama
- petar.external.galpy (9 options)
  - path: /home/lwang/bin/petar.external.galpy
- petar.galpy.help (5 options)
  - path: /home/lwang/bin/petar.galpy.help
- petar.read.par (4 options)
  - path: /home/lwang/bin/petar.read.par
- petar.select (7 options)
  - path: /home/lwang/bin/petar.select

## Detected Configure-Feature Suffix Tokens

- .64b
  - group: architecture
  - binaries: 1
- .agama
  - group: external-potential
  - binaries: 2
- .avx2
  - group: architecture
  - binaries: 22
- .bse
  - group: interrupt-mode
  - binaries: 9
- .dsm
  - group: interrupt-mode
  - binaries: 3
- .g
  - group: debug
  - binaries: 2
- .galpy
  - group: external-potential
  - binaries: 8
- .gasdrag
  - group: external-hard
  - binaries: 7
- .gpu
  - group: parallel-runtime
  - binaries: 5
- .kdk
  - group: step-mode
  - binaries: 1
- .kdkdk4
  - group: step-mode
  - binaries: 1
- .mp
  - group: mpfrc
  - binaries: 1
- .mpi
  - group: parallel-runtime
  - binaries: 22
- .omp
  - group: parallel-runtime
  - binaries: 22
- .pnall
  - group: post-newtonian
  - binaries: 1
- .pnhermite
  - group: post-newtonian
  - binaries: 1

## Common Core Options

- --ar-ds-scale
- --ar-max-error
- --ar-max-nstep
- --ar-slowdown-factor
- --ar-sym-order
- --center-id
- --disable-print-info
- --domain-nstep
- --domain-weight-mode
- --dt-soft-kepler-nstep
- --dt-soft-sigma-factor
- --energy-err-hard
- --help
- --hermite-acc-offset-sq
- --hermite-de-crit
- --hermite-dm-crit
- --hermite-dt-max
- --hermite-dt-min-index
- --hermite-eta
- --hermite-eta-init
- --hermite-n-neighbor-max
- --id-offset
- --kdtree-n-particles-min
- --n-sample-average
- --r-escape
- --r-group
- --r-ratio
- --r-search-group
- --r-search-min
- --r-search-peri-factor
- --r-search-vel-factor
- --record-id-end-one
- --record-id-end-two
- --record-id-start-one
- --record-id-start-two
- --soft-eps
- --tree-ngroup-limit
- --tree-nleaf-limit
- --tree-nstep-mklist
- --tt-nstep
- --tt-switch
- --write-group-info
- -1
- -6
- -G
- -T
- -a
- -b
- -f
- -h
- -i
- -n
- -o
- -p
- -r
- -s
- -t
- -u
- -w

## Feature-specific Options By Binary

### petar.mpi.omp.avx2.64b

- (no extra options beyond core)

### petar.mpi.omp.avx2.agama

- --agama-conf-file
- --agama-rscale
- --agama-vscale

### petar.mpi.omp.avx2.bse.agama

- --agama-conf-file
- --agama-rscale
- --agama-vscale
- --bse-alpha
- --bse-beta
- --bse-bhflag
- --bse-bhwacc
- --bse-bwind
- --bse-ceflag
- --bse-ecflag
- --bse-eddfac
- --bse-epsnov
- --bse-gamma
- --bse-hewind
- --bse-kmech
- --bse-lambda
- --bse-metallicity
- --bse-mscale
- --bse-neta
- --bse-nsflag
- --bse-psflag
- --bse-pts1
- --bse-pts2
- --bse-pts3
- --bse-rscale
- --bse-sigma
- --bse-tflag
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.bse.galpy.gasdrag

- --bse-alpha
- --bse-beta
- --bse-bhflag
- --bse-bhwacc
- --bse-bwind
- --bse-ceflag
- --bse-ecflag
- --bse-eddfac
- --bse-epsnov
- --bse-gamma
- --bse-hewind
- --bse-kmech
- --bse-lambda
- --bse-metallicity
- --bse-mscale
- --bse-neta
- --bse-nsflag
- --bse-psflag
- --bse-pts1
- --bse-pts2
- --bse-pts3
- --bse-rscale
- --bse-sigma
- --bse-tflag
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --ext-hard-switch
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
- --gdf-K
- --gdf-coulomb-log
- --gdf-gamma
- --gdf-gaspot-index
- --gdf-hard-mode
- --gdf-ifunc-mach-lower
- --gdf-ifunc-mach-upper
- --gdf-ifunc-smooth-order
- --gdf-scale-density
- --gdf-sound-speed
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.bse.galpy

- --bse-alpha
- --bse-beta
- --bse-bhflag
- --bse-bhwacc
- --bse-bwind
- --bse-ceflag
- --bse-ecflag
- --bse-eddfac
- --bse-epsnov
- --bse-gamma
- --bse-hewind
- --bse-kmech
- --bse-lambda
- --bse-metallicity
- --bse-mscale
- --bse-neta
- --bse-nsflag
- --bse-psflag
- --bse-pts1
- --bse-pts2
- --bse-pts3
- --bse-rscale
- --bse-sigma
- --bse-tflag
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.bse.gasdrag

- --bse-alpha
- --bse-beta
- --bse-bhflag
- --bse-bhwacc
- --bse-bwind
- --bse-ceflag
- --bse-ecflag
- --bse-eddfac
- --bse-epsnov
- --bse-gamma
- --bse-hewind
- --bse-kmech
- --bse-lambda
- --bse-metallicity
- --bse-mscale
- --bse-neta
- --bse-nsflag
- --bse-psflag
- --bse-pts1
- --bse-pts2
- --bse-pts3
- --bse-rscale
- --bse-sigma
- --bse-tflag
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --ext-hard-switch
- --gdf-K
- --gdf-coulomb-log
- --gdf-decay-time
- --gdf-gamma
- --gdf-gas-density
- --gdf-hard-mode
- --gdf-ifunc-mach-lower
- --gdf-ifunc-mach-upper
- --gdf-ifunc-smooth-order
- --gdf-sound-speed
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.bse

- --bse-alpha
- --bse-beta
- --bse-bhflag
- --bse-bhwacc
- --bse-bwind
- --bse-ceflag
- --bse-ecflag
- --bse-eddfac
- --bse-epsnov
- --bse-gamma
- --bse-hewind
- --bse-kmech
- --bse-lambda
- --bse-metallicity
- --bse-mscale
- --bse-neta
- --bse-nsflag
- --bse-psflag
- --bse-pts1
- --bse-pts2
- --bse-pts3
- --bse-rscale
- --bse-sigma
- --bse-tflag
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.bse.pnhermite

- --bse-alpha
- --bse-beta
- --bse-bhflag
- --bse-bhwacc
- --bse-bwind
- --bse-ceflag
- --bse-ecflag
- --bse-eddfac
- --bse-epsnov
- --bse-gamma
- --bse-hewind
- --bse-kmech
- --bse-lambda
- --bse-metallicity
- --bse-mscale
- --bse-neta
- --bse-nsflag
- --bse-psflag
- --bse-pts1
- --bse-pts2
- --bse-pts3
- --bse-rscale
- --bse-sigma
- --bse-tflag
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --pn-c
- --pn-crit-hermite
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.dsm.galpy.gasdrag

- --detect-interrupt
- --dsm-dt-factor
- --dsm-epsilon-bh
- --dsm-epsilon-he
- --dsm-epsilon-mdot
- --dsm-he-disk
- --dsm-lambda0
- --dsm-medd
- --dsm-merger-dm
- --dsm-merger-tdelay
- --dsm-new-star-mode
- --dsm-rstar-power
- --dsm-rstar-scale
- --dsm-salpeter-time
- --dsm-seed-mass
- --dsm-speed-of-light
- --ext-hard-switch
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
- --gdf-K
- --gdf-coulomb-log
- --gdf-gamma
- --gdf-gaspot-index
- --gdf-hard-mode
- --gdf-ifunc-mach-lower
- --gdf-ifunc-mach-upper
- --gdf-ifunc-smooth-order
- --gdf-scale-density
- --gdf-sound-speed
- --rand-seed
- --rand-seedfile

### petar.mpi.omp.avx2.dsm.gasdrag

- --detect-interrupt
- --dsm-dt-factor
- --dsm-epsilon-bh
- --dsm-epsilon-he
- --dsm-epsilon-mdot
- --dsm-he-disk
- --dsm-lambda0
- --dsm-medd
- --dsm-merger-dm
- --dsm-merger-tdelay
- --dsm-new-star-mode
- --dsm-rstar-power
- --dsm-rstar-scale
- --dsm-salpeter-time
- --dsm-seed-mass
- --dsm-speed-of-light
- --ext-hard-switch
- --gdf-K
- --gdf-coulomb-log
- --gdf-decay-time
- --gdf-gamma
- --gdf-gas-density
- --gdf-hard-mode
- --gdf-ifunc-mach-lower
- --gdf-ifunc-mach-upper
- --gdf-ifunc-smooth-order
- --gdf-sound-speed
- --rand-seed
- --rand-seedfile

### petar.mpi.omp.avx2.galpy.mp

- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale

### petar.mpi.omp.avx2.galpy

- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale

### petar.mpi.omp.avx2.gasdrag

- --ext-hard-switch
- --gdf-K
- --gdf-coulomb-log
- --gdf-decay-time
- --gdf-gamma
- --gdf-gas-density
- --gdf-hard-mode
- --gdf-ifunc-mach-lower
- --gdf-ifunc-mach-upper
- --gdf-ifunc-smooth-order
- --gdf-sound-speed

### petar.mpi.omp.avx2.gpu.bse.galpy.gasdrag

- --bse-alpha
- --bse-beta
- --bse-bhflag
- --bse-bhwacc
- --bse-bwind
- --bse-ceflag
- --bse-ecflag
- --bse-eddfac
- --bse-epsnov
- --bse-gamma
- --bse-hewind
- --bse-kmech
- --bse-lambda
- --bse-metallicity
- --bse-mscale
- --bse-neta
- --bse-nsflag
- --bse-psflag
- --bse-pts1
- --bse-pts2
- --bse-pts3
- --bse-rscale
- --bse-sigma
- --bse-tflag
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --ext-hard-switch
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
- --gdf-K
- --gdf-coulomb-log
- --gdf-gamma
- --gdf-gaspot-index
- --gdf-hard-mode
- --gdf-ifunc-mach-lower
- --gdf-ifunc-mach-upper
- --gdf-ifunc-smooth-order
- --gdf-scale-density
- --gdf-sound-speed
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.gpu.bse.galpy

- --bse-alpha
- --bse-beta
- --bse-bhflag
- --bse-bhwacc
- --bse-bwind
- --bse-ceflag
- --bse-ecflag
- --bse-eddfac
- --bse-epsnov
- --bse-gamma
- --bse-hewind
- --bse-kmech
- --bse-lambda
- --bse-metallicity
- --bse-mscale
- --bse-neta
- --bse-nsflag
- --bse-psflag
- --bse-pts1
- --bse-pts2
- --bse-pts3
- --bse-rscale
- --bse-sigma
- --bse-tflag
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.gpu.bse

- --bse-alpha
- --bse-beta
- --bse-bhflag
- --bse-bhwacc
- --bse-bwind
- --bse-ceflag
- --bse-ecflag
- --bse-eddfac
- --bse-epsnov
- --bse-gamma
- --bse-hewind
- --bse-kmech
- --bse-lambda
- --bse-metallicity
- --bse-mscale
- --bse-neta
- --bse-nsflag
- --bse-psflag
- --bse-pts1
- --bse-pts2
- --bse-pts3
- --bse-rscale
- --bse-sigma
- --bse-tflag
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.gpu.dsm.galpy.gasdrag

- --detect-interrupt
- --dsm-dt-factor
- --dsm-epsilon-bh
- --dsm-epsilon-he
- --dsm-epsilon-mdot
- --dsm-he-disk
- --dsm-lambda0
- --dsm-medd
- --dsm-merger-dm
- --dsm-merger-tdelay
- --dsm-new-star-mode
- --dsm-rstar-power
- --dsm-rstar-scale
- --dsm-salpeter-time
- --dsm-seed-mass
- --dsm-speed-of-light
- --ext-hard-switch
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
- --gdf-K
- --gdf-coulomb-log
- --gdf-gamma
- --gdf-gaspot-index
- --gdf-hard-mode
- --gdf-ifunc-mach-lower
- --gdf-ifunc-mach-upper
- --gdf-ifunc-smooth-order
- --gdf-scale-density
- --gdf-sound-speed
- --rand-seed
- --rand-seedfile

### petar.mpi.omp.avx2.gpu

- (no extra options beyond core)

### petar.mpi.omp.avx2.kdk.g

- (no extra options beyond core)

### petar.mpi.omp.avx2.kdkdk4.g

- (no extra options beyond core)

### petar.mpi.omp.avx2

- (no extra options beyond core)

### petar.mpi.omp.avx2.pnall

- --pn-c
- --pn-crit-ar
- --pn-crit-hermite


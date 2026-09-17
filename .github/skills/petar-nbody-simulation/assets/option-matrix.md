# PeTar Option Matrix (Auto-generated)

This file is generated from commands discovered in PATH via <binary> -h option extraction.

Discovery rule:
- scan PATH for executable names starting with 'petar'
- classify helper tools by suffix (*.hard.debug, *.format.transfer)
- classify solver binaries by required core options (-u, -t, -o)

## Solver Binaries Included

- petar.mpi.omp.avx2.agama (63 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.agama
- petar.mpi.omp.avx2.bse.agama (95 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama
- petar.mpi.omp.avx2.bse.galpy.gasdrag (108 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.bse.galpy (97 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- petar.mpi.omp.avx2.bse (92 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse
- petar.mpi.omp.avx2.bseEmp.galpy (99 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bseEmp.galpy
- petar.mpi.omp.avx2.btlogh.bseEmp.galpy (100 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.bseEmp.galpy
- petar.mpi.omp.avx2.btlogh (62 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh
- petar.mpi.omp.avx2.dsm (78 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm
- petar.mpi.omp.avx2.galpy.gasdrag (77 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.gasdrag
- petar.mpi.omp.avx2.galpy (66 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy
- petar.mpi.omp.avx2.merger (63 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.merger
- petar.mpi.omp.avx2.mp (60 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.mp
- petar.mpi.omp.avx2 (62 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2
- petar.omp.avx2.64b (58 options)
  - path: /home/lwang/bin/petar.omp.avx2.64b
- petar.omp.avx2 (58 options)
  - path: /home/lwang/bin/petar.omp.avx2

## Requirement-Driven Solver Filters

No fixed priority across physics scenarios. Select filters by user requirements first, then choose by performance.

Configure mapping:
- --with-interrupt -> suffix tokens: merger | bse | mobse | bseEmp | dsm
- --with-external -> suffix tokens: galpy | agama
- --with-external-hard -> suffix token: gasdrag
- --with-pn -> suffix tokens: pn*

### Interrupt Module (--with-interrupt)

- petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 108
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.btlogh.bseEmp.galpy
  - perf-tags: avx2,omp,mpi
  - options: 100
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.bseEmp.galpy
- petar.mpi.omp.avx2.bseEmp.galpy
  - perf-tags: avx2,omp,mpi
  - options: 99
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bseEmp.galpy
- petar.mpi.omp.avx2.bse.galpy
  - perf-tags: avx2,omp,mpi
  - options: 97
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- petar.mpi.omp.avx2.bse.agama
  - perf-tags: avx2,omp,mpi
  - options: 95
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama
- petar.mpi.omp.avx2.bse
  - perf-tags: avx2,omp,mpi
  - options: 92
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse
- petar.mpi.omp.avx2.dsm
  - perf-tags: avx2,omp,mpi
  - options: 78
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm
- petar.mpi.omp.avx2.merger
  - perf-tags: avx2,omp,mpi
  - options: 63
  - path: /home/lwang/bin/petar.mpi.omp.avx2.merger

### Stellar Evolution (bse/mobse/bseEmp)

- petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 108
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.btlogh.bseEmp.galpy
  - perf-tags: avx2,omp,mpi
  - options: 100
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.bseEmp.galpy
- petar.mpi.omp.avx2.bseEmp.galpy
  - perf-tags: avx2,omp,mpi
  - options: 99
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bseEmp.galpy
- petar.mpi.omp.avx2.bse.galpy
  - perf-tags: avx2,omp,mpi
  - options: 97
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- petar.mpi.omp.avx2.bse.agama
  - perf-tags: avx2,omp,mpi
  - options: 95
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama
- petar.mpi.omp.avx2.bse
  - perf-tags: avx2,omp,mpi
  - options: 92
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse

### Long-timescale External: Galpy (--with-external=galpy)

- petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 108
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.btlogh.bseEmp.galpy
  - perf-tags: avx2,omp,mpi
  - options: 100
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.bseEmp.galpy
- petar.mpi.omp.avx2.bseEmp.galpy
  - perf-tags: avx2,omp,mpi
  - options: 99
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bseEmp.galpy
- petar.mpi.omp.avx2.bse.galpy
  - perf-tags: avx2,omp,mpi
  - options: 97
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- petar.mpi.omp.avx2.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 77
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.gasdrag
- petar.mpi.omp.avx2.galpy
  - perf-tags: avx2,omp,mpi
  - options: 66
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy

### Long-timescale External: Agama (--with-external=agama)

- petar.mpi.omp.avx2.bse.agama
  - perf-tags: avx2,omp,mpi
  - options: 95
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama
- petar.mpi.omp.avx2.agama
  - perf-tags: avx2,omp,mpi
  - options: 63
  - path: /home/lwang/bin/petar.mpi.omp.avx2.agama

### Short-timescale External Hard: Gas Drag (--with-external-hard=gasdrag)

- petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 108
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- petar.mpi.omp.avx2.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 77
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.gasdrag

### Post-Newtonian Relativity (--with-pn)

- (no matching solver detected)

### Special-purpose Features (Non-default Recommendations)

These are for special user requests and should not be selected by default.

### High-precision Tree Force (--enable-64b)

- petar.omp.avx2.64b
  - perf-tags: avx2,omp
  - options: 58
  - path: /home/lwang/bin/petar.omp.avx2.64b

### High-precision Position Representation (--enable-mpfrc)

- petar.mpi.omp.avx2.mp
  - perf-tags: avx2,omp,mpi
  - options: 60
  - path: /home/lwang/bin/petar.mpi.omp.avx2.mp

### Debug Build: g (--with-debug=g)

- (no matching solver detected)

### Debug Build: assert (--with-debug=assert)

- (no matching solver detected)


## Recommended Solver Candidates (Performance-first)

Priority order: gpu > avx512 > avx2 > omp > mpi
Apply this ordering after choosing a requirement-driven filter.
- #1 petar.mpi.omp.avx2.bse.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 108
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag
- #2 petar.mpi.omp.avx2.btlogh.bseEmp.galpy
  - perf-tags: avx2,omp,mpi
  - options: 100
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.bseEmp.galpy
- #3 petar.mpi.omp.avx2.bseEmp.galpy
  - perf-tags: avx2,omp,mpi
  - options: 99
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bseEmp.galpy
- #4 petar.mpi.omp.avx2.bse.galpy
  - perf-tags: avx2,omp,mpi
  - options: 97
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy
- #5 petar.mpi.omp.avx2.bse.agama
  - perf-tags: avx2,omp,mpi
  - options: 95
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama
- #6 petar.mpi.omp.avx2.bse
  - perf-tags: avx2,omp,mpi
  - options: 92
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse
- #7 petar.mpi.omp.avx2.dsm
  - perf-tags: avx2,omp,mpi
  - options: 78
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm
- #8 petar.mpi.omp.avx2.galpy.gasdrag
  - perf-tags: avx2,omp,mpi
  - options: 77
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.gasdrag
- #9 petar.mpi.omp.avx2.galpy
  - perf-tags: avx2,omp,mpi
  - options: 66
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy
- #10 petar.mpi.omp.avx2.merger
  - perf-tags: avx2,omp,mpi
  - options: 63
  - path: /home/lwang/bin/petar.mpi.omp.avx2.merger
- #11 petar.mpi.omp.avx2.agama
  - perf-tags: avx2,omp,mpi
  - options: 63
  - path: /home/lwang/bin/petar.mpi.omp.avx2.agama
- #12 petar.mpi.omp.avx2.btlogh
  - perf-tags: avx2,omp,mpi
  - options: 62
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh

## Helper Tool Binaries Excluded From Simulation Commands

- petar.mpi.omp.avx2.agama.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.agama.format.transfer
- petar.mpi.omp.avx2.agama.hard.debug (63 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.agama.hard.debug
- petar.mpi.omp.avx2.bse.agama.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama.format.transfer
- petar.mpi.omp.avx2.bse.agama.hard.debug (95 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.agama.hard.debug
- petar.mpi.omp.avx2.bse.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.format.transfer
- petar.mpi.omp.avx2.bse.galpy.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.format.transfer
- petar.mpi.omp.avx2.bse.galpy.gasdrag.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag.format.transfer
- petar.mpi.omp.avx2.bse.galpy.gasdrag.hard.debug (108 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag.hard.debug
- petar.mpi.omp.avx2.bse.galpy.hard.debug (97 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.hard.debug
- petar.mpi.omp.avx2.bse.hard.debug (92 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.hard.debug
- petar.mpi.omp.avx2.bseEmp.galpy.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bseEmp.galpy.format.transfer
- petar.mpi.omp.avx2.bseEmp.galpy.hard.debug (99 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bseEmp.galpy.hard.debug
- petar.mpi.omp.avx2.btlogh.bseEmp.galpy.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.bseEmp.galpy.format.transfer
- petar.mpi.omp.avx2.btlogh.bseEmp.galpy.hard.debug (100 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.bseEmp.galpy.hard.debug
- petar.mpi.omp.avx2.btlogh.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.format.transfer
- petar.mpi.omp.avx2.btlogh.hard.debug (62 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.hard.debug
- petar.mpi.omp.avx2.dsm.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.format.transfer
- petar.mpi.omp.avx2.dsm.hard.debug (78 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dsm.hard.debug
- petar.mpi.omp.avx2.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.format.transfer
- petar.mpi.omp.avx2.galpy.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.format.transfer
- petar.mpi.omp.avx2.galpy.gasdrag.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.gasdrag.format.transfer
- petar.mpi.omp.avx2.galpy.gasdrag.hard.debug (77 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.gasdrag.hard.debug
- petar.mpi.omp.avx2.galpy.hard.debug (66 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.hard.debug
- petar.mpi.omp.avx2.hard.debug (62 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.hard.debug
- petar.mpi.omp.avx2.merger.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.merger.format.transfer
- petar.mpi.omp.avx2.merger.hard.debug (63 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.merger.hard.debug
- petar.mpi.omp.avx2.mp.format.transfer (9 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.mp.format.transfer
- petar.mpi.omp.avx2.mp.hard.debug (60 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.mp.hard.debug
- petar.omp.avx2.64b.format.transfer (9 options)
  - path: /home/lwang/bin/petar.omp.avx2.64b.format.transfer
- petar.omp.avx2.64b.hard.debug (58 options)
  - path: /home/lwang/bin/petar.omp.avx2.64b.hard.debug
- petar.omp.avx2.format.transfer (9 options)
  - path: /home/lwang/bin/petar.omp.avx2.format.transfer
- petar.omp.avx2.hard.debug (58 options)
  - path: /home/lwang/bin/petar.omp.avx2.hard.debug

## Script Tools (Non-solver Workflow Commands)

- petar.data.clear (14 options)
  - path: /home/lwang/bin/petar.data.clear
- petar.data.gether (5 options)
  - path: /home/lwang/bin/petar.data.gether
- petar.data.process (44 options)
  - path: /home/lwang/bin/petar.data.process
- petar.external.pot.movie (12 options)
  - path: /home/lwang/bin/petar.external.pot.movie
- petar.find.dt (12 options)
  - path: /home/lwang/bin/petar.find.dt
- petar.format.transfer.post (16 options)
  - path: /home/lwang/bin/petar.format.transfer.post
- petar.galev.process (20 options)
  - path: /home/lwang/bin/petar.galev.process
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
- petar.bseEmp (45 options)
  - path: /home/lwang/bin/petar.bseEmp
- petar.external.agama (7 options)
  - path: /home/lwang/bin/petar.external.agama
- petar.external.galpy (9 options)
  - path: /home/lwang/bin/petar.external.galpy
- petar.galpy.help (5 options)
  - path: /home/lwang/bin/petar.galpy.help
- petar.get.init.binary.bse (3 options)
  - path: /home/lwang/bin/petar.get.init.binary.bse
- petar.hard.test (40 options)
  - path: /home/lwang/bin/petar.hard.test
- petar.mpi.omp.avx2.bse.galpy.gasdrag.dump2test (39 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bse.galpy.gasdrag.dump2test
- petar.mpi.omp.avx2.bseEmp.galpy.dump2test (42 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.bseEmp.galpy.dump2test
- petar.mpi.omp.avx2.btlogh.bseEmp.galpy.dump2test (43 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.bseEmp.galpy.dump2test
- petar.mpi.omp.avx2.btlogh.dump2test (41 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.btlogh.dump2test
- petar.mpi.omp.avx2.dump2test (44 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.dump2test
- petar.mpi.omp.avx2.galpy.dump2test (40 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.dump2test
- petar.mpi.omp.avx2.galpy.gasdrag.dump2test (40 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.galpy.gasdrag.dump2test
- petar.mpi.omp.avx2.merger.dump2test (42 options)
  - path: /home/lwang/bin/petar.mpi.omp.avx2.merger.dump2test
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
  - binaries: 16
- .bse
  - group: interrupt-mode
  - binaries: 4
- .bseEmp
  - group: interrupt-mode
  - binaries: 2
- .dsm
  - group: interrupt-mode
  - binaries: 1
- .galpy
  - group: external-potential
  - binaries: 6
- .gasdrag
  - group: external-hard
  - binaries: 2
- .merger
  - group: interrupt-mode
  - binaries: 1
- .mp
  - group: mpfrc
  - binaries: 1
- .mpi
  - group: parallel-runtime
  - binaries: 14
- .omp
  - group: parallel-runtime
  - binaries: 16

## Common Core Options

- --ar-ds-scale
- --ar-max-error
- --ar-max-nstep
- --ar-slowdown-factor
- --ar-sym-order
- --center-id
- --disable-print-info
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
- --keep-tmp-on-startup
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

### petar.mpi.omp.avx2.agama

- --agama-conf-file
- --agama-rscale
- --agama-vscale
- --domain-nstep
- --domain-weight-mode

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
- --domain-nstep
- --domain-weight-mode
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
- --domain-nstep
- --domain-weight-mode
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
- --domain-nstep
- --domain-weight-mode
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
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
- --domain-nstep
- --domain-weight-mode
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.bseEmp.galpy

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
- --bse-trackmode
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --domain-nstep
- --domain-weight-mode
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
- --r-search-group-safety
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.btlogh.bseEmp.galpy

- --ar-g-func
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
- --bse-trackmode
- --bse-tscale
- --bse-vscale
- --bse-wdflag
- --bse-xi
- --detect-interrupt
- --domain-nstep
- --domain-weight-mode
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
- --r-search-group-safety
- --rand-seed
- --rand-seedfile
- --stellar-evolution

### petar.mpi.omp.avx2.btlogh

- --ar-g-func
- --domain-nstep
- --domain-weight-mode
- --r-search-group-safety

### petar.mpi.omp.avx2.dsm

- --detect-interrupt
- --domain-nstep
- --domain-weight-mode
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
- --rand-seed
- --rand-seedfile

### petar.mpi.omp.avx2.galpy.gasdrag

- --domain-nstep
- --domain-weight-mode
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
- --r-search-group-safety

### petar.mpi.omp.avx2.galpy

- --domain-nstep
- --domain-weight-mode
- --galpy-conf-file
- --galpy-rscale
- --galpy-set
- --galpy-type-arg
- --galpy-vscale
- --r-search-group-safety

### petar.mpi.omp.avx2.merger

- --ar-g-func
- --detect-interrupt
- --domain-nstep
- --domain-weight-mode
- --r-search-group-safety

### petar.mpi.omp.avx2.mp

- --domain-nstep
- --domain-weight-mode

### petar.mpi.omp.avx2

- --ar-g-func
- --domain-nstep
- --domain-weight-mode
- --r-search-group-safety

### petar.omp.avx2.64b

- (no extra options beyond core)

### petar.omp.avx2

- (no extra options beyond core)


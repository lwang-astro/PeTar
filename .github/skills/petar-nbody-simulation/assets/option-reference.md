# PeTar CLI Option Reference (Source-Generated)

All possible PeTar command-line options across every configure variant.
Generated from source headers — no binary compilation needed.

> Generated: 2026-06-09 01:37 UTC
> Regenerate: `python3 .github/skills/petar-nbody-simulation/assets/generate_option_reference.py`

See [`option-matrix.md`](option-matrix.md) for currently installed binaries.
Use `<binary> -h` for exact runtime option validation.

---

## Summary

| Configure Flag Group | Option Count |
|-----------------------|-------------|
| Core (always available) | 76 |
| BSE Interrupt (`--with-interrupt=bse|mobse|bseEmp`) | 61 |
| DSM Interrupt (`--with-interrupt=dsm`) | 17 |
| Galpy External (`--with-external=galpy`) | 10 |
| Agama External (`--with-external=agama`) | 3 |
| Gas Drag (`--with-external-hard=gasdrag`) | 10 |
| External Hard (`--with-external-hard`) | 2 |
| **Total** | **179** |

---

## Core (always available)

### `--G`

- **Type**: `PS::F64`
- **Default**: `1.0`
- **Description**: Gravitational constant, if -u 1, G = 0.00449830997959438 pc^3/(Msun*Myr^2)
- **Source**: `src/petar.hpp`:155

### `--T`

- **Type**: `PS::F64`
- **Default**: `0.3`
- **Description**: Particle-tree opening angle theta
- **Source**: `src/petar.hpp`:143

### `--a`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Data file output mode; 0: overwrite files except object dump files, include header lines; 1: append files except snapshots, no header line
- **Source**: `src/petar.hpp`:173

### `--ar-ds-scale`

- **Type**: `PS::F64`
- **Default**: `1.0`
- **Description**: Scale factor for SDAR step size calculation
- **Source**: `src/hard.hpp`:140

### `--ar-max-error`

- **Type**: `PS::F64`
- **Default**: `1e-8`
- **Description**: Maximum energy error allowed for the SDAR integrator
- **Source**: `src/hard.hpp`:137

### `--ar-max-nstep`

- **Type**: `PS::S64`
- **Default**: `1000000`
- **Description**: Maximum step allowed for the SDAR sym integrator
- **Source**: `src/hard.hpp`:138

### `--ar-slowdown-factor`

- **Type**: `PS::F64`
- **Default**: `1e-4`
- **Description**: Slowdown perturbation criterion
- **Source**: `src/hard.hpp`:141

### `--ar-sym-order`

- **Type**: `PS::S64`
- **Default**: `-6`
- **Description**: Order of the symplectic integrator for SDAR, should be even number; -6,-8: Yoshida 2nd symplectic method; 4,6,8,...: Yoshida 1st symplectic method
- **Source**: `src/hard.hpp`:139

### `--b`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Number of primordial binaries (n_bin) for initialization (assuming the binaries' IDs are 1,2*n_bin)
- **Source**: `src/petar.hpp`:152

### `--center-id`

- **Type**: `PS::S64`
- **Default**: `-1`
- **Description**: id of the central object for a system like a stellar disk
- **Source**: `src/hard.hpp`:132

### `--detect-interrupt`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Stellar evolution of binaries in SDAR integration; 0: switch off; 1: using BSE based code (if '--stellar-evolution != 0)
- **Source**: `src/hard.hpp`:147
- **Guards**: `BSE_BASE` → --with-interrupt=base, `STELLAR_EVOLUTION` → --with-interrupt

### `--domain-nstep`

- **Type**: `PS::S64`
- **Default**: `16`
- **Description**: Number of steps between domain decompositions
- **Source**: `src/petar.hpp`:177
- **Guards**: `PARTICLE_SIMULATOR_MPI_PARALLEL` → MPI build

### `--domain-weight-mode`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Domain decomposition weight mode for MPI parallel; 0: equal weight for each MPI processor; 1: use force calculation time as weight to obtain better load balance with losing simulation reproducibility
- **Source**: `src/petar.hpp`:176
- **Guards**: `PARTICLE_SIMULATOR_MPI_PARALLEL` → MPI build

### `--dt-soft-kepler-nstep`

- **Type**: `PS::F64`
- **Default**: `16.0`
- **Description**: Factor 'nstep' to determine dt_soft by P(r_in)/nstep, see option '-s' and '-r'
- **Source**: `src/petar.hpp`:161
- **Guards**: `KDKDK_4TH` → 4th-order KDKDK step mode

### `--dt-soft-sigma-factor`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Factor 'alpha' to determine dt_soft by alpha*r_in/sigma_3D, see option '-s' and '-r'; = 0: not used, apply --dt-soft-kepler-nstep; > 0: use this option instead of '--dt-soft-kepler-nstep'
- **Source**: `src/petar.hpp`:165

### `--energy-err-hard`

- **Type**: `PS::F64`
- **Default**: `1e-4`
- **Description**: Maximum energy error allowed for the hard integrator
- **Source**: `src/hard.hpp`:117
- **Guards**: `HARD_CHECK_ENERGY` → development debug

### `--f`

- **Type**: `std::string`
- **Default**: `"data"`
- **Description**: Prefix of filenames for output data: [prefix].**
- **Source**: `src/petar.hpp`:179

### `--gdf-decay-time`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: gas density decay time scale in units of PeTar input, if 0, no decay
- **Source**: `src/gas_drag.hpp`:57
- **Guards**: `GALPY` → --with-external=galpy (negated)

### `--gdf-gas-density`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: gas density in units of PeTar input
- **Source**: `src/gas_drag.hpp`:56
- **Guards**: `GALPY` → --with-external=galpy (negated)

### `--gdf-gaspot-index`

- **Type**: `long long int`
- **Default**: `-1`
- **Description**: galpy potential set index for gas component, used for obtaining gas density
- **Source**: `src/gas_drag.hpp`:53
- **Guards**: `GALPY` → --with-external=galpy

### `--gdf-scale-density`

- **Type**: `double`
- **Default**: `1/G_ASTRO`
- **Description**: scale factor for galpy potential density
- **Source**: `src/gas_drag.hpp`:54
- **Guards**: `GALPY` → --with-external=galpy

### `--hermite-acc-offset-sq`

- **Type**: `PS::F64`
- **Default**: `-1.0`
- **Description**: Square acceleration offset for Hermite time step calculation to avoid too small step when weak acceleration exists; = -1: calculate from mean mass <m> and r_out (G*<m>/r_out^2)^2; = 0: no offset; > 0: custom offset value
- **Source**: `src/hard.hpp`:125

### `--hermite-de-crit`

- **Type**: `PS::F64`
- **Default**: `1e-4`
- **Description**: Ekin change rate criterion for reinitializing hermite time step
- **Source**: `src/hard.hpp`:143

### `--hermite-dm-crit`

- **Type**: `PS::F64`
- **Default**: `1e-4`
- **Description**: Mass change rate criterion for reinitializing hermite time step
- **Source**: `src/hard.hpp`:142

### `--hermite-dt-max`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Maximum hermite timestep
- **Source**: `src/hard.hpp`:135

### `--hermite-dt-min-index`

- **Type**: `PS::S64`
- **Default**: `40`
- **Description**: Power index n for the smallest timestep (0.5^n) allowed in the Hermite integrator
- **Source**: `src/hard.hpp`:136

### `--hermite-eta`

- **Type**: `PS::F64`
- **Default**: `0.1`
- **Description**: Hermite timestep coefficient eta
- **Source**: `src/hard.hpp`:133

### `--hermite-eta-init`

- **Type**: `PS::F64`
- **Default**: `0.001`
- **Description**: Hermite timestep coefficient eta for initial step in 2nd order
- **Source**: `src/hard.hpp`:134

### `--hermite-n-neighbor-max`

- **Type**: `PS::S64`
- **Default**: `300`
- **Description**: Maximum number of group neighbors to be stored
- **Source**: `src/hard.hpp`:144

### `--i`

- **Type**: `PS::S64`
- **Default**: `2`
- **Description**: Data file reading and writing format; snapshots, status and escaper outputs follow the write mode selected here; 0: read and write in BINARY; 1: read and write in ASCII; 2: read in ASCII, write in BINARY; 3: read in BINARY, write in ASCII
- **Source**: `src/petar.hpp`:171

### `--id-offset`

- **Type**: `PS::S64`
- **Default**: `-1`
- **Description**: Starting ID for artificial particles, total number of real particles must always be smaller than this
- **Source**: `src/hard.hpp`:131

### `--kdtree-n-particles-min`

- **Type**: `PS::S64`
- **Default**: `32`
- **Description**: Minimum number of particles + groups for building kdtree to speed up neighbor search in Hermite-only neighbor force calculation
- **Source**: `src/hard.hpp`:170
- **Guards**: `HERMITE_ONLY_CALC_NEIGHBOR_FORCE` → --enable-hermite-only-calc-neighbor-force

### `--keep-tmp-on-startup`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Startup tmp handling for transactional outputs; 0: auto-remove residual tmp files before run (default); 1: keep residual tmp files and only print warning
- **Source**: `src/petar.hpp`:174

### `--n`

- **Type**: `PS::S64`
- **Default**: `100000`
- **Description**: Total number of particles, used only when the input data filename is __Plummer
- **Source**: `src/petar.hpp`:156

### `--n-sample-average`

- **Type**: `PS::S64`
- **Default**: `100`
- **Description**: Average target number of sample particles per process
- **Source**: `src/petar.hpp`:151

### `--o`

- **Type**: `PS::F64`
- **Default**: `1.0`
- **Description**: Output time interval for particle dataset snapshots
- **Source**: `src/petar.hpp`:170

### `--os-nsplit`

- **Type**: `PS::S64`
- **Default**: `4`
- **Description**: Number of binary sample points for tree perturbation force using orbit-sampling method
- **Source**: `src/hard.hpp`:129
- **Guards**: `ORBIT_SAMPLING` → --enable-orbit-sampling

### `--p`

- **Type**: `std::string`
- **Default**: `"input.par"`
- **Description**: Input parameter file (this option should be used first before any other options)
- **Source**: `src/petar.hpp`:180

### `--pn-c`

- **Type**: `PS::F64`
- **Default**: `1`
- **Description**: speed of light value for Post Newtonian; if -u 1 is used, auto determined
- **Source**: `src/hard.hpp`:161

### `--pn-crit-ar`

- **Type**: `PS::F64`
- **Default**: `1e-6`
- **Description**: AR speed criterion to switch on PN terms, min (v/c)^2
- **Source**: `src/hard.hpp`:167
- **Guards**: `SDAR_PN` → --with-pn (SDAR)

### `--pn-crit-hermite`

- **Type**: `PS::F64`
- **Default**: `1e-6`
- **Description**: Hermite speed criterion to switch on PN terms, min (v/c)^2
- **Source**: `src/hard.hpp`:164
- **Guards**: `HERMITE_PN` → --with-pn (Hermite)

### `--r`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Outer changeover radius (r_out); > 0: custom r_out value and check '-s dt_soft';      dt_soft = 0: calculate dt_soft and then adjust r_out by dt_soft;      dt_soft > 0: use custom r_out directly; = 0 (default): check '--dt-soft-sigma-factor alpha':;      alpha > 0: r_out = alpha*dt_soft*sigma_3D/r-ratio;          sigma_3D: global 3D velocity dispersion;          r-ratio: defined by --r-ratio;      alpha = 0 (default): r_out = a(r_in)/r-ratio;          a(r_in): the binary semi-major axis with the period of nstep*dt_soft;          nstep: defined by --dt-soft-kepler-nstep
- **Source**: `src/petar.hpp`:158

### `--r-escape`

- **Type**: `PS::F64`
- **Default**: `PS::LARGE_FLOAT`
- **Description**: Object escape radius criterion; < 0: remove objects when r>-r_escape; >= 0: remove objects when r>r_escape and energy>0
- **Source**: `src/petar.hpp`:169

### `--r-group`

- **Type**: `PS::F64`
- **Default**: `-1.0`
- **Description**: Tidal tensor box size and the radial criterion for detecting multiple groups (binaries, triples, etc.); = -1: auto-determine by 0.8*r_search_group; = 0: switch off SDAR; > 0: custom criterion value
- **Source**: `src/hard.hpp`:121

### `--r-ratio`

- **Type**: `PS::F64`
- **Default**: `0.1`
- **Description**: Ratio between inner (r_in) and outer (r_out) changeover radii
- **Source**: `src/petar.hpp`:159

### `--r-search-group`

- **Type**: `PS::F64`
- **Default**: `-1.0`
- **Description**: The radial criterion for detecting multiple group candidates; = -1: auto-determine by 1.0*r_in; = 0: switch off SDAR; > 0: custom criterion value
- **Source**: `src/hard.hpp`:122

### `--r-search-min`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Minimum neighbor search radius for hard clusters; = 0: auto-determine by max(search-vel-factor*sigma_1D*dt_soft + rout, 1.2 r_out); > 0: custom search radius value
- **Source**: `src/petar.hpp`:168

### `--r-search-peri-factor`

- **Type**: `PS::F64`
- **Default**: `1.5`
- **Description**: Neighbor search coefficient for periapsis check
- **Source**: `src/petar.hpp`:167

### `--r-search-vel-factor`

- **Type**: `PS::F64`
- **Default**: `3.0`
- **Description**: Neighbor search coefficient for velocity check (v*dt)
- **Source**: `src/petar.hpp`:166

### `--rand-seed`

- **Type**: `long long int`
- **Default**: `0`
- **Description**: Random number seed (positive integer); suppressed when --rand-seedfile is provided; if used, the single integer seed is used to generate multiple seeds for each pair of OpenMP thread and MPI processor
- **Source**: `parallel-random/rand_io.hpp`:18

### `--rand-seedfile`

- **Type**: `std::string`
- **Default**: `"__NONE__"`
- **Description**: Name for a file contain random seeds of all threads and MPI processors; For restart the simulation, the randseeds file can be used to restore all seeds
- **Source**: `parallel-random/rand_io.hpp`:19

### `--record-id-end-one`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Ending of the first id range for hard dump; notice that the ending id is not included in hard dump
- **Source**: `src/hard.hpp`:157

### `--record-id-end-two`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Ending of the 2nd id range for hard dump; notice that the ending id is not included in hard dump
- **Source**: `src/hard.hpp`:159

### `--record-id-start-one`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Starting of the first id range for hard dump recording every tree step, save into files object_[id]
- **Source**: `src/hard.hpp`:156

### `--record-id-start-two`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Starting of the 2nd id range for hard dump recording every tree step
- **Source**: `src/hard.hpp`:158

### `--s`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Tree timestep (dt_soft); > 0: custom dt_soft value, regularized to 0.5^n, where n is an integer; = 0: check '-r r_out':;      r_out = 0 (default): dt_soft = 2.6E-4*GM/sigma_3D^3, and is regularized to 0.5^n;          sigma_3D: global 3D velocity dispersion;      r_out > 0: check '--dt-soft-sigma-factor alpha':;          alpha > 0: dt_soft = alpha*r_in/(sqrt(3)*sigma);              r_in: determined by --r-ratio and r_out;          alpha = 0 (default): dt_soft = P(r_in)/nstep;              P(r_in): the binary period with the semi-major axis of r_in;              nstep: defined by --dt-soft-kepler-nstep
- **Source**: `src/petar.hpp`:157

### `--snap-filename`

- **Type**: `std::string`
- **Default**: `"__NONE__"`
- **Description**: Input data file
- **Source**: `src/petar.hpp`:181

### `--soft-eps`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Softening epsilon
- **Source**: `src/hard.hpp`:120

### `--stellar-evolution`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Stellar evolution of stars in Hermite and SDAR integration; 0: switch off; >=1: using SSE/BSE based codes; 2: activate dynamical tide and hyperbolic GW radiation
- **Source**: `src/hard.hpp`:148
- **Guards**: `BSE_BASE` → --with-interrupt=base, `STELLAR_EVOLUTION` → --with-interrupt

### `--t`

- **Type**: `PS::F64`
- **Default**: `10.0`
- **Description**: End time of simulation
- **Source**: `src/petar.hpp`:153

### `--tree-ngroup-limit`

- **Type**: `PS::S64`
- **Default**: `1024`
- **Description**: Particle-tree group number limit; Optimal value for x86-AVX512 is 1024
- **Source**: `src/petar.hpp`:146
- **Guards**: `USE__AVX512`

### `--tree-nleaf-limit`

- **Type**: `PS::S64`
- **Default**: `20`
- **Description**: Particle-tree leaf number limit; Optimal value should be slightly >= artificial particle number (tidal tensor 8 + anti-force sample 3) + 2 (binary member) + 1 (binary c.m.)
- **Source**: `src/petar.hpp`:144

### `--tree-nstep-mklist`

- **Type**: `PS::S64`
- **Default**: `2`
- **Description**: Particle-tree make-list interval in number of soft-step
- **Source**: `src/petar.hpp`:150

### `--tt-nstep`

- **Type**: `PS::S64`
- **Default**: `4`
- **Description**: Number of steps per slow-down binary orbits (period/dt_soft) for isolated binaries; also the maximum criterion for activating tidal tensor method
- **Source**: `src/hard.hpp`:126

### `--tt-switch`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Tidal tensor calculation for (counter-)perturbation (from)on binaries: 0: off, 1: on
- **Source**: `src/hard.hpp`:127

### `--u`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Input data unit; 0: based on the value of G; 1: mass:Msun, length:pc, time:Myr, velocity:pc/Myr, modify G to fit this unit set
- **Source**: `src/petar.hpp`:154

### `--w`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Data file writing style; 0: no output; 1: write all files separately; 2. write snapshots in status files in one line per step (no MPI support); 3. write files except snapshots
- **Source**: `src/petar.hpp`:172

### `--write-group-info`

- **Type**: `PS::S64`
- **Default**: `2`
- **Description**: Write information of new and end groups; 0: no output; 1: ascii output; 2: binary output, files are [data filename prefix].group.[MPI rank].n[N_member]
- **Source**: `src/hard.hpp`:154
- **Guards**: `ADJUST_GROUP_PRINT` → --enable-adjust-group-print

## BSE Interrupt (`--with-interrupt=bse|mobse|bseEmp`)

### `--bse-alpha`

- **Type**: `double`
- **Default**: `3.0`
- **Description**: Common-envelope efficiency parameter
- **Source**: `bse-interface/bse_interface.h`:868

### `--bse-beta`

- **Type**: `double`
- **Default**: `0.125`
- **Description**: wind velocity factor: proportional to vwind**2
- **Source**: `bse-interface/bse_interface.h`:870

### `--bse-bhflag`

- **Type**: `long long int`
- **Default**: `2`
- **Description**: BH kick option: 0: no kick; 1: same as NS; 2: scaled by fallback
- **Source**: `bse-interface/bse_interface.h`:881

### `--bse-bhwacc`

- **Type**: `double`
- **Default**: `1.5`
- **Description**: Bondi-Hoyle wind accretion factor
- **Source**: `bse-interface/bse_interface.h`:872

### `--bse-bwind`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: Binary enhanced mass loss parameter, inactive for single
- **Source**: `bse-interface/bse_interface.h`:865

### `--bse-ceflag`

- **Type**: `long long int`
- **Default**: `0`
- **Description**: if =3, activates de Kool common-envelope model
- **Source**: `bse-interface/bse_interface.h`:877

### `--bse-ecflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, ECS is switched on
- **Source**: `bse-interface/bse_interface.h`:885

### `--bse-eddfac`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Eddington limit factor for mass transfer
- **Source**: `bse-interface/bse_interface.h`:874

### `--bse-epsnov`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: The fraction of accreted matter retained in nova eruption
- **Source**: `bse-interface/bse_interface.h`:873

### `--bse-gamma`

- **Type**: `double`
- **Default**: `-1.0`
- **Description**: Angular momentum factor for mass lost during Roche
- **Source**: `bse-interface/bse_interface.h`:875

### `--bse-hewind`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Helium star mass loss factor
- **Source**: `bse-interface/bse_interface.h`:866

### `--bse-ifflag`

- **Type**: `long long int`
- **Default**: `2`
- **Description**: if > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800
- **Source**: `bse-interface/bse_interface.h`:879

### `--bse-kmech`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: Kick mechanism; 1: standard momentum-conserving; 2: convection-asymmetry-driven; 3: collapse-asymmerty-driven; 4: neutrino driven
- **Source**: `bse-interface/bse_interface.h`:884

### `--bse-lambda`

- **Type**: `double`
- **Default**: `0.5`
- **Description**: Binding energy factor for common envelope evolution
- **Source**: `bse-interface/bse_interface.h`:869

### `--bse-metallicity`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: Metallicity Z, ranging from 0 to 0.03; when Z<0.0001, using EMP track, please make a symbolic link in the simulation directory to the track directory according to the --bse-trackmode option
- **Source**: `bse-interface/bse_interface.h`:897
- **Guards**: `BSEEMP` → --with-interrupt=bseEmp

### `--bse-mscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Mass scale factor from input data unit (IN) to Msun (m[Msun]=m[IN]*mscale)
- **Source**: `bse-interface/bse_interface.h`:894

### `--bse-mxns`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Helium star mass loss factor
- **Source**: `bse-interface/bse_interface.h`:867

### `--bse-neta`

- **Type**: `double`
- **Default**: `0.5`
- **Description**: Reimers mass-loss coefficent [neta*4x10^-13]
- **Source**: `bse-interface/bse_interface.h`:864

### `--bse-nsflag`

- **Type**: `long long int`
- **Default**: `3`
- **Description**: NS/BH foramtion options; 0: original SSE; 1: Belczynski (2002); 2: Belczynski (2008); 3: Fryer (2012) rapid SN; 4: Fryer (2012) delayed SN; 5: Eldridge & Tout (2004)
- **Source**: `bse-interface/bse_interface.h`:882

### `--bse-psflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: PPSN condition (Belczynski 2016); 0: no PPSN; 1: strong (Leung 2019); 2: moderate; 3: weak
- **Source**: `bse-interface/bse_interface.h`:883

### `--bse-pts1`

- **Type**: `double`
- **Default**: `0.05`
- **Description**: time step of MS
- **Source**: `bse-interface/bse_interface.h`:886

### `--bse-pts2`

- **Type**: `double`
- **Default**: `0.01`
- **Description**: time step of GB, CHeB, AGB, HeGB
- **Source**: `bse-interface/bse_interface.h`:887

### `--bse-pts3`

- **Type**: `double`
- **Default**: `0.02`
- **Description**: time step of HG, HeMS
- **Source**: `bse-interface/bse_interface.h`:888

### `--bse-rscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Radius scale factor from input data unit (IN) to Rsun (r[Rsun]=r[IN]*rscale)
- **Source**: `bse-interface/bse_interface.h`:893

### `--bse-sigma`

- **Type**: `double`
- **Default**: `265.0`
- **Description**: Dispersion in the Maxwellian for the SN kick speed [km/s]
- **Source**: `bse-interface/bse_interface.h`:876

### `--bse-tflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, activates tidal circularisation
- **Source**: `bse-interface/bse_interface.h`:878

### `--bse-trackmode`

- **Type**: `long long int`
- **Default**: `2`
- **Description**: star evolution option, need to make a soft link to the data directory in bse-interface/bseEmp/emptrack/: 1: L model (larger overshoot; directory name: ffbonn); 2: M model (smaller overshoot; directory name: ffgeneva); See details in Appendix A of Tanikawa et al. (2022, ApJ, 926, 83)
- **Source**: `bse-interface/bse_interface.h`:890
- **Guards**: `BSEEMP` → --with-interrupt=bseEmp

### `--bse-tscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Time scale factor from input data unit (IN) to Myr (time[Myr]=time[IN]*tscale)
- **Source**: `bse-interface/bse_interface.h`:892

### `--bse-vscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Velocity scale factor from input data unit(IN) to km/s (v[km/s]=v[IN]*vscale)
- **Source**: `bse-interface/bse_interface.h`:895

### `--bse-wdflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, uses WD IFMR of HPE, 1995, MNRAS, 272, 800
- **Source**: `bse-interface/bse_interface.h`:880

### `--bse-xi`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: wind accretion efficiency factor
- **Source**: `bse-interface/bse_interface.h`:871

### `--mobse-alpha`

- **Type**: `double`
- **Default**: `3.0`
- **Description**: Common-envelope efficiency parameter
- **Source**: `bse-interface/bse_interface.h`:909

### `--mobse-beta`

- **Type**: `double`
- **Default**: `0.125`
- **Description**: wind velocity factor: proportional to vwind**2
- **Source**: `bse-interface/bse_interface.h`:911

### `--mobse-bhflag`

- **Type**: `long long int`
- **Default**: `3`
- **Description**: BH kick option; 0: no kick; 1: same as NS; 2: scaled by fallback; 3: Giacobbo&Mapelli (2020)
- **Source**: `bse-interface/bse_interface.h`:923

### `--mobse-bwacc`

- **Type**: `double`
- **Default**: `1.5`
- **Description**: Bondi-Hoyle wind accretion factor
- **Source**: `bse-interface/bse_interface.h`:913

### `--mobse-cflag`

- **Type**: `long long int`
- **Default**: `0`
- **Description**: if =3, activates de Kool common-envelope model
- **Source**: `bse-interface/bse_interface.h`:919

### `--mobse-eddfac`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Eddington limit factor for mass transfer
- **Source**: `bse-interface/bse_interface.h`:915

### `--mobse-epsnov`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: The fraction of accreted matter retained in nova eruption
- **Source**: `bse-interface/bse_interface.h`:914

### `--mobse-gamma`

- **Type**: `double`
- **Default**: `-1.0`
- **Description**: Angular momentum factor for mass lost during Roche
- **Source**: `bse-interface/bse_interface.h`:916

### `--mobse-hewind`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Helium star mass loss factor
- **Source**: `bse-interface/bse_interface.h`:907

### `--mobse-lambda`

- **Type**: `double`
- **Default**: `0.1`
- **Description**: Binding energy factor for common envelope evolution
- **Source**: `bse-interface/bse_interface.h`:910

### `--mobse-metallicity`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: Metallicity
- **Source**: `bse-interface/bse_interface.h`:936

### `--mobse-msclae`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Mass scale factor from input data unit (IN) to Msun (m[Msun]=m[IN]*mscale)
- **Source**: `bse-interface/bse_interface.h`:934

### `--mobse-neta`

- **Type**: `double`
- **Default**: `0.5`
- **Description**: Reimers mass-loss coefficent [neta*4x10^-13]
- **Source**: `bse-interface/bse_interface.h`:905

### `--mobse-nsflag`

- **Type**: `long long int`
- **Default**: `3`
- **Description**: NS/BH formation options; 0: original SSE; 1: Belczynski (2008); 2: Fryer (2012) rapid SN; 3: Fryer (2012) delayed SN; 4: Belczynski (2008); 5: no SN explosion
- **Source**: `bse-interface/bse_interface.h`:924

### `--mobse-piflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: PPSN condition (Spera et al. 2015)
- **Source**: `bse-interface/bse_interface.h`:925

### `--mobse-pts1`

- **Type**: `double`
- **Default**: `0.05`
- **Description**: time step of MS
- **Source**: `bse-interface/bse_interface.h`:929

### `--mobse-pts2`

- **Type**: `double`
- **Default**: `0.01`
- **Description**: time step of GB, CHeB, AGB, HeGB
- **Source**: `bse-interface/bse_interface.h`:930

### `--mobse-pts3`

- **Type**: `double`
- **Default**: `0.02`
- **Description**: time step of HG, HeMS
- **Source**: `bse-interface/bse_interface.h`:931

### `--mobse-rscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Radius scale factor from input data unit (IN) to Rsun (r[Rsun]=r[IN]*rscale)
- **Source**: `bse-interface/bse_interface.h`:933

### `--mobse-sigma1`

- **Type**: `double`
- **Default**: `265.0`
- **Description**: Dispersion in the Maxwellian for the CCSN kick speed [km/s]
- **Source**: `bse-interface/bse_interface.h`:917

### `--mobse-sigma2`

- **Type**: `double`
- **Default**: `265.0`
- **Description**: Dispersion in the Maxwellian for the ECSN kick speed [km/s]
- **Source**: `bse-interface/bse_interface.h`:918

### `--mobse-tflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, activates tidal circularisation
- **Source**: `bse-interface/bse_interface.h`:920

### `--mobse-tscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Time scale factor from input data unit (IN) to Myr (time[Myr]=time[IN]*tscale)
- **Source**: `bse-interface/bse_interface.h`:932

### `--mobse-vsclae`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Velocity scale factor from input data unit(IN) to km/s (v[km/s]=v[IN]*vscale)
- **Source**: `bse-interface/bse_interface.h`:935

### `--mobse-wdflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, uses WD IFMR of HPE, 1995, MNRAS, 272, 800
- **Source**: `bse-interface/bse_interface.h`:922

### `--mobse-wind`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: Binary enhanced mass loss parameter, inactive for single
- **Source**: `bse-interface/bse_interface.h`:906

### `--mobse-xi`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: wind accretion efficiency factor
- **Source**: `bse-interface/bse_interface.h`:912

### `--p`

- **Type**: `std::string`
- **Default**: `"input.par"`
- **Description**: Input parameter file for sse/bse (this option should be used first before any other options)
- **Source**: `bse-interface/bse_interface.h`:901

## DSM Interrupt (`--with-interrupt=dsm`)

### `--G`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: gravitational constant
- **Source**: `src/disk_star_merger.hpp`:59

### `--dsm-dt-factor`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: time step factor for mass change calculation
- **Source**: `src/disk_star_merger.hpp`:61

### `--dsm-epsilon-bh`

- **Type**: `double`
- **Default**: `0.06`
- **Description**: the kenetic energy to radiation conversion efficiency of Eddington-limited accretion for BH; = 0: no accretion growth for BH
- **Source**: `src/disk_star_merger.hpp`:57

### `--dsm-epsilon-he`

- **Type**: `double`
- **Default**: `0.006`
- **Description**: helium enrichment scaling factor, used to calculate helium enrichment timescale; = 0: no helium enrichment
- **Source**: `src/disk_star_merger.hpp`:56

### `--dsm-epsilon-mdot`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: mass-change scaling factor, used to calculate mass-change timescale; = 0: no growth or mass loss
- **Source**: `src/disk_star_merger.hpp`:58

### `--dsm-he-disk`

- **Type**: `double`
- **Default**: `0.28`
- **Description**: helium fraction in the disk, used to calculate the equilbrium mass
- **Source**: `src/disk_star_merger.hpp`:54

### `--dsm-lambda0`

- **Type**: `double`
- **Default**: `0.75`
- **Description**: fraction of star's intrinsic luminosity over the Eddington luminosity without merger
- **Source**: `src/disk_star_merger.hpp`:53

### `--dsm-medd`

- **Type**: `double`
- **Default**: `253.3124306069483`
- **Description**: initial equilbrium mass of star
- **Source**: `src/disk_star_merger.hpp`:52

### `--dsm-merger-dm`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: mass loss rate for merger
- **Source**: `src/disk_star_merger.hpp`:47

### `--dsm-merger-tdelay`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: time delay for merger to increase mass
- **Source**: `src/disk_star_merger.hpp`:50

### `--dsm-new-star-mode`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: redistribute star mode, 0: no redistribute; 1: redistribute star position and velocity in random position along a circular orbit with the semi-major axis being the distance to the center; 2: redistribute star by choosing next type 3 star
- **Source**: `src/disk_star_merger.hpp`:62

### `--dsm-rstar-power`

- **Type**: `double`
- **Default**: `0.6`
- **Description**: stellar radius power index 'n', rs = s M^n
- **Source**: `src/disk_star_merger.hpp`:48

### `--dsm-rstar-scale`

- **Type**: `double`
- **Default**: `0.0046`
- **Description**: stellar radius scale 's', rs = s M^n
- **Source**: `src/disk_star_merger.hpp`:49

### `--dsm-salpeter-time`

- **Type**: `double`
- **Default**: `0`
- **Description**: salpeter timescale, for the star to reach equilbrium, if zero, no stellar evolution
- **Source**: `src/disk_star_merger.hpp`:55

### `--dsm-seed-mass`

- **Type**: `double`
- **Default**: `10.0`
- **Description**: initial mass of star seed
- **Source**: `src/disk_star_merger.hpp`:51

### `--dsm-speed-of-light`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: speed of light
- **Source**: `src/disk_star_merger.hpp`:60

### `--p`

- **Type**: `std::string`
- **Default**: `"input.par"`
- **Description**: Input parameter file for external force (this option should be used first before any other options)
- **Source**: `src/disk_star_merger.hpp`:63

## Galpy External (`--with-external=galpy`)

### `--galpy-conf-file`

- **Type**: `std::string`
- **Default**: `"__NONE__"`
- **Description**: A configure file of the time- and space-dependent potentials; Use petar.galpy.help to check how to generate configure file.
- **Source**: `galpy-interface/galpy_interface.h`:36

### `--galpy-fscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Acceleration scale factor (vscale^2/rscale) from unit of the input particle data (IN) to Galpy acceleration unit (1.0)
- **Source**: `galpy-interface/galpy_interface.h`:41

### `--galpy-pscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Potential scale factor (vscale^2) from unit of the input particle data (IN) to Galpy potential unit (1.0)
- **Source**: `galpy-interface/galpy_interface.h`:42

### `--galpy-rscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Radius scale factor from unit of the input particle data (IN) to Galpy distance unit (1.0)
- **Source**: `galpy-interface/galpy_interface.h`:38

### `--galpy-set`

- **Type**: `std::string`
- **Default**: `"__NONE__"`
- **Description**: Add a pre-defined potential set to the potential list, options are: MWPotential2014
- **Source**: `galpy-interface/galpy_interface.h`:35

### `--galpy-tscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Time scale factor (rscale/vscale) from unit of the input particle data (IN) to Galpy time (1.0)
- **Source**: `galpy-interface/galpy_interface.h`:39

### `--galpy-type-arg`

- **Type**: `std::string`
- **Default**: `"__NONE__"`
- **Description**: Add potential types and arguments to the potential list in the center of the galactic reference frame; Use petar.galpy.help to check how to setup types and arguments.
- **Source**: `galpy-interface/galpy_interface.h`:34

### `--galpy-units`

- **Type**: `std::string`
- **Default**: `"unscale"`
- **Description**: Units conversion set: 'unscale': no conversion; all scale factors are 1.0; 'bovy': radial unit: 8 kpc; velocity unit: 220 km/s
- **Source**: `galpy-interface/galpy_interface.h`:37

### `--galpy-vscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Velocity scale factor from unit of the input particle data (IN) to Galpy velocity unit (1.0)
- **Source**: `galpy-interface/galpy_interface.h`:40

### `--p`

- **Type**: `std::string`
- **Default**: `"input.par"`
- **Description**: Input parameter file for Galpy (this option should be used first before any other options)
- **Source**: `galpy-interface/galpy_interface.h`:43

## Agama External (`--with-external=agama`)

### `--agama-conf-file`

- **Type**: `std::string`
- **Default**: `"__NONE__"`
- **Description**: A configure file of agama potential
- **Source**: `agama-interface/agama_interface.h`:23

### `--agama-rscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Radius scale factor from unit of the input particle data (IN) to Agama distance unit
- **Source**: `agama-interface/agama_interface.h`:24

### `--agama-vscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Velocity scale factor from unit of the input particle data (IN) to Agama velocity unit
- **Source**: `agama-interface/agama_interface.h`:25

## Gas Drag (`--with-external-hard=gasdrag`)

### `--G`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Gravitational constant
- **Source**: `src/gas_drag.hpp`:59

### `--gdf-K`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Polytropic constant in units of PeTar input, used to evaluate Pressure P = K rho^gamma
- **Source**: `src/gas_drag.hpp`:47

### `--gdf-coulomb-log`

- **Type**: `double`
- **Default**: `3.1`
- **Description**: coulomb logarithm
- **Source**: `src/gas_drag.hpp`:46

### `--gdf-gamma`

- **Type**: `double`
- **Default**: `4.0/3.0`
- **Description**: Polytropic exponent, used to evaluate Pressure P = K rho^gamma
- **Source**: `src/gas_drag.hpp`:48

### `--gdf-hard-mode`

- **Type**: `long long int`
- **Default**: `0`
- **Description**: GDF external force for hard integration; 0, not used; 1, gas dynamical friction (Ostriker 1999, Rozner et al. 2022); 2, gas dynamical friction only in radial direction
- **Source**: `src/gas_drag.hpp`:44

### `--gdf-ifunc-mach-lower`

- **Type**: `double`
- **Default**: `0.9`
- **Description**: lower Mach boundary for Ifunc Hermite interpolation; recommended range: about 0.85 (or smaller) for derivative order = 3, typical 0.9 for order = 2
- **Source**: `src/gas_drag.hpp`:49

### `--gdf-ifunc-mach-upper`

- **Type**: `double`
- **Default**: `1.1`
- **Description**: upper Mach boundary for Ifunc Hermite interpolation; recommended range: about 1.15 (or larger) for derivative order = 3, typical 1.1 for order = 2
- **Source**: `src/gas_drag.hpp`:50

### `--gdf-ifunc-smooth-order`

- **Type**: `long long int`
- **Default**: `2`
- **Description**: Ifunc Hermite smooth derivative order at boundaries (2 or 3); recommended: 2 for robustness, 3 with a wider Mach range (e.g. 0.85-1.15 or wider)
- **Source**: `src/gas_drag.hpp`:51

### `--gdf-sound-speed`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: sound speed in units of PeTar input, if given 0, calculate by sqrt(P/rho), based on hydrostatic equilibrium
- **Source**: `src/gas_drag.hpp`:45

### `--p`

- **Type**: `std::string`
- **Default**: `"input.par"`
- **Description**: Input parameter file for external force (this option should be used first before any other options)
- **Source**: `src/gas_drag.hpp`:60

## External Hard (`--with-external-hard`)

### `--ext-hard-switch`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: switch of external hard force; 0: off, 1: on
- **Source**: `src/external_hard.hpp`:30

### `--p`

- **Type**: `std::string`
- **Default**: `"input.par"`
- **Description**: Input parameter file for external force (this option should be used first before any other options)
- **Source**: `src/external_hard.hpp`:31


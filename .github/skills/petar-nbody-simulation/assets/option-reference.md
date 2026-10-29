# PeTar CLI Option Reference (Source-Generated)

All possible PeTar command-line options across every configure variant.
Generated from source headers — no binary compilation needed.

> Generated: 2026-10-29 05:48 UTC
> Regenerate: `python3 .github/skills/petar-nbody-simulation/assets/generate_option_reference.py`

Use `petar.select` to list installed binary families and `<binary> -h` for exact runtime option validation.

---

## Summary

| Configure Flag Group | Option Count |
|-----------------------|-------------|
| Core (always available) | 78 |
| BSE Interrupt (`--with-interrupt=bse|mobse|bseEmp`) | 63 |
| DSM Interrupt (`--with-interrupt=dsm`) | 19 |
| Galpy External (`--with-external=galpy`) | 10 |
| Agama External (`--with-external=agama`) | 3 |
| Gas Drag (`--with-external-hard=gasdrag`) | 10 |
| External Hard (`--with-external-hard`) | 2 |
| **Total** | **185** |

---

## Core (always available)

### `--G`

- **Type**: `PS::F64`
- **Default**: `1.0`
- **Description**: Gravitational constant, if -u 1, G =
- **Source**: `src/petar.hpp`:164

### `--T`

- **Type**: `PS::F64`
- **Default**: `0.3`
- **Description**: Particle-tree opening angle theta
- **Source**: `src/petar.hpp`:152

### `--a`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Data file output mode; 0: overwrite files except object dump files, include header lines; 1: append files except snapshots, no header line
- **Source**: `src/petar.hpp`:182

### `--ar-ds-scale`

- **Type**: `PS::F64`
- **Default**: `1.0`
- **Description**: Scale factor for SDAR step size calculation
- **Source**: `src/hard.hpp`:153

### `--ar-g-func`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: SDAR g-function (time transformation) method; 0: standard LogH; 1: BTLogH
- **Source**: `src/hard.hpp`:156
- **Guards**: `AR_G_FUNC`

### `--ar-max-error`

- **Type**: `PS::F64`
- **Default**: `1e-8`
- **Description**: Maximum energy error allowed for the SDAR integrator
- **Source**: `src/hard.hpp`:150

### `--ar-max-nstep`

- **Type**: `PS::S64`
- **Default**: `1000000`
- **Description**: Maximum step allowed for the SDAR sym integrator
- **Source**: `src/hard.hpp`:151

### `--ar-slowdown-factor`

- **Type**: `PS::F64`
- **Default**: `1e-4`
- **Description**: Slowdown perturbation criterion
- **Source**: `src/hard.hpp`:154

### `--ar-sym-order`

- **Type**: `PS::S64`
- **Default**: `-6`
- **Description**: Order of the symplectic integrator for SDAR, should be even number; -6,-8: Yoshida 2nd symplectic method; 4,6,8,...: Yoshida 1st symplectic method
- **Source**: `src/hard.hpp`:152

### `--b`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Number of primordial binaries (n_bin) for initialization (assuming the binaries' IDs are 1,2*n_bin)
- **Source**: `src/petar.hpp`:161

### `--center-id`

- **Type**: `PS::S64`
- **Default**: `-1`
- **Description**: id of the central object for a system like a stellar disk
- **Source**: `src/hard.hpp`:145

### `--detect-interrupt`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Stellar evolution of binaries in SDAR integration; 0: switch off; 1: using BSE based code (if '--stellar-evolution != 0)
- **Source**: `src/hard.hpp`:163
- **Guards**: `BSE_BASE` → --with-interrupt (bse/mobse/bseEmp), `STELLAR_EVOLUTION` → --with-interrupt

### `--domain-nstep`

- **Type**: `PS::S64`
- **Default**: `16`
- **Description**: Number of steps between domain decompositions
- **Source**: `src/petar.hpp`:186
- **Guards**: `PARTICLE_SIMULATOR_MPI_PARALLEL` → MPI build

### `--domain-weight-mode`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Domain decomposition weight mode for MPI parallel; 0: equal weight for each MPI processor; 1: use force calculation time as weight to obtain better load balance with losing simulation reproducibility
- **Source**: `src/petar.hpp`:185
- **Guards**: `PARTICLE_SIMULATOR_MPI_PARALLEL` → MPI build

### `--dt-soft-kepler-nstep`

- **Type**: `PS::F64`
- **Default**: `16.0`
- **Description**: Factor 'nstep' to determine dt_soft by P(r_in)/nstep, see option '-s' and '-r'
- **Source**: `src/petar.hpp`:170
- **Guards**: `KDKDK_4TH` → 4th-order KDKDK step mode

### `--dt-soft-sigma-factor`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Factor 'alpha' to determine dt_soft by alpha*r_in/sigma_3D, see option '-s' and '-r'; = 0: not used, apply --dt-soft-kepler-nstep; > 0: use this option instead of '--dt-soft-kepler-nstep'
- **Source**: `src/petar.hpp`:174

### `--energy-err-hard`

- **Type**: `PS::F64`
- **Default**: `1e-4`
- **Description**: Maximum energy error allowed for the hard integrator
- **Source**: `src/hard.hpp`:129
- **Guards**: `HARD_CHECK_ENERGY` → development debug

### `--f`

- **Type**: `std::string`
- **Default**: `"data"`
- **Description**: Prefix of filenames for output data: [prefix].**
- **Source**: `src/petar.hpp`:188

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
- **Description**: Square acceleration offset for Hermite time step calculation to avoid too small step when weak acceleration exists; = -1: calculate from mean mass <m> and r_out (G*<m>/r_out^2)^2; = 0: no offset; > 0: custom offset value; Note: on restart with -s or -r given and this option not, the stored value is re-determined
- **Source**: `src/hard.hpp`:138

### `--hermite-de-crit`

- **Type**: `PS::F64`
- **Default**: `1e-4`
- **Description**: Ekin change rate criterion for reinitializing hermite time step
- **Source**: `src/hard.hpp`:159

### `--hermite-dm-crit`

- **Type**: `PS::F64`
- **Default**: `1e-4`
- **Description**: Mass change rate criterion for reinitializing hermite time step
- **Source**: `src/hard.hpp`:158

### `--hermite-dt-max`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Maximum hermite timestep; = 0: auto-determine by one tree drift step;      KDKDK tree step (4th/2nd order): 0.5*dt_soft;      KDK tree step (2nd order): dt_soft; > 0: custom value, must not exceed one tree drift step; Note: on restart with -s or -r given and this option not, the stored value is re-determined
- **Source**: `src/hard.hpp`:148

### `--hermite-dt-min-index`

- **Type**: `PS::S64`
- **Default**: `40`
- **Description**: Power index n for the smallest timestep (0.5^n) allowed in the Hermite integrator
- **Source**: `src/hard.hpp`:149

### `--hermite-eta`

- **Type**: `PS::F64`
- **Default**: `0.1`
- **Description**: Hermite timestep coefficient eta
- **Source**: `src/hard.hpp`:146

### `--hermite-eta-init`

- **Type**: `PS::F64`
- **Default**: `0.001`
- **Description**: Hermite timestep coefficient eta for initial step in 2nd order
- **Source**: `src/hard.hpp`:147

### `--hermite-n-neighbor-max`

- **Type**: `PS::S64`
- **Default**: `300`
- **Description**: Initial capacity of the group neighbor address lists (auto-grown on demand, not a limit)
- **Source**: `src/hard.hpp`:160

### `--i`

- **Type**: `PS::S64`
- **Default**: `2`
- **Description**: Data file reading and writing format; snapshots, status and escaper outputs follow the write mode selected here; 0: read and write in BINARY; 1: read and write in ASCII; 2: read in ASCII, write in BINARY; 3: read in BINARY, write in ASCII
- **Source**: `src/petar.hpp`:180

### `--id-offset`

- **Type**: `PS::S64`
- **Default**: `-1`
- **Description**: Starting ID for artificial particles, total number of real particles must always be smaller than this
- **Source**: `src/hard.hpp`:144

### `--kdtree-n-particles-min`

- **Type**: `PS::S64`
- **Default**: `32`
- **Description**: Minimum number of particles + groups for building kdtree to speed up neighbor search in Hermite-only neighbor force calculation
- **Source**: `src/hard.hpp`:186
- **Guards**: `HERMITE_ONLY_CALC_NEIGHBOR_FORCE` → --enable-hermite-only-calc-neighbor-force

### `--keep-tmp-on-startup`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Startup tmp handling for transactional outputs; 0: auto-remove residual tmp files before run (default); 1: keep residual tmp files and only print warning
- **Source**: `src/petar.hpp`:183

### `--n`

- **Type**: `PS::S64`
- **Default**: `100000`
- **Description**: Total number of particles, used only when the input data filename is __Plummer
- **Source**: `src/petar.hpp`:165

### `--n-sample-average`

- **Type**: `PS::S64`
- **Default**: `100`
- **Description**: Average target number of sample particles per process
- **Source**: `src/petar.hpp`:160

### `--o`

- **Type**: `PS::F64`
- **Default**: `1.0`
- **Description**: Output time interval for particle dataset snapshots
- **Source**: `src/petar.hpp`:179

### `--os-nsplit`

- **Type**: `PS::S64`
- **Default**: `4`
- **Description**: Number of binary sample points for tree perturbation force using orbit-sampling method
- **Source**: `src/hard.hpp`:142
- **Guards**: `ORBIT_SAMPLING` → --enable-orbit-sampling

### `--p`

- **Type**: `std::string`
- **Default**: `"data.par"`
- **Description**: Input parameter file (this option should be used first before any other options); The auto-determined filename is '[prefix of output filename].par', in default, the prefix is 'data'.
- **Source**: `src/petar.hpp`:189

### `--pn-c`

- **Type**: `PS::F64`
- **Default**: `1`
- **Description**: speed of light value for Post Newtonian; if -u 1 is used, auto determined
- **Source**: `src/hard.hpp`:177

### `--pn-crit-ar`

- **Type**: `PS::F64`
- **Default**: `1e-6`
- **Description**: AR speed criterion to switch on PN terms, min (v/c)^2
- **Source**: `src/hard.hpp`:183
- **Guards**: `SDAR_PN` → --with-pn (SDAR)

### `--pn-crit-hermite`

- **Type**: `PS::F64`
- **Default**: `1e-6`
- **Description**: Hermite speed criterion to switch on PN terms, min (v/c)^2
- **Source**: `src/hard.hpp`:180
- **Guards**: `HERMITE_PN` → --with-pn (Hermite)

### `--r`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Outer changeover radius (r_out); > 0: custom r_out value and check '-s dt_soft';      dt_soft = 0: calculate dt_soft and then adjust r_out by dt_soft;      dt_soft > 0: use custom r_out directly; = 0 (default): check '--dt-soft-sigma-factor alpha':;      alpha > 0: r_out = alpha*dt_soft*sigma_3D/r-ratio;          sigma_3D: global 3D velocity dispersion;          r-ratio: defined by --r-ratio;      alpha = 0 (default): r_out = a(r_in)/r-ratio;          a(r_in): the binary semi-major axis with the period of nstep*dt_soft;          nstep: defined by --dt-soft-kepler-nstep; Note: on restart with this option given, r-search-min, r-search-group, r-group,;      hermite-acc-offset-sq and dt_soft (if -s not given) are re-determined
- **Source**: `src/petar.hpp`:167

### `--r-escape`

- **Type**: `PS::F64`
- **Default**: `PS::LARGE_FLOAT`
- **Description**: Object escape radius criterion; < 0: remove objects when r>-r_escape; >= 0: remove objects when r>r_escape and energy>0
- **Source**: `src/petar.hpp`:178

### `--r-group`

- **Type**: `PS::F64`
- **Default**: `-1.0`
- **Description**: Tidal tensor box size and the radial criterion for detecting multiple groups (binaries, triples, etc.); = -1: auto-determine by 0.8*r_search_group; = 0: switch off SDAR; > 0: custom criterion value; Note: on restart with -s, -r, --r-ratio or --r-search-group-safety given and this option not, the stored value is re-determined
- **Source**: `src/hard.hpp`:133

### `--r-ratio`

- **Type**: `PS::F64`
- **Default**: `0.1`
- **Description**: Ratio between inner (r_in) and outer (r_out) changeover radii; Note: on restart with this option given, r-search-group and r-group are re-determined
- **Source**: `src/petar.hpp`:168

### `--r-search-group`

- **Type**: `PS::F64`
- **Default**: `-1.0`
- **Description**: The radial criterion for detecting multiple group candidates; = -1: auto-determine by 1.0*r_in (or by sigma_1D based hard-soft boundary if sigma_1D is provided); = 0: switch off SDAR; > 0: custom criterion value; Note: on restart with -s, -r, --r-ratio or --r-search-group-safety given and this option not, the stored value is re-determined
- **Source**: `src/hard.hpp`:134

### `--r-search-group-safety`

- **Type**: `PS::F64`
- **Default**: `1.0`
- **Description**: Safety factor for physics-based SDAR group search: r_search_group = min(r_in, safety * G * m_avg / sigma_1D^2); only used when r-search-group=-1 and sigma_1D is provided; Note: on restart with this option given, r-search-group and r-group are re-determined
- **Source**: `src/hard.hpp`:135

### `--r-search-min`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Minimum neighbor search radius for hard clusters; = 0: auto-determine by max(search-vel-factor*sigma_1D*dt_soft + rout, 1.2 r_out); > 0: custom search radius value; Note: on restart with -s, -r or --r-search-vel-factor given, the stored value is re-determined
- **Source**: `src/petar.hpp`:177

### `--r-search-peri-factor`

- **Type**: `PS::F64`
- **Default**: `1.5`
- **Description**: Neighbor search coefficient for periapsis check
- **Source**: `src/petar.hpp`:176

### `--r-search-vel-factor`

- **Type**: `PS::F64`
- **Default**: `3.0`
- **Description**: Neighbor search coefficient for velocity check (v*dt)
- **Source**: `src/petar.hpp`:175

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
- **Source**: `src/hard.hpp`:173

### `--record-id-end-two`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Ending of the 2nd id range for hard dump; notice that the ending id is not included in hard dump
- **Source**: `src/hard.hpp`:175

### `--record-id-start-one`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Starting of the first id range for hard dump recording every tree step, save into files <prefix>.object_[id]
- **Source**: `src/hard.hpp`:172

### `--record-id-start-two`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Starting of the 2nd id range for hard dump recording every tree step
- **Source**: `src/hard.hpp`:174

### `--s`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Tree timestep (dt_soft); > 0: custom dt_soft value, regularized to 0.5^n, where n is an integer; = 0: check '-r r_out':;      r_out = 0 (default): dt_soft = 2.6E-4*GM/sigma_3D^3, and is regularized to 0.5^n;          sigma_3D: global 3D velocity dispersion;      r_out > 0: check '--dt-soft-sigma-factor alpha':;          alpha > 0: dt_soft = alpha*r_in/(sqrt(3)*sigma);              r_in: determined by --r-ratio and r_out;          alpha = 0 (default): dt_soft = P(r_in)/nstep;              P(r_in): the binary period with the semi-major axis of r_in;              nstep: defined by --dt-soft-kepler-nstep; Note: on restart with this option given, hermite-dt-max, r-search-min and r_out (if -r not given) are re-determined
- **Source**: `src/petar.hpp`:166

### `--snap-filename`

- **Type**: `std::string`
- **Default**: `"__NONE__"`
- **Description**: Input data file
- **Source**: `src/petar.hpp`:190

### `--soft-eps`

- **Type**: `PS::F64`
- **Default**: `0.0`
- **Description**: Softening epsilon
- **Source**: `src/hard.hpp`:132

### `--stellar-evolution`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Stellar evolution of stars in Hermite and SDAR integration; 0: switch off; >=1: using SSE/BSE based codes; 2: activate dynamical tide and hyperbolic GW radiation
- **Source**: `src/hard.hpp`:164
- **Guards**: `BSE_BASE` → --with-interrupt (bse/mobse/bseEmp), `STELLAR_EVOLUTION` → --with-interrupt

### `--t`

- **Type**: `PS::F64`
- **Default**: `10.0`
- **Description**: End time of simulation
- **Source**: `src/petar.hpp`:162

### `--tree-ngroup-limit`

- **Type**: `PS::S64`
- **Default**: `1024`
- **Description**: Particle-tree group number limit; Optimal value for x86-AVX512 is 1024
- **Source**: `src/petar.hpp`:155
- **Guards**: `USE__AVX512`

### `--tree-nleaf-limit`

- **Type**: `PS::S64`
- **Default**: `20`
- **Description**: Particle-tree leaf number limit; Optimal value should be slightly >= artificial particle number (tidal tensor 8 + anti-force sample 3) + 2 (binary member) + 1 (binary c.m.)
- **Source**: `src/petar.hpp`:153

### `--tree-nstep-mklist`

- **Type**: `PS::S64`
- **Default**: `2`
- **Description**: Particle-tree make-list interval in number of soft-step
- **Source**: `src/petar.hpp`:159

### `--tt-nstep`

- **Type**: `PS::S64`
- **Default**: `4`
- **Description**: Number of steps per slow-down binary orbits (period/dt_soft) for isolated binaries; also the maximum criterion for activating tidal tensor method
- **Source**: `src/hard.hpp`:139

### `--tt-switch`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Tidal tensor calculation for (counter-)perturbation (from)on binaries: 0: off, 1: on
- **Source**: `src/hard.hpp`:140

### `--u`

- **Type**: `PS::S64`
- **Default**: `0`
- **Description**: Input data unit; 0: based on the value of G; 1: mass:Msun, length:pc, time:Myr, velocity:pc/Myr, modify G to fit this unit set
- **Source**: `src/petar.hpp`:163

### `--w`

- **Type**: `PS::S64`
- **Default**: `1`
- **Description**: Data file writing style; 0: no output; 1: write all files separately; 2. write snapshots in status files in one line per step (no MPI support); 3. write files except snapshots
- **Source**: `src/petar.hpp`:181

### `--write-group-info`

- **Type**: `PS::S64`
- **Default**: `2`
- **Description**: Write information of new and end groups; 0: no output; 1: ascii output; 2: binary output, files are [data filename prefix].group.[MPI rank].n[N_member]
- **Source**: `src/hard.hpp`:170
- **Guards**: `ADJUST_GROUP_PRINT` → --enable-adjust-group-print

## BSE Interrupt (`--with-interrupt=bse|mobse|bseEmp`)

### `--bse-alpha`

- **Type**: `double`
- **Default**: `3.0`
- **Description**: Common-envelope efficiency parameter
- **Source**: `bse-interface/bse_interface.h`:1057
- **Guards**: `!SEVN` (negated)

### `--bse-beta`

- **Type**: `double`
- **Default**: `0.125`
- **Description**: wind velocity factor: proportional to vwind**2
- **Source**: `bse-interface/bse_interface.h`:1059
- **Guards**: `!SEVN` (negated)

### `--bse-bhflag`

- **Type**: `long long int`
- **Default**: `2`
- **Description**: BH kick option: 0: no kick; 1: same as NS; 2: scaled by fallback
- **Source**: `bse-interface/bse_interface.h`:1070
- **Guards**: `!SEVN` (negated)

### `--bse-bhwacc`

- **Type**: `double`
- **Default**: `1.5`
- **Description**: Bondi-Hoyle wind accretion factor
- **Source**: `bse-interface/bse_interface.h`:1061
- **Guards**: `!SEVN` (negated)

### `--bse-bwind`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: Binary enhanced mass loss parameter, inactive for single
- **Source**: `bse-interface/bse_interface.h`:1054
- **Guards**: `!SEVN` (negated)

### `--bse-ceflag`

- **Type**: `long long int`
- **Default**: `0`
- **Description**: if =3, activates de Kool common-envelope model
- **Source**: `bse-interface/bse_interface.h`:1066
- **Guards**: `!SEVN` (negated)

### `--bse-ecflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, ECS is switched on
- **Source**: `bse-interface/bse_interface.h`:1074

### `--bse-eddfac`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Eddington limit factor for mass transfer
- **Source**: `bse-interface/bse_interface.h`:1063
- **Guards**: `!SEVN` (negated)

### `--bse-epsnov`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: The fraction of accreted matter retained in nova eruption
- **Source**: `bse-interface/bse_interface.h`:1062
- **Guards**: `!SEVN` (negated)

### `--bse-gamma`

- **Type**: `double`
- **Default**: `-1.0`
- **Description**: Angular momentum factor for mass lost during Roche
- **Source**: `bse-interface/bse_interface.h`:1064
- **Guards**: `!SEVN` (negated)

### `--bse-hewind`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Helium star mass loss factor
- **Source**: `bse-interface/bse_interface.h`:1055
- **Guards**: `!SEVN` (negated)

### `--bse-ifflag`

- **Type**: `long long int`
- **Default**: `2`
- **Description**: if > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800
- **Source**: `bse-interface/bse_interface.h`:1068
- **Guards**: `!SEVN` (negated)

### `--bse-kmech`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: Kick mechanism; 1: standard momentum-conserving; 2: convection-asymmetry-driven; 3: collapse-asymmerty-driven; 4: neutrino driven
- **Source**: `bse-interface/bse_interface.h`:1073

### `--bse-lambda`

- **Type**: `double`
- **Default**: `0.5`
- **Description**: Binding energy factor for common envelope evolution
- **Source**: `bse-interface/bse_interface.h`:1058
- **Guards**: `!SEVN` (negated)

### `--bse-metallicity`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: Metallicity Z, ranging from 0 to 0.03; when Z<0.0001, using EMP track, please make a symbolic link in the simulation directory to the track directory according to the --bse-trackmode option
- **Source**: `bse-interface/bse_interface.h`:1086
- **Guards**: `BSEEMP` → --with-interrupt=bseEmp

### `--bse-mscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Mass scale factor from input data unit (IN) to Msun (m[Msun]=m[IN]*mscale)
- **Source**: `bse-interface/bse_interface.h`:1083

### `--bse-mxns`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Helium star mass loss factor
- **Source**: `bse-interface/bse_interface.h`:1056
- **Guards**: `!SEVN` (negated)

### `--bse-neta`

- **Type**: `double`
- **Default**: `0.5`
- **Description**: Reimers mass-loss coefficent [neta*4x10^-13]
- **Source**: `bse-interface/bse_interface.h`:1053
- **Guards**: `!SEVN` (negated)

### `--bse-nsflag`

- **Type**: `long long int`
- **Default**: `3`
- **Description**: NS/BH foramtion options; 0: original SSE; 1: Belczynski (2002); 2: Belczynski (2008); 3: Fryer (2012) rapid SN; 4: Fryer (2012) delayed SN; 5: Eldridge & Tout (2004)
- **Source**: `bse-interface/bse_interface.h`:1071
- **Guards**: `!SEVN` (negated)

### `--bse-psflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: PPSN condition (Belczynski 2016); 0: no PPSN; 1: strong (Leung 2019); 2: moderate; 3: weak
- **Source**: `bse-interface/bse_interface.h`:1072

### `--bse-pts1`

- **Type**: `double`
- **Default**: `0.05`
- **Description**: time step of MS
- **Source**: `bse-interface/bse_interface.h`:1075
- **Guards**: `!SEVN` (negated)

### `--bse-pts2`

- **Type**: `double`
- **Default**: `0.01`
- **Description**: time step of GB, CHeB, AGB, HeGB
- **Source**: `bse-interface/bse_interface.h`:1076
- **Guards**: `!SEVN` (negated)

### `--bse-pts3`

- **Type**: `double`
- **Default**: `0.02`
- **Description**: time step of HG, HeMS
- **Source**: `bse-interface/bse_interface.h`:1077
- **Guards**: `!SEVN` (negated)

### `--bse-rscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Radius scale factor from input data unit (IN) to Rsun (r[Rsun]=r[IN]*rscale)
- **Source**: `bse-interface/bse_interface.h`:1082

### `--bse-sigma`

- **Type**: `double`
- **Default**: `265.0`
- **Description**: Dispersion in the Maxwellian for the SN kick speed [km/s]
- **Source**: `bse-interface/bse_interface.h`:1065

### `--bse-tflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, activates tidal circularisation
- **Source**: `bse-interface/bse_interface.h`:1067
- **Guards**: `!SEVN` (negated)

### `--bse-trackmode`

- **Type**: `long long int`
- **Default**: `2`
- **Description**: star evolution option, need to make a soft link to the data directory in bse-interface/bseEmp/emptrack/: 1: L model (larger overshoot; directory name: ffbonn); 2: M model (smaller overshoot; directory name: ffgeneva); See details in Appendix A of Tanikawa et al. (2022, ApJ, 926, 83)
- **Source**: `bse-interface/bse_interface.h`:1079
- **Guards**: `BSEEMP` → --with-interrupt=bseEmp

### `--bse-tscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Time scale factor from input data unit (IN) to Myr (time[Myr]=time[IN]*tscale)
- **Source**: `bse-interface/bse_interface.h`:1081

### `--bse-vscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Velocity scale factor from input data unit(IN) to km/s (v[km/s]=v[IN]*vscale)
- **Source**: `bse-interface/bse_interface.h`:1084

### `--bse-wdflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, uses WD IFMR of HPE, 1995, MNRAS, 272, 800
- **Source**: `bse-interface/bse_interface.h`:1069
- **Guards**: `!SEVN` (negated)

### `--bse-xi`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: wind accretion efficiency factor
- **Source**: `bse-interface/bse_interface.h`:1060
- **Guards**: `!SEVN` (negated)

### `--mobse-alpha`

- **Type**: `double`
- **Default**: `3.0`
- **Description**: Common-envelope efficiency parameter
- **Source**: `bse-interface/bse_interface.h`:1098
- **Guards**: `!SEVN` (negated)

### `--mobse-beta`

- **Type**: `double`
- **Default**: `0.125`
- **Description**: wind velocity factor: proportional to vwind**2
- **Source**: `bse-interface/bse_interface.h`:1100
- **Guards**: `!SEVN` (negated)

### `--mobse-bhflag`

- **Type**: `long long int`
- **Default**: `3`
- **Description**: BH kick option; 0: no kick; 1: same as NS; 2: scaled by fallback; 3: Giacobbo&Mapelli (2020)
- **Source**: `bse-interface/bse_interface.h`:1112
- **Guards**: `!SEVN` (negated)

### `--mobse-bwacc`

- **Type**: `double`
- **Default**: `1.5`
- **Description**: Bondi-Hoyle wind accretion factor
- **Source**: `bse-interface/bse_interface.h`:1102
- **Guards**: `!SEVN` (negated)

### `--mobse-cflag`

- **Type**: `long long int`
- **Default**: `0`
- **Description**: if =3, activates de Kool common-envelope model
- **Source**: `bse-interface/bse_interface.h`:1108
- **Guards**: `!SEVN` (negated)

### `--mobse-eddfac`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Eddington limit factor for mass transfer
- **Source**: `bse-interface/bse_interface.h`:1104
- **Guards**: `!SEVN` (negated)

### `--mobse-epsnov`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: The fraction of accreted matter retained in nova eruption
- **Source**: `bse-interface/bse_interface.h`:1103
- **Guards**: `!SEVN` (negated)

### `--mobse-gamma`

- **Type**: `double`
- **Default**: `-1.0`
- **Description**: Angular momentum factor for mass lost during Roche
- **Source**: `bse-interface/bse_interface.h`:1105
- **Guards**: `!SEVN` (negated)

### `--mobse-hewind`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Helium star mass loss factor
- **Source**: `bse-interface/bse_interface.h`:1096
- **Guards**: `!SEVN` (negated)

### `--mobse-lambda`

- **Type**: `double`
- **Default**: `0.1`
- **Description**: Binding energy factor for common envelope evolution
- **Source**: `bse-interface/bse_interface.h`:1099
- **Guards**: `!SEVN` (negated)

### `--mobse-metallicity`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: Metallicity
- **Source**: `bse-interface/bse_interface.h`:1125

### `--mobse-msclae`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Mass scale factor from input data unit (IN) to Msun (m[Msun]=m[IN]*mscale)
- **Source**: `bse-interface/bse_interface.h`:1123

### `--mobse-neta`

- **Type**: `double`
- **Default**: `0.5`
- **Description**: Reimers mass-loss coefficent [neta*4x10^-13]
- **Source**: `bse-interface/bse_interface.h`:1094
- **Guards**: `!SEVN` (negated)

### `--mobse-nsflag`

- **Type**: `long long int`
- **Default**: `3`
- **Description**: NS/BH formation options; 0: original SSE; 1: Belczynski (2008); 2: Fryer (2012) rapid SN; 3: Fryer (2012) delayed SN; 4: Belczynski (2008); 5: no SN explosion
- **Source**: `bse-interface/bse_interface.h`:1113
- **Guards**: `!SEVN` (negated)

### `--mobse-piflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: PPSN condition (Spera et al. 2015)
- **Source**: `bse-interface/bse_interface.h`:1114

### `--mobse-pts1`

- **Type**: `double`
- **Default**: `0.05`
- **Description**: time step of MS
- **Source**: `bse-interface/bse_interface.h`:1118
- **Guards**: `!SEVN` (negated)

### `--mobse-pts2`

- **Type**: `double`
- **Default**: `0.01`
- **Description**: time step of GB, CHeB, AGB, HeGB
- **Source**: `bse-interface/bse_interface.h`:1119
- **Guards**: `!SEVN` (negated)

### `--mobse-pts3`

- **Type**: `double`
- **Default**: `0.02`
- **Description**: time step of HG, HeMS
- **Source**: `bse-interface/bse_interface.h`:1120
- **Guards**: `!SEVN` (negated)

### `--mobse-rscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Radius scale factor from input data unit (IN) to Rsun (r[Rsun]=r[IN]*rscale)
- **Source**: `bse-interface/bse_interface.h`:1122

### `--mobse-sigma1`

- **Type**: `double`
- **Default**: `265.0`
- **Description**: Dispersion in the Maxwellian for the CCSN kick speed [km/s]
- **Source**: `bse-interface/bse_interface.h`:1106

### `--mobse-sigma2`

- **Type**: `double`
- **Default**: `265.0`
- **Description**: Dispersion in the Maxwellian for the ECSN kick speed [km/s]
- **Source**: `bse-interface/bse_interface.h`:1107

### `--mobse-tflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, activates tidal circularisation
- **Source**: `bse-interface/bse_interface.h`:1109
- **Guards**: `!SEVN` (negated)

### `--mobse-tscale`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Time scale factor from input data unit (IN) to Myr (time[Myr]=time[IN]*tscale)
- **Source**: `bse-interface/bse_interface.h`:1121

### `--mobse-vsclae`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: Velocity scale factor from input data unit(IN) to km/s (v[km/s]=v[IN]*vscale)
- **Source**: `bse-interface/bse_interface.h`:1124

### `--mobse-wdflag`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: if >0, uses WD IFMR of HPE, 1995, MNRAS, 272, 800
- **Source**: `bse-interface/bse_interface.h`:1111
- **Guards**: `!SEVN` (negated)

### `--mobse-wind`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: Binary enhanced mass loss parameter, inactive for single
- **Source**: `bse-interface/bse_interface.h`:1095
- **Guards**: `!SEVN` (negated)

### `--mobse-xi`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: wind accretion efficiency factor
- **Source**: `bse-interface/bse_interface.h`:1101
- **Guards**: `!SEVN` (negated)

### `--p`

- **Type**: `std::string`
- **Default**: `"data.par"`
- **Description**: Parameter file prefix for sse/bse; the parameter file is [prefix].bse (BSE) or [prefix].bseEmp (BSEEMP); this option should be used first before any other options; The auto-determined prefix is '[prefix of output filename].par'; in default, it is 'data.par' ('data' is the default output filename prefix)
- **Source**: `bse-interface/bse_interface.h`:1090

### `--z`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: Metallicity
- **Source**: `bse-interface/bse_interface.h`:1138

## DSM Interrupt (`--with-interrupt=dsm`)

### `--G`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: gravitational constant
- **Source**: `src/disk_star_merger.hpp`:63

### `--dsm-dt-factor`

- **Type**: `double`
- **Default**: `0.001`
- **Description**: time step factor for mass change calculation
- **Source**: `src/disk_star_merger.hpp`:65

### `--dsm-epsilon-bh`

- **Type**: `double`
- **Default**: `0.06`
- **Description**: the kenetic energy to radiation conversion efficiency of Eddington-limited accretion for BH; = 0: no accretion growth for BH
- **Source**: `src/disk_star_merger.hpp`:59

### `--dsm-epsilon-he`

- **Type**: `double`
- **Default**: `0.006`
- **Description**: helium enrichment scaling factor, used to calculate helium enrichment timescale; = 0: no helium enrichment
- **Source**: `src/disk_star_merger.hpp`:58

### `--dsm-epsilon-mdot`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: mass-change scaling factor, used to calculate mass-change timescale; = 0: no growth or mass loss
- **Source**: `src/disk_star_merger.hpp`:60

### `--dsm-he-disk`

- **Type**: `double`
- **Default**: `0.28`
- **Description**: helium fraction in the disk, used to calculate the equilbrium mass
- **Source**: `src/disk_star_merger.hpp`:56

### `--dsm-lambda0`

- **Type**: `double`
- **Default**: `0.75`
- **Description**: fraction of star's intrinsic luminosity over the Eddington luminosity without merger
- **Source**: `src/disk_star_merger.hpp`:55

### `--dsm-medd`

- **Type**: `double`
- **Default**: `253.3124306069483`
- **Description**: initial equilbrium mass of star
- **Source**: `src/disk_star_merger.hpp`:54

### `--dsm-merger-dm`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: mass loss rate for merger
- **Source**: `src/disk_star_merger.hpp`:49

### `--dsm-merger-tdelay`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: time delay for merger to increase mass
- **Source**: `src/disk_star_merger.hpp`:52

### `--dsm-new-star-mode`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: redistribute star mode, 0: no redistribute; 1: redistribute star position and velocity in random position along a circular orbit with the semi-major axis being the distance to the center; 2: redistribute star by choosing next type 3 star
- **Source**: `src/disk_star_merger.hpp`:66

### `--dsm-rstar-power`

- **Type**: `double`
- **Default**: `0.6`
- **Description**: stellar radius power index 'n', rs = s M^n
- **Source**: `src/disk_star_merger.hpp`:50

### `--dsm-rstar-scale`

- **Type**: `double`
- **Default**: `0.0046`
- **Description**: stellar radius scale 's', rs = s M^n
- **Source**: `src/disk_star_merger.hpp`:51

### `--dsm-salpeter-mass-index`

- **Type**: `double`
- **Default**: `0.0`
- **Description**: power-law index for mass-dependent Salpeter timescale, tau_eff = tau * (threshold/mass)^index when mass > threshold; = 0: no mass dependence
- **Source**: `src/disk_star_merger.hpp`:62

### `--dsm-salpeter-mass-thres`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: mass threshold multiplier for Salpeter timescale, threshold = equilbrium_mass * this value
- **Source**: `src/disk_star_merger.hpp`:61

### `--dsm-salpeter-time`

- **Type**: `double`
- **Default**: `0`
- **Description**: salpeter timescale, for the star to reach equilbrium, if zero, no stellar evolution
- **Source**: `src/disk_star_merger.hpp`:57

### `--dsm-seed-mass`

- **Type**: `double`
- **Default**: `10.0`
- **Description**: initial mass of star seed
- **Source**: `src/disk_star_merger.hpp`:53

### `--dsm-speed-of-light`

- **Type**: `double`
- **Default**: `1.0`
- **Description**: speed of light
- **Source**: `src/disk_star_merger.hpp`:64

### `--p`

- **Type**: `std::string`
- **Default**: `"data.par"`
- **Description**: Parameter file prefix for disk star merger; the parameter file is [prefix].disk_star_merger; this option should be used first before any other options; The auto-determined prefix is '[prefix of output filename].par'; in default, it is 'data.par' ('data' is the default output filename prefix)
- **Source**: `src/disk_star_merger.hpp`:67

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
- **Default**: `"data.par"`
- **Description**: Parameter file prefix for Galpy; the parameter file is [prefix].galpy; this option should be used first before any other options; The auto-determined prefix is '[prefix of output filename].par'; in default, it is 'data.par' ('data' is the default output filename prefix)
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
- **Default**: `"data.par"`
- **Description**: Parameter file prefix for gas drag; the parameter file is [prefix].gasdrag; this option should be used first before any other options; The auto-determined prefix is '[prefix of output filename].par'; in default, it is 'data.par' ('data' is the default output filename prefix)
- **Source**: `src/gas_drag.hpp`:60

## External Hard (`--with-external-hard`)

### `--ext-hard-switch`

- **Type**: `long long int`
- **Default**: `1`
- **Description**: switch of external hard force; 0: off, 1: on
- **Source**: `src/external_hard.hpp`:30

### `--p`

- **Type**: `std::string`
- **Default**: `"data.par"`
- **Description**: Parameter file prefix for external hard force; the parameter file is [prefix].exthard; this option should be used first before any other options; The auto-determined prefix is '[prefix of output filename].par'; in default, it is 'data.par' ('data' is the default output filename prefix)
- **Source**: `src/external_hard.hpp`:31


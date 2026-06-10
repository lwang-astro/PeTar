# PeTar Script Tool Inventory

These tools are installed from `install_script_tool` in `Makefile.in` and are part of the PeTar workflow surface.

## Initialization and setup

- `petar.init`
  Convert raw particle tables into PeTar input snapshots.

- `petar.find.dt`
  Search for a suitable tree time step for a snapshot and launch setup.

- `petar.update.par`
  Update legacy parameter files from old PeTar versions.

## Output management and restart support

- `petar.data.gether` (legacy)
  Legacy tool for consolidating per-rank MPI output files and splitting mixed stellar-evolution event files from old runs.
  Current PeTar (with the tmp file mechanism) automatically generates `data.snap.lst` and committed shared output files — no gathering step is needed for current runs.
  Only use this tool for old simulation runs that lack `snap.lst`.

- `petar.data.clear`
  Remove events after a specified time before restarting from an intermediate snapshot.

## Post-processing and extraction

- `petar.data.process`
  Detect binaries/multiples, compute core and Lagrangian properties, and identify escapers.

- `petar.get.object.snap`
  Extract selected objects or binaries across snapshots into time-series files.

  **Command template**:
  ```bash
  petar.get.object.snap -i <interrupt_mode> -t <external_mode> -p <out_prefix> -f <snap_format> -m id <object_id> data.snap.lst
  ```

  | Flag | Purpose | Example |
  |------|---------|---------|
  | `-i` | Interrupt mode (none/merger/bse) | `-i bse` |
  | `-t` | External mode (none/galpy/agama) | `-t galpy` |
  | `-p` | Output prefix | `-p object` |
  | `-f` | Snapshot format (origin/post) | `-f origin` |
  | `-m id <N>` | Select object by ID | `-m id 1` |

  Omit `-i`/`-t` when the solver uses no interrupt/external mode.

- `petar.format.transfer.post`
  Convert post-processed snapshot formats among ascii, binary, and npy.

  **Command template**:
  ```bash
  petar.format.transfer.post -d single -s binary -o npy -i <interrupt_mode> -t <external_mode> data.snap.last.lst
  ```

  | Flag | Purpose | Example |
  |------|---------|---------|
  | `-d` | Data kind (single/binary) | `-d single` |
  | `-s` | Source format (ascii/binary/npy) | `-s binary` |
  | `-o` | Output format (ascii/binary/npy) | `-o npy` |
  | `-i` | Interrupt mode | `-i bse` |
  | `-t` | External mode | `-t galpy` |

  The input list file contains snapshot filenames (one per line), typically `data.snap.last.lst` for a single snapshot. Mode flags `-i`/`-t` must match the producing solver family.

- `petar.galev.process`
  Generate or convert inputs for Galev workflows.

## Visualization

- `petar.movie`
  Generate movies of particle distributions, HR diagrams, binaries, and Lagrangian evolution.

- `petar.external.pot.movie`
  Generate movies for external potential map evolution.

- `petar.galpy.pot.movie`
  Generate movies for Galpy potential map evolution.

## Additional utility scripts

- `petar.get.init.binary`
  Generate initial binary lists for BSE-style initialization workflows.

- `petar.galpy.help`
  Query Galpy potential families and argument/config help.

- `petar.external.galpy` / `petar.external.agama`
  Generate external potential map snapshots from run parameter files for visualization workflows.

## Skill usage rule

When the user asks for one of these tasks, the skill should suggest the corresponding tool command, not only the main solver binary.

## External potential map workflow

When users ask to generate external potential maps before movies/overlays:

1. Create a map config file (for example `pot_conf`) describing time and grid ranges.
2. Generate map snapshots:
  - `petar.external.galpy -p data.par -m pot_conf`
  - `petar.external.agama -p data.par -m pot_conf`
3. Visualize with `petar.external.pot.movie` or overlay with `petar.movie --ext-pot`.

## Parallel sizing quick rule

Keep launch recommendations consistent with `sample/` scripts:

- Non-binaries, small-`N` (`N~10^3`) examples:
  default to single thread (`OMP_NUM_THREADS=1`) unless benchmark shows benefit.
- Binaries-rich, small-`N` (`N~10^3` with many primordial binaries) examples:
  do not hard-cap to one thread; suggest trying moderate OpenMP (for example `OMP_NUM_THREADS=2-4`) and benchmarking.
- In all cases:
  set `OMP_STACKSIZE=128M`, set `OMP_NUM_THREADS` explicitly when user asks for a concrete launch layout, and avoid leaving thread count implicit in final tuned commands.

## Snapshot reader warning rule

For Python-based readers such as `petar.movie`, `petar.data.process`, `petar.format.transfer.post`, `petar.galev.process`, and `petar.get.object.snap`:

- `Binary file size is not aligned with dtype itemsize` means binary snapshot reading is likely misconfigured.
- `The reading data shape or the number of columns mismatches the number of columns` means ascii/schema reading is likely misconfigured.
- `utf-8` decode errors usually mean binary data is being read as ascii/text.

When these appear, stop the current processing command, correct format/schema parameters, and retry before trusting outputs.

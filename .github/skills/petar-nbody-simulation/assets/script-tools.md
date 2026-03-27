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

- `petar.data.gether`
  Gather MPI outputs, split SSE/BSE event files, and generate snapshot lists.

- `petar.data.clear`
  Remove events after a specified time before restarting from an intermediate snapshot.

## Post-processing and extraction

- `petar.data.process`
  Detect binaries/multiples, compute core and Lagrangian properties, and identify escapers.

- `petar.get.object.snap`
  Extract selected objects or binaries across snapshots into time-series files.

- `petar.format.transfer.post`
  Convert post-processed snapshot formats among ascii, binary, and npy.

- `petar.galev.process`
  Generate or convert inputs for Galev workflows.

## Visualization

- `petar.movie`
  Generate movies of particle distributions, HR diagrams, binaries, and Lagrangian evolution.

- `petar.external.pot.movie`
  Generate movies for external potential map evolution.

## Skill usage rule

When the user asks for one of these tasks, the skill should suggest the corresponding tool command, not only the main solver binary.

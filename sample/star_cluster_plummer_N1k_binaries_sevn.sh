set -e

if ! command -v mcluster >/dev/null 2>&1; then
	echo "Error: mcluster is not found in PATH. Please install mcluster first." >&2
	exit 1
fi

if ! command -v petar.select >/dev/null 2>&1; then
	echo "Error: petar.select is not found in PATH. Please run 'make install' first." >&2
	exit 1
fi

if ! command -v petar.find.dt >/dev/null 2>&1; then
	echo "Error: petar.find.dt is not found in PATH. Please run 'make install' first." >&2
	exit 1
fi

# The SEVN library installer does NOT copy the lookup tables: they remain in
# the SEVN source checkout. The table location is therefore resolved after
# petar.select below, from the build-time default baked into the selected
# binary (always pointing at the SEVN source of that build), so this sample
# does not depend on any fixed local path. Override by exporting SEVN_TABLES
# (full path to a table set) or SEVN_DIR (SEVN source checkout).
# The 'gold' MIST tables are used: with SEVNtracks_MIST_AGBrobust_gold the
# usable grid at the metallicity used below is 0.7-80 Msun. The non-gold MIST
# tables have sparse coverage below ~1 Msun (stellar initialisation aborts
# there); the default PARSEC tables start at 2.2 Msun.

# SEVN runtime options used below:
# - '--tables' selects the MIST gold table set (see above).
# - '--tabuse_rhe false --tabuse_rco false --tabuse_envconv false': the MIST
#   tables ship without these optional files.
# - '-z 0.00142857' must match a metallicity directory of the table set exactly
#   (available Z: 0.0000254039, 0.0000451753, 0.0000803343, 0.000142857,
#    0.000254039, 0.000451753, 0.00142857, 0.00451753).
# - '-m 0.7 -m 80' in the mcluster call restricts stellar masses (primaries and
#   secondaries) to the table grid range: masses outside abort in stellar
#   initialisation.
# NOTE: SEVN random processes (e.g. SN kicks) are not yet seeded by PeTar, so
# runs are not bit-wise reproducible (see README, "Random numbers in stellar
# evolution").
# NOTE: '--hermite-n-neighbor-max 2000' raises the hard-integrator group
# neighbor list (default 300): with this IMF the most massive stars reach
# r_out ~ 10 pc and exceed 300 neighbors, aborting the run otherwise.
# NOTE (2026-10-04, experiment branch): primordial-binary ICs like this one can
# trigger a 'dump_large_step_h4' hard-integrator dump at t=0; this reproduces
# identically with the bse build (pre-existing branch issue, see
# doc/sevn_integration_plan.md §3-4). Without primordial binaries
# (drop '-b 0.95' and '-b 500'), the same script runs to completion.

# use mcluster to generate a star cluster with the initial condition: N=1000 Kroupa (2001) IMF, 95% binary (Kroupa 1995 a,b Sana 2012 ..., see mcluster manual)
# The initial condition for NBODY6++GPU is created with -C 5, this is used to generated initial condtion for PeTar
mcluster -N 1000 -b 0.95 -m 0.7 -m 80 -C 5 -u 1 >mc.log

# use petar.init to create initial data for petar.
# the mcluster option '-u 1' generate data in astronomical unit (Msun, pc, km/s), but petar requires a self-consistent unit of velocity: pc/Myr, '-v kms2pcmyr' will do this.
# the stellar evolution is switched on '-s sevn'
petar.init -s sevn -v kms2pcmyr -f input test.dat.10

# Switch the active petar binary family to the SEVN stellar evolution build.
# To switch on SEVN stellar evolution, the option '--with-interrupt=sevn' is needed during configuration of petar.
# Parallel hint: for N~10^3 with many primordial binaries (this sample), multi-threading can help.
# If needed, set 'OMP_NUM_THREADS=[number of threads]' (e.g. 2-4) and benchmark on your machine.
petar.select --require sevn --optional mpi,omp,avx512,avx2

# Resolve the SEVN tables location (order: SEVN_TABLES > baked default of the
# selected binary > SEVN_DIR), then verify the gold table set exists.
if [ -z "$SEVN_TABLES" ]; then
	SEVN_TABLES_DEFAULT=$(petar -h 2>&1 | awk '/^  --tables /{print $NF}')
	SEVN_TABLES="$(dirname "$SEVN_TABLES_DEFAULT")/SEVNtracks_MIST_AGBrobust_gold"
fi
if [ ! -d "$SEVN_TABLES" ] && [ -n "$SEVN_DIR" ]; then
	SEVN_TABLES="$SEVN_DIR/tables/SEVNtracks_MIST_AGBrobust_gold"
fi
if [ ! -d "$SEVN_TABLES" ]; then
	echo "Error: SEVN tables not found at '$SEVN_TABLES'." >&2
	echo "The SEVN library installer does not copy tables; export SEVN_TABLES (path to a table set) or SEVN_DIR (SEVN source directory)." >&2
	exit 1
fi
echo "SEVN tables: $SEVN_TABLES"

# Optimise tree time step with petar.find.dt for best performance.
# The SEVN runtime options are repeated via '-a' so that the test steps use the same stellar evolution setup.
# NOTE: petar.find.dt only tests the first 6 steps, so the recommended dt
# may degrade later. Always halve it for the production run.
# See SKILL.md "Performance Optimisation" for details.
dt_rec=$(petar.find.dt -a "-u 1 -b 500 -z 0.00142857 --tables $SEVN_TABLES --tabuse_rhe false --tabuse_rco false --tabuse_envconv false --hermite-n-neighbor-max 2000" -i 1 input 2>/dev/null | \
    grep "Best performance choice" | sed 's/.*tree step: //' | sed 's/,.*//')
dt_use=$(python3 -c "print(float('$dt_rec') / 2.0)")
echo "--- Production tree time step (half of recommended): $dt_use ---"

# Use PeTar to execute the simulation with SEVN stellar evolution.
# Use '-t 100.0' to run the simulation for 100 Myr.
# Use '-o 5' to generate output snapshots every 5 Myr.
# Use '-u 1' to set the units to astronomical units (Msun, pc, pc/Myr).
# Use '-b 500' to specify the number of primordial binaries as 500.
# Use '-z 0.00142857' to set the metallicity (must match a Z directory of the table set).
# Use '-s $dt_use' to use the optimised tree time step from petar.find.dt.
SEVN_ARGS="-z 0.00142857 --tables $SEVN_TABLES --tabuse_rhe false --tabuse_rco false --tabuse_envconv false --hermite-n-neighbor-max 2000"
OMP_STACKSIZE=128M petar -u 1 -b 500 $SEVN_ARGS -t 100.0 -o 5.0 -s "$dt_use" input &>output

# after mode finished, gether the output data and do post-data process to detect binaries, obtain Lagrangian and core radii and corresponding properties.
petar.data.gether data
petar.data.process -i sevn data.snap.lst

# use petar.movie to generate an animation with x-y, binary semi-ecc, and Lagrangian-radius panels.
# For SEVN cases, use '-c logtemp' to color particles by logarithmic temperature.
# - x-y panel uses '-m x-y -R 10'.
# - semi-ecc panel is enabled by '-b'.
# - Lagrangian panel uses '-L data.lagr --rlagr-min 0 --rlagr-max 5'.
petar.movie -i sevn -t none -m x-y -R 10 -c logtemp -b -L data.lagr --rlagr-min 0 --rlagr-max 5 data.snap.lst

set -e

if ! command -v mcluster >/dev/null 2>&1; then
	echo "Error: mcluster is not found in PATH. Please install mcluster first." >&2
	exit 1
fi

if ! command -v petar.select >/dev/null 2>&1; then
	echo "Error: petar.select is not found in PATH. Please run 'make install' first." >&2
	exit 1
fi

# use mcluster to generate a star cluster with the initial condtion: 
# N=1000
# Kroupa (2001) IMF
# half-mass radius of 1 pc
# a Plummer density profile
# units in astronomical units (Msun, pc, km/s)
# The option '-C 5' is used to create initial conditions for NBODY6++GPU, which will be used to generate initial conditions for PeTar.
mcluster -N 1000 -R 1 -C 5 -u 1 >mc.log

# use petar.init to create initial data for petar.
# the mcluster option '-u 1' generate data in astronomical unit (Msun, pc, km/s), but petar requires a self-consistent unit of velocity: pc/Myr, '-v kms2pcmyr' will do this.
petar.init -v kms2pcmyr -f input test.dat.10

# Use PeTar to execute the simulation.
# Use '-t 100.0' to run the simulation for 100 Myr.
# Use '-o 1.0' to generate output snapshots every 1 Myr.
# Use '-u 1' to set the units to astronomical units (Msun, pc, pc/Myr).
# Parallel hint: this is a non-binaries N~10^3 sample, where one thread is usually efficient.
# If needed, set 'OMP_NUM_THREADS=[number of threads]' and benchmark on your machine.
# set 'OMP_STACKSIZE' to ensure sufficient stack memory for each thread, otherwise segmentation faults may occur.
petar.select --optional mpi,omp,avx512,avx2
OMP_NUM_THREADS=1 OMP_STACKSIZE=128M petar -u 1 -t 100.0 -o 1.0 input &>output

# after mode finished, gether the output data and do post-data process to detect binaries, obtain Lagrangian and core radii and corresponding properties.
# To maintain consistent units during post-processing, use '-G 0.00449830997959438' to set the gravitational constant to astronomical units.
petar.data.gether data
petar.data.process -G 0.00449830997959438 data.snap.lst

# use petar.movie to generate an animation with one x-y panel and one Lagrangian-radius panel.
# - x-y panel uses '-m x-y -R 10'.
# - Lagrangian panel uses '-L data.lagr --rlagr-min 0 --rlagr-max 5'.
petar.movie -i none -t none -m x-y -R 10 -L data.lagr --rlagr-min 0 --rlagr-max 5 data.snap.lst

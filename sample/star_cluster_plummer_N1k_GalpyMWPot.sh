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
# the external potential using galpy is switched on by '-t'
# the central position and velocity of cluster '-c x[pc],y,z,vx[km/s],vy,vz' in the galactic central frame is needed
petar.init -c 8000,0,0,0,220,0 -t -v kms2pcmyr -f input test.dat.10

# Use PeTar to execute the simulation.
# Use '-t 100.0' to run the simulation for 100 Myr.
# Use '-o 1.0' to generate output snapshots every 1 Myr.
# Use '-u 1' to set the units to astronomical units (Msun, pc, pc/Myr).
# Use '--galpy-set MWPotential2014' to switch on MilkyWay potential from Galpy model (see Bovy 2015 for details).
# Parallel hint: this is a non-binaries N~10^3 sample, where one thread is usually efficient.
# If needed, set 'OMP_NUM_THREADS=[number of threads]' and benchmark on your machine.
# set 'OMP_STACKSIZE' to ensure sufficient stack memory for each thread, otherwise segmentation faults may occur.
# To switch on Galpy potential package, the option '--with-external=galpy' is needed during configuration of petar.
petar.select --require galpy --optional mpi,omp,avx512,avx2

# Optimise tree time step with petar.find.dt for best performance.
# NOTE: petar.find.dt only tests the first 6 steps, so the recommended dt
# may degrade later. Always halve it for the production run.
# See SKILL.md "Performance Optimisation" for details.
dt_rec=$(petar.find.dt -a "-u 1 --galpy-set MWPotential2014" -i 1 input 2>/dev/null | \
    grep "Best performance choice" | sed 's/.*tree step: //' | sed 's/,.*//')
dt_use=$(python3 -c "print(float('$dt_rec') / 2.0)")
echo "--- Production tree time step (half of recommended): $dt_use ---"

# Use PeTar to execute the simulation.
# '-s $dt_use' uses the optimised tree time step from petar.find.dt.
OMP_NUM_THREADS=1 OMP_STACKSIZE=128M petar -u 1 --galpy-set MWPotential2014 -t 100.0 -o 1.0 -s "$dt_use" input &>output

# after mode finished, gether the output data and do post-data process to detect binaries, obtain Lagrangian and core radii and corresponding properties.
# To maintain consistent units during post-processing, use '-G 0.00449830997959438' to set the gravitational constant to astronomical units.
# to obtain the estimation of tidal radius, using '--r-escape tidal'
# the files data.lagr, data.core, data.tidal are generated, see README of PeTar for details
petar.data.gether data
petar.data.process -t galpy --r-escape tidal -G 0.00449830997959438 data.snap.lst

# use petar.movie to generate an animation with two x-y panels and one Lagrangian-radius panel.
# - first x-y panel focuses on the cluster with '-R 10'.
# - second x-y panel shows the Galactic scale with '-R 10000'.
# - Lagrangian panel uses '-L data.lagr --rlagr-min 0 --rlagr-max 5'.
petar.movie -i none -t galpy -m x-y,x-y -R 10,10000 --cm-mode core,none --marker-scale 1,0.1 -L data.lagr --rlagr-min 0 --rlagr-max 5 data.snap.lst

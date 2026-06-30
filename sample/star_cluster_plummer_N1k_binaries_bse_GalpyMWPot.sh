set -e

# use mcluster to generate a star cluster with the initial condition: N=1000 Kroupa (2001) IMF, 95% binary (Kroupa 1995 a,b Sana 2012 ..., see mcluster manual) and place in the Milkyway potential
# The initial condition for NBODY6++GPU is created with -C 5, this is used to generated initial condtion for PeTar
if ! command -v mcluster >/dev/null 2>&1; then
    echo "Error: mcluster is not found in PATH. Please install mcluster first." >&2
    exit 1
fi

if ! command -v petar.select >/dev/null 2>&1; then
	echo "Error: petar.select is not found in PATH. Please run 'make install' for PeTar first." >&2
	exit 1
fi

if ! command -v petar.find.dt >/dev/null 2>&1; then
	echo "Error: petar.find.dt is not found in PATH. Please run 'make install' first." >&2
	exit 1
fi
mcluster -N 1000 -b 0.95 -C 5 -u 1 >mc.log

# use petar.init to create initial data for petar.
# the mcluster option '-u 1' generate data in astronomical unit (Msun, pc, km/s), but petar requires a self-consistent unit of velocity: pc/Myr, '-v kms2pcmyr' will do this.
# the stellar evolution is switched on by '-s bse' 
# the external potential using galpy is switched on by '-t'
# the central position and velocity of cluster '-c x[pc],y,z,vx[km/s],vy,vz' in the galactic central frame is needed
petar.init -c 8000,0,0,0,220,0 -t -s bse -v kms2pcmyr -f input test.dat.10

# Use PeTar to execute the simulation with stellar evolution and MilkyWay potential.
# Use '-t 100.0' to run the simulation for 100 Myr.
# Use '-o 5' to generate output snapshots every 5 Myr.
# Use '-u 1' to set the units to astronomical units (Msun, pc, pc/Myr).
# Use '-b 500' to specify the number of primordial binaries as 500.
# Use '--bse-metallicity 0.02' to set the metallicity of star as Z=0.02.
# Use '--galpy-set MWPotential2014' to switch on MilkyWay potential from Galpy model (see Bovy 2015 for details).
# To switch on SSE/BSE stellar evolution package and Galpy potential package, the option '--with-interrupt=bse --with-external=galpy' is needed during configuration of petar.
# Parallel hint: for N~10^3 with many primordial binaries (this sample), multi-threading can help.
# If needed, set 'OMP_NUM_THREADS=[number of threads]' (e.g. 2-4) and benchmark on your machine.
petar.select --require bse,galpy --optional mpi,omp,avx512,avx2

# Optimise tree time step with petar.find.dt for best performance.
# NOTE: petar.find.dt only tests the first 6 steps, so the recommended dt
# may degrade later. Always halve it for the production run.
# See SKILL.md "Performance Optimisation" for details.
dt_rec=$(petar.find.dt -a "-u 1 -b 500 --galpy-set MWPotential2014 --bse-metallicity 0.02" -i 1 input 2>/dev/null | \
    grep "Best performance choice" | sed 's/.*tree step: //' | sed 's/,.*//')
dt_use=$(python3 -c "print(float('$dt_rec') / 2.0)")
echo "--- Production tree time step (half of recommended): $dt_use ---"

# Use PeTar to execute the simulation with stellar evolution and MilkyWay potential.
# '-s $dt_use' uses the optimised tree time step from petar.find.dt.
OMP_STACKSIZE=128M petar -u 1 -b 500 --galpy-set MWPotential2014 --bse-metallicity 0.02 -t 100.0 -o 5.0 -s "$dt_use" input &>output

# after mode finished, gether the output data and do post-data process to detect binaries, obtain Lagrangian and core radii and corresponding properties.
# to obtain the estimation of tidal radius, using '--r-escape tidal'
# the files data.lagr, data.core, data.tidal are generated, see README of PeTar for details
petar.data.gether data
petar.data.process -i bse -t galpy --r-escape tidal data.snap.lst

# use petar.movie to generate an animation with two x-y panels, one binary semi-ecc panel, and one Lagrangian-radius panel.
# For BSE cases, use '-c logtemp' to color particles by logarithmic temperature.
# - first x-y panel focuses on the cluster with '-R 10'.
# - second x-y panel shows the Galactic scale with '-R 10000'.
# - semi-ecc panel is enabled by '-b'.
# - Lagrangian panel uses '-L data.lagr --rlagr-min 0 --rlagr-max 5'.
petar.movie -i bse -t galpy -m x-y,x-y -R 10,10000 --cm-mode core,none --marker-scale 1,0.1 -c logtemp,logtemp -b -L data.lagr --rlagr-min 0 --rlagr-max 5 data.snap.lst


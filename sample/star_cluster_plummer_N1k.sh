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
# By default, OpenMP utilizes all CPU threads. For small N<=1000, one CPU is sufficient, 
# use 'OMP_NUM_THREADS=[number of threads]' to limit the number of threads.
# set 'OMP_STACKSIZE' to ensure sufficient stack memory for each thread, otherwise segmentation faults may occur.
OMP_NUM_THREADS=1 OMP_STACKSIZE=128M petar -u 1 -t 100.0 -o 1.0 input &>output

# after mode finished, gether the output data and do post-data process to detect binaries, obtain Lagrangian and core radii and corresponding properties.
# To maintain consistent units during post-processing, use '-G 0.00449830997959438' to set the gravitational constant to astronomical units.
petar.data.gether data
petar.data.process -G 0.00449830997959438 data.snap.lst

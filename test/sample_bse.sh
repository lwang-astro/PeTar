# use mcluster to generate a star cluster with the initial condtion: N=1000 Kroupa (2001) IMF, 95% binary (Kroupa 1995 a,b Sana 2012 ..., see mcluster manual) 
# The initial condition for NBODY6++GPU is created with -C 5, this is used to generated initial condtion for PeTar
mcluster -N 1000 -b 0.95 -C 5 -u 1 >mc.log

# use petar.init to create initial data for petar.
# the mcluster option '-u 1' generate data in astronomical unit (Msun, pc, km/s), but petar requires a self-consistent unit of velocity: pc/Myr, '-v kms2pcmyr' will do this.
# the stellar evolution is switched on '-s bse' 
petar.init -s bse -v kms2pcmyr -f input test.dat.10

# Use PeTar to execute the simulation with stellar evolution.
# Use '-t 100.0' to run the simulation for 100 Myr.
# Use '-o 5' to generate output snapshots every 5 Myr.
# Use '-u 1' to set the units to astronomical units (Msun, pc, pc/Myr).
# Use '-b 500' to specify the number of primordial binaries as 500.
# Use '--bse-metallicity 0.02' to set the metallicity of star as Z=0.02.
# To switch on SSE/BSE stellar evolution package, the option '--with-interrupt=bse' is needed during configuration of petar.
# By default, OpenMP utilizes all CPU threads. If you wish to use a specific number of threads, please add 'OMP_NUM_THREADS=[number of threads]'."
OMP_NUM_THREADS=8 OMP_STACKSIZE=128M petar -u 1 -b 500 --bse-metallicity 0.02 -t 100.0 -o 5.0 input &>output

# after mode finished, gether the output data and do post-data process to detect binaries, obtain Lagrangian and core radii and corresponding properties.
petar.data.gether data
petar.data.process -i bse data.snap.lst

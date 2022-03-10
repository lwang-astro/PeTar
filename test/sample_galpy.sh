# use mcluster to generate a star cluster with the initial condition: N=1000 Kroupa (2001) IMF, 95% binary (Kroupa 1995 a,b Sana 2012 ..., see mcluster manual) and place in the Milkyway potential
# The initial condition for NBODY6++GPU is created with -C 5, this is used to generated initial condtion for PeTar
mcluster -N 1000 -b 0.95 -C 5 -u 1 >mc.log

# use petar.init to create initial data for petar.
# the mcluster option '-u 1' generate data in astronomical unit (Msun, pc, km/s), but petar requires a self-consistent unit of velocity: pc/Myr, '-v kms2pcmyr' will do this.
# the stellar evolution is switched on by '-s bse' 
# the external potential using galpy is switched on by '-t'
# in configure of petar, '--with-interrupt=[bse/mobse] --with-external=galpy' is needed
# the central position and velocity of cluster '-c x[pc],y,z,vx[km/s],vy,vz' in the galactic central frame is needed
petar.init -c 8000,0,0,0,220,0 -t -s bse -v kms2pcmyr -f input test.dat.10

# use petar to run the model for 100 Myr
# the OpenMP number needs to be adjusted based on number of CPU cores
# the Galactic potential MWPotential2014 from galpy is used. 
# the metallicity of stars is Z=0.02
OMP_NUM_THREADS=8 OMP_STACKSIZE=128M petar -u 1 -b 500 --galpy-set MWPotential2014 --bse-metallicity 0.02 -t 100.0 -o 5.0 input &>output

# after mode finished, gether the output data and do post-data process to detect binaries, obtain Lagrangian and core radii and corresponding properties.
# to obtain the estimation of tidal radius, using '--r-escape tidal'
petar.data.gether data
petar.data.process -i bse -t galpy --r-escape tidal data.snap.lst
# the files data.lagr, data.core, data.tidal are generated, see README of PeTar for details


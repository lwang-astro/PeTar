# P3TA
Particle-particle \& Particle-tree with Algorithmic regularization code for collisional system

# Dependence: 
FDPS: https://github.com/FDPS/FDPS
ARC: https://github.com/lwang-astro/ARC

# Current developement status:
## Schedule: 
1. soft tree_for_force with artificial particles
2. Search clusters with all artificial particles and force correction for cutoff function
3. Hard find ARC systems: because all artificial particles are inside the same cluster with the corresponding binaries, we only need to find which are binaries, which are corresponding artificial particles
4. Hard integration
5. Resolve all binaries to singles and remove all artificial particles
6. Find new binaries based on new positions of all particles in the same cluster and create new artificial particles
7. Domain + exchange particle (include artificial particles)
8. Go to 1

## Issue:
1. Energy conservation for massive body encounters
2. Unstable time step for specific ARC systems

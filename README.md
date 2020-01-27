# PeTar
Particle-particle \& Particle-tree with Algorithmic regularization code for collisional system

# Dependence: 
FDPS: https://github.com/FDPS/FDPS
SDAR: https://github.com/lwang-astro/SDAR

# Current developement status:
## Schedule: 
1. tree_for_force_short_symmetry for real particles neighbor number ( no spj) 
2. Search clusters with real particles (no force correction)
3. For each cluster, we find new binaries and create artificial particles by adding to the end of particle system locally
4. Do soft tree_for_force
5. Go to hard: find all artificial particles to corresponding binaries (all artificial particles address are registered in the binary particle information), 
6. force correction for cutoff function
7. Hard integration 
8. Remove artificial particles by setNumberOfParticalLocal
9. Performance weighted Domain + exchange particle (only real particles)
10. Go to 1

## Issue:
1. Energy conservation for tidal tensor potential
2. Some perturbation jump happen after implementation of varying changeover functions

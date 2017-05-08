# P3T+ARC
N-body with P3T for long-range force and Algorithmic regularization chain for short range force.

# Dependence: 
FDPS: https://github.com/FDPS/FDPS
ARC: https://github.com/lwang-astro/ARC

# Current developement status:

1. Long-range interaction calculation
   Using FDPS
2. Short-range interaction calculation (few-body systems)
   1. Using ARC for multi-body systems (N>2)
   2. Using kepler solver for isolated two-body systems (N=2)
   3. Using hyperbolic solver for hyperbolic encounter
	  *Not yet implemented, currently ARC is used*.
   4. **Need solution for stable triple (quadruple) systems**
3. Neighbor radius criterion
   1. Multi-body systems use outer-most binary apocenter distance 
   2. Binaries use apocenter distance
   3. Singles use r\_in
   4. If above radius is less than r\_in, use r\_in
4. Clustering method (**need to be implementeds**)
   1. Method 1:
	  1. symmetric search in FPDS is necessary
	  2. Use asymmetric radius (only j-particle depended) and modify current clustering method
		 Need to add isolated particle detection in FPDS force calculation
   2. Method 2:
	  1. symmetric search in FDPS
	  2. Use i-j depended searching radius and save largest radius for each particle

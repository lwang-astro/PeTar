### How to generate agama Miklyway Potential with flat NFW halo potential configure file 

First, Modify Agama offical example file: example_mw_potential_hunter24.py 

Change Line 90, 91 from

```python
    params_dark = dict(type='Spheroid', densitynorm=2.774e11, gamma=0, beta=0, alpha=1,
        outerCutoffRadius=8.682e-6, cutoffStrength = 0.1704)
```

to

```python
    params_dark = dict(type='Spheroid', densitynorm=1.21e7, gamma=1, beta=3, alpha=1,
        outerCutoffRadius=1e16, cutoffStrength = 5, axisRatioZ=0.84, scaleRadius=14.39)
```

Then, run example_mw_potential_hunter24.py to generate "MWPotentialHunter24_full.ini", which corresponds to the NFW model with q = 0.84 for the dark matter halo. You can directly modify the axisRatioZ parameter to change q. The densitynorm and scaleRadius values are taken from: Huang, Y., Liu, X. W., Yuan, H. B., et al. 2016, MNRAS, 463, 2623.

Subsequently, use MWPotentialHunter24_full.ini (bar + other axisymmetric or spherical components) and MWPotentialHunter24_spiral.ini (spiral components) as files, and add rotation (rotation) as needed.


Parameters:

`type='Spheroid'`: This specifies that the density profile is a spheroidal distribution, which is a common choice for modeling dark matter halos.
$$
\rho = \rho_0 \left(\frac{\tilde{r}}{a}\right)^{-\gamma} \left[1 + \left(\frac{\tilde{r}}{a}\right)^\alpha \right]^{\frac{\gamma - \beta}{\alpha}} \times \exp \left[-\left(\frac{\tilde{r}}{r_{\text{cut}}}\right)^\xi \right]
$$
When $\gamma=1$，$\alpha=1$， $\beta=3$，$q=0.84$，$r_{\text{cut}}$ uses a large value, $\xi$ is positive, let exponential part approach 1, then the above equation becomes
$$
{\rho }_{\text{halo }}\left( r\right)  = \frac{{\rho }_{0}}{\left( {r/{r}_{\mathrm{s}}}\right) {\left( 1 + r/{r}_{\mathrm{s}}\right) }^{2}}
$$
This is the NFW model. For more details, refer to page 16 of the Agama official documentation.



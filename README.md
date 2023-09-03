# Radiative transfer of Lyman-α photons at cosmic dawn with realistic gas physics
Here you will find the necessary files to run the cosmological initial conditions generator code called `MUSIC` and the namelist file for `RAMSES` to reproduce the main results of [Mittal et al (2023)](www.arxiv.com).

`21cmsig.py` is the main code that takes the gas density, temperature and neutral hydrogen fraction from `RAMSES` and p<sub>α</sub> from `RASCAS` and calculates various quantities such as Lyα coupling, 21-cm signal in 3D space. The outputs are in a 3D array. 

`T_rho.py` will plot the temperature vs density scatter plot. Also, outputs the ideal-gas-law adiabatic index. An example figure is uploaded `Tvsrho.pdf`.

`haloesinfo.py` runs a halo finder on your `RAMSES` snapshot. Then it will generate an unformatted sequential binary file to be used by the code `PhotonsFromHaloes.f90`. The first level is the number of haloes (N), the second to (N+1)th level are the positions in code units and finally, we have the box luminosity in units of number of photons per second.

`hmf.py` runs a halo finder on your RAMSES snapshot (if not done already), makes a halo mass function plot and compares it with an analytical function such as Press & Schechter (1974).

`gas_analysis.py` can be used to plot the gas density, temperature or neutral hydrogen fraction slice/projections of your `RAMSES` output. This gives more flexibility in plotting so better to use this one rather than `gas_analysis_old.py`.

`dm_analysis.py` can be used to plot the DM density slice/projections of your `RAMSES` output. It can be run on both types of simulations, hydro or DM-only. This gives more flexibility in plotting so better to use this one rather than `dm_analysis_old.py`.

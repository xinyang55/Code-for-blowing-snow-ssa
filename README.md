# Code-for-blowing-snow-ssa

Fortran code for blowing-snow sea salt aerosol (SSA) production only (no bromine).

## Citation
Please cite (always resolves to the latest version):

Xin Yang. Code for blowing-snow SSA production. Zenodo. https://doi.org/10.5281/zenodo.17735152

## Versions
- v1.1.0 (2026): https://doi.org/10.5281/zenodo.23559494
  Corrects the unit conversion on line 141 (/100 -> /1000; 1 mm/day SWE = 1/86400 kg/m²/s).
  v1.0.0 overestimated the sublimation and SSA fluxes by a factor of 10.
  This bug affects only this standalone script, not the related papers or model simulations.
- v1.0.0 (2025): https://doi.org/10.5281/zenodo.17735153 (superseded; contains the bug above)

# Code-for-blowing-snow-ssa

Fortran code for blowing-snow sea salt aerosol (SSA) production only (no bromine).

## Citation
Please cite (always resolves to the latest version):

Xin Yang. Code for blowing-snow SSA production. Zenodo. https://doi.org/10.5281/zenodo.17735152

## Versions
- v1.3.0 (2026): DOI to be added
  Optional extra snow-particle bin covering 0.3-0.5 um (centred at 0.4 um), as used in the SIP
  (ice-chip) work; switch iextra (1 = on, default; 0 = off). Adds Aitken-mode SSA
  (total SSA number x1.5-1.7 at 12 m/s) with a negligible change in SSA mass.
  With iextra = 0 the results are identical to v1.2.0.
- v1.2.0 (2026): https://doi.org/10.5281/zenodo.23570639 (superseded)
  Saturation vapour pressure (est) now calculated over ice (Jacobson, 2005, eq. 2.64),
  consistent with RH wrt ice in Dery and Yau (1999, 2001). RH input must be relative to ice.
- v1.1.0 (2026): https://doi.org/10.5281/zenodo.23559494 (superseded)
  Corrects the unit conversion (/100 -> /1000; 1 mm/day SWE = 1/86400 kg/m²/s).
  v1.0.0 overestimated the sublimation and SSA fluxes by a factor of 10.
  This bug affects only this standalone script, not the related papers or model simulations.
- v1.0.0 (2025): https://doi.org/10.5281/zenodo.17735153 (superseded; contains the bug above)

# spectre
S.P.E.C.T.R.E. (Spectral Energy Tendencies and Integrals in the Ocean &amp; Atmosphere)

A set of Python/Cython functions and scripts to calculate the transfer of kinetic/potential energy 
in the atmosphere and ocean. 
The code was used by Kjellsson & Zanna (Fluids, 2017) http://dx.doi.org/10.3390/fluids2030045
to calculate kinetic energy spectra and tendencies from a set of global ocean simulations with NEMO. 

# structure
| File name           | Usage                                                                                      |
|:-------------------:| ------------------------------------------------------------------------------------------ |
| settings.py         | User settings (simulation name, time period, regions of interest etc.)                     |
| spectral_flux.py    | Main code to run the analysis                                                              | 
| plot_stored_data.py | Script to plot results of analysis                                                         |
| psd_functions.py    | Functions for spectral transforms, e.g. setup wavenumbers, or calculate spectral fluxes.   |
| plot_functions_2.py | Functions to plot results                                                                  |
| read_stored_data.py | Functions to read the results of analysis                                                  |

# note of caution
The code for available potential energy (APE) and vertical energy transfers are still experimental. 
The methods for evaluating the momentum tendencies are specific to NEMO (e.g. vector invariant form etc.). They may need to be changed for other model outputs. 
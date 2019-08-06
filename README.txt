GalDNA is an IDL code that uses an MCMC algorithm to sample 4-dimensional parameter space (log(Z_neb/Z_sun), log(N/O), log(U), and log(Z_star/Z_sun)) and compare predictions from a photoionization model grid to meaurements of common nebular emission lines. As output, GalDNA produces a text file with the full MCMC chains that can then be analyzed as desired.

The 'pro' directory contains the main GalDNA code, as well as routines that GalDNA needs to run. The GSFC astronomy routines from the IDL Astronomy User's Library (https://idlastro.gsfc.nasa.gov/) also need to be in the IDL path for GalDNA to run.

The 'example' directory contains a simple code that shows how to use GalDNA ('galdna_example.pro'), along with a default photoionization model grid ('cloudy_strom2018.fits') and measurements from an example galaxy ('kbss_lm1_lines.fits').

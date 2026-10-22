extend code for AOD/AAOD in
/work/noaa/wrf-chem/pagowski/jedi/code/jedi-bundle/crtm/
branch feature/aodfm_mzp
to calculate fine-mode AOD/AAOD
in the following way:

1. to functions  CRTM_AOD, CRTM_AOD_TL, CRTM_AOD_AD, CRTM_AOD_K
introduce optional logical argument FINE_MODE

2. When this argument FINE_MODE is present and TRUE
modify definition of AOD and AAOD to include only fractions
of dust bins determined by factors fmass1 in table

/work/noaa/wrf-chem/pagowski/gmao_aero_tables/optics_DU.v2.2.1.legacy.nc4

of seasalt bins determined by factors fmass1 in table

/work/noaa/wrf-chem/pagowski/gmao_aero_tables/optics_SS.v2.2.1.legacy.nc4

and of nitrate bins determined by factors fmass1 in table

/work/noaa/wrf-chem/pagowski/gmao_aero_tables/optics_NI.v2.2.1.legacy.nc4

fmass1(bin, rh) is a function of bin number and relative humidity rh

That means for calculating AOD/AAOD for a given bin, use
relative humidity at given point (i,j,k) and determine the actual value
of the factor fmass1 by interpolating in rh from the tables for dust, seasalt, and nitrate 

Then multiply mixing ratio of dust/seasalt/nitrate in bin 1 by this factor
for bin 1, multiply mixing ratio of dust/seasalt/nitrate in bin 2 by this factor
for bin 2, etc. for all given bins.

For scaling other species without bins:
BC,OC, and sulfate with fmass1 use tables
optics_BC.v2.2.1.legacy.nc4
optics_OC.v2.2.1.legacy.nc4
optics_SU.v2.2.1.legacy.nc4

Make corresponding modifications for K-matrix , tangent linear, adjoint:
CRTM_AOD_K, CRTM_AOD_TL, CRTM_AOD_AD


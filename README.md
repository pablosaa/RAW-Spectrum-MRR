# RAW-Spectrum-MRR
_Micro-Rain-Radar processor by using the Raw spectral reflectivity power_

A stand-alone program written mainly in Fortran 90 and C, with interface as MEX function for _MATLAB_.
The program optionally can be ran from _MATLAB_ workspace when the compiled MEX function is used, or can be ran from the Linux command line to directly generate _NetCDF_ archives when the stand-alone version is used.

The program uses the Micro-Rain-Radar (MRR) RAW data archives to estimate a set of precipitating variables like reflectivity, liquid water content, rain rates, and other moments related to the drop-size-distribution.
The code is specially usefull in cases when strong convective precipitation is observed,  it is when the MRR firmware cannot resolve properly issues with the spectral power, e.g. misestimated noise level, folding and aliasing. When this issues happpen,  

Within the project ADMIRARI, a Micro Rain Radar (MRR) has been used to enrich the radiometer m
easurements with a reflectivity profile, mainly with a slant configuration. For strong convect
ive cases the MRR Firmware is troublesome since the spectrum suffers from  Therefore only the MRR raw data is meaningful to be processed 
from scratch. This code not only solves the problem of aliasing and folding but also improves 
the spectrum noise level estimation thereafter the instrument sensitivity.

# RAW-Spectrum-MRR
_Micro-Rain-Radar processor by using the Raw spectral reflectivity power_

A stand-alone program written mainly in Fortran 90 and C, with interface as MEX function for _MATLAB_.
The program optionally can be ran from _MATLAB_ workspace when the compiled MEX function is used, or can be ran from the Linux command line to directly generate _NetCDF_ archives when the stand-alone version is used.

The program uses the Micro-Rain-Radar (MRR) RAW data archives to estimate a set of precipitating variables like reflectivity, liquid water content, rain rates, and other moments related to the drop-size-distribution.
The code is specially usefull in cases when strong convective precipitation is observed, it is when the MRR firmware has shown poor performance to resolve properly issues with the spectral power, e.g. misestimated noise level, folding and aliasing. When this issues happpen all the radar variables need to be re-processed from scratch by exploiting the MRR RAW data archives. The same is true when the MRR is used -unconventionally- at different configuration as the conventional vertical pointing. 

This code was developed within the project ADMIRARI, where a Micro Rain Radar (MRR) has been used to enrich the radiometer measurements with a reflectivity profile. Both instruments measured in sinergy and with a set-up configuration for slant field of view.
This code not only solves the problem of aliasing and folding but also improves the spectrum noise level estimation thereafter the instrument sensitivity.

The code gives also the possibility to use a set of input parameters in order to tune the calculation of radar variables suited for liquid of solid precipitation.


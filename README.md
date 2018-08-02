# RAW-Spectrum-MRR
_Micro-Rain-Radar processor by using the Raw spectral reflectivity power_

## Description
A stand-alone program written mainly in Fortran 90 and C, with functional interface as MEX function for _MATLAB_ and _GNU/Octave_.

The program optionally can be used as a _MATLAB_ function from workspace when the compiled MEX function is available, or the stand-alone version can be used directly from the Linux command line to process the data and generate _NetCDF_ archives.

The program is a processor for the Micro-Rain-Radar (MRR) RAW data archives and its goal is to estimate a set of precipitating variables like reflectivity, liquid water content, rain rates, and other spectrum related moments as the drop-size-distribution, spectral falling velocity.

The code is specially usefull in cases when:
* strong convective precipitation is observed, where the MRR firmware has shown poor performance to resolve properly issues related with the spectral power, e.g. misestimated noise level, folding and aliasing. When that happpens all the radar variables need to be re-processed from scratch by exploiting the unprocessed MRR RAW data archives.

* higher sensitivity is needed, e.g. shallow rain, virga detection, snowfall, etc., by improving the instrument noise level estimation thereafter the instrument's minimum detection.

* In some special cases, when the MRR is used -unconventionally- at different configuration unlike the standard vertical pointing. 

Additional capabilities have been added to the processing methods, like the possibility to use a set of alternative input parameters in order to tune the calculation of radar variables to suit specific observations like liquid of solid precipitation.

# Code Structure
The following is the tree of folders where modules, object files and source files are located:

    RAW-Spectrum-MRR
    ├── RAW-Spectrum-MRR/input
    ├── RAW-Spectrum-MRR/LICENSE
    ├── RAW-Spectrum-MRR/Makefile
    ├── RAW-Spectrum-MRR/mod
    │   ├── RAW-Spectrum-MRR/mod/extras.mod
    │   ├── RAW-Spectrum-MRR/mod/psg_nc.mod
    │   └── RAW-Spectrum-MRR/mod/variables.mod
    ├── RAW-Spectrum-MRR/MRR4ADMI
    ├── RAW-Spectrum-MRR/obj
    │   ├── RAW-Spectrum-MRR/obj/extras.o
    │   ├── RAW-Spectrum-MRR/obj/mrrmain_tux.o
    │   ├── RAW-Spectrum-MRR/obj/psg_netcdf.o
    │   ├── RAW-Spectrum-MRR/obj/read_RAW_Cfile.o
    │   └── RAW-Spectrum-MRR/obj/var_init.o
    ├── RAW-Spectrum-MRR/README.md
    └── RAW-Spectrum-MRR/src
        ├── RAW-Spectrum-MRR/src/extras.F90
        ├── RAW-Spectrum-MRR/src/mrrmain_tux.F90
        ├── RAW-Spectrum-MRR/src/psg_netcdf.F90
        ├── RAW-Spectrum-MRR/src/read_RAW_Cfile.c
        └── RAW-Spectrum-MRR/src/var_init.F90


# Compilation
To compile with GNU compilers (i.e. g++ and gfortran), use the Makefile provided
    
    user@linux:~/MIUB/scripts/RAW-Spectrum-MRR> make
    /usr/bin/gcc -Wall -O3 -c src/read_RAW_Cfile.c -o obj/read_RAW_Cfile.o
    /usr/bin/gfortran obj/mrrmain_tux.o obj/read_RAW_Cfile.o obj/psg_netcdf.o obj/extras.o obj/var_init.o -o MRR4ADMI -L/usr/lib64 -lnetcdff

this will create an executable names `MRR4ADMI` which can be run as shown below.
In case other compilers are used or other platform, the `Makefile` needs to be adapted according to the system, mainly the folloring variables in the `Makefile` are subject to be changed:

* `FCC = /usr/bin/gfortran` -> full path to FORTRAN 90 compiler
* `GCC = /usr/bin/gcc`      -> full path to C/C++ compiler
* `LIBS = /usr/lib64`       -> full path to NetCDF libraries and modules
* `INCS = /usr/include`     -> full path to NetCDF header files
* `EXEC = MRR4ADMI`         -> (optional) name of the final executable file

# Usage
In Linux terminal:

    > ./MRR4ADMI FULL_PATH_TO_INPUT_FILE/InputFileMRR.raw [OUTPUT_PATH_NETCDF]


# Examples
The following example process the RAW data file `RawSpectra-20160601.raw` and after processing the output NetCDF file is storaged at the same folder as the input file, since the second option is onmited

    user@linux:~/MIUB/scripts/RAW-Spectrum-MRR> ./MRR4ADMI ../../data/mrr1/RawSpectra-20160601.raw 
    Working at ../../data/mrr1/RawSpectra-20160601.raw
    Reflectivity calculated
    Creating NetCDF file: ../../data/mrr1/RawSpectra-20160601.nc
    NetCDF file closed succesfully


# Output

# Acknowledges
This code has been developed under the _Deutsche Forschungsgemainsam_ project ADMIRARI, where a Micro Rain Radar (MRR) has been used to enrich the radiometer retrieval capavilities with a reflectivity profile. Both instruments have been measuring in sinergy and with a set-up configuration for slant viewing angle, e.g. 30° elevation.

See LICENSE.TXT

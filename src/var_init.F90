! ===============================================================
! Module containing the general variables characterizing the MRR
! The variables are based for the MRR-2 last version (2008)
! For older version the variables could be adjusted.
! ===============================================================
! This module is part of the program MRR4ADMI to read the RAW data from the MRR
! installed aside ADMIRARI radiometer.
!
!   Copyright 2011-2014 Pablo Saavedra Garfias (pablosaa@uni-bonn.de)
!   See LICENSE.TXT
!

module variables
  use, intrinsic :: iso_c_binding, only: c_int
  implicit none
  ! ++++++++++++++++++++++++++++++++++++++++
  ! MRR data file specific parameters:
  ! ++++++++++++++++++++++++++++++++++++++++
  ! Dimentional parameters:
  integer, parameter :: N_Dh = 32
  integer, parameter :: Nspc = 64
  integer, parameter :: NData = 9000
  integer, parameter :: MAXLINE = 295
  ! Number of levels considered in the RAW data file:
  integer(kind=c_int), bind(c, name='Nh') :: Nh = N_Dh
  ! Number of spectral bins in the RAW data:
  integer(kind=c_int), bind(c, name='Nspec') :: Nspec = Nspc
  ! Max number of profiles per day: (at 10sec interval, max data is 8640)
  integer(kind=c_int), bind(c, name='Nt') :: Nt = NData
  ! size of every profile element in file:
  integer(kind=c_int), bind(c, name='PFIELD') :: PFIELD = 9
  ! size of every header in row (old version=6, new ver=3)
  integer(kind=c_int), bind(c, name='PHEAD') :: PHEAD=6
  ! Max number of characters in a like (at least PHEAD+PFIELD*Nh):
  integer(kind=c_int), bind(c, name='NLINE') :: NLINE = MAXLINE

  ! Average MRR receiver Noise level: 
  !(Warning! this is a specific table for the instrument SN:xxxxxx,
  ! do not use as a general parameter. If unkown, set all the values to zero
  ! and the result will be calculated taking into account only the noise level
  ! estimation algorithm.)
!!$  real, dimension(N_Dh), parameter :: ANL = 0.0;

!! The following is the average noise level sampled from CABAUW:
  real, dimension(N_Dh) :: ANL = &
       (/2.8488,5.0059,9.8412,12.7857,14.3621,&
         15.0165,15.2412,&
         15.2537,15.0547,14.7303,14.3248,13.8327,13.2462,12.6714,&
         12.0964,11.5600,11.0181,10.3279,9.7308,9.1801,8.6684,&
          8.1300,7.5689,7.0221,6.5785,6.0946,5.6615,5.3373,4.9115,&
          4.5105,3.7008,2.1814/)
!! The folloing is the average noise level sampled from CHUVA:
!!$  real, dimension(N_Dh), parameter :: ANL = &
!!$       (/2.8503, 4.6175, 9.6915, 12.7419, 14.2210, 14.8260, 14.8908, 14.7491,&
!!$       14.4107,13.9662,13.4251,12.8595,12.2717,11.6751,11.0953,10.5369,&
!!$       9.9602,9.3474,8.8178,8.2790,7.7683,7.2825,6.8066, 6.3597,&
!!$       5.9421, 5.5258, 5.1627, 4.7732, 4.4356, 4.0769, 3.4600, 1.9814/)

  ! ++++++++++++++++++++++++++++++++++++++++
  ! MRR working Frequency: 
  real :: MRR_freq = 24.1E9   ! Hz
  ! Kappa2 = |m^2-1/m^2+1|^2  where m is complex refractive index
  ! for water Kappa2 = 0.92, for ice Kappa2 = 0.18
  real, parameter :: Kappa2 = 0.92
  real, parameter :: Kappa_ice = 0.18
  ! +++++++++++++++++++++++++++++++++++++++++
  ! Physics and Math constants
  real, parameter :: PI = 3.141592
  ! speed of light in vacuum (def) m/s
  real, parameter :: c_li = 0.299792458e9
  ! water density [gr/cm^3]
  real, parameter :: RHOw = 1.0

  ! Backscattering coefficient estimated based on MIE calculations for MRR_freq
  ! and spherical raindrops with diameter Dn. With sigma_B having Nspec elements:
  real :: sigma_B(Nspc) = &
       &(/1.928934e-08, 9.265726e-08, 3.310414e-07, 9.674124e-07, 2.442663e-06,&
       & 5.527846e-06, 1.151002e-05, 2.245290e-05, 4.144296e-05, 7.313812e-05,&
       & 1.242299e-04, 2.043987e-04, 3.274050e-04, 5.115989e-04, 7.834597e-04,&
       & 1.178168e-03, 1.742847e-03, 2.544000e-03, 3.671177e-03, 5.242944e-03,&
       & 7.418144e-03, 1.042946e-02, 1.459371e-02, 2.034690e-02, 2.833417e-02,&
       & 3.950635e-02, 5.525775e-02, 7.770427e-02, 1.100831e-01, 1.572617e-01,&
       & 2.261786e-01, 3.270029e-01, 4.733466e-01, 6.830392e-01, 9.791517e-01,&
       & 1.390520e+00, 1.953149e+00, 2.712240e+00, 3.723007e+00, 5.055644e+00,&
       & 6.793417e+00, 9.041031e+00, 1.192839e+01, 1.558110e+01, 2.006057e+01,&
       & 2.519536e+01, 3.032725e+01, 3.401733e+01, 3.370516e+01, 2.482624e+01,&
       & 5.270436e+00, 6.869670e+01, 6.869670e+01, 6.869670e+01, 6.869670e+01,&
       & 6.869670e+01, 6.869670e+01, 6.869670e+01, 6.869670e+01, 6.869670e+01,&
       & 6.869670e+01, 6.869670e+01, 6.869670e+01, 6.869670e+01/)

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Precipitation Terminal Fall speed model:
  ! Vnmodel is a character string to specify the model to use
  ! in order to estimate the particle size from the falling velocity.
  ! It could be a model for raindrops, snow particles or hail.
  ! Options:
  !     * 'metek' as eq. (2.5) by METEK Physical Basis
  !     * ''
  character(len=15) :: Vnmodel = 'metek'  ! default (see Metek documentation)
  
end module variables

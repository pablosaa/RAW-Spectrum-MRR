! ------------------------------------------------------------------------------------
! This is a program to process the RAW data files from a MicroRainRadar (MRR) by METEK
! This is the main FORTRAN code 'mrrmain_tux.F90' for the program named MRR4ADMI which
! is hosted at the github.com/pablosaa/RAW-MRR-spectrum repository.
!
! USAGE:
! > ./MRR4ADMI /full_path_to_input_mrr/raw_datafile.mrr [/full_path/output_netcdf_file/]
! WHERE:
! * first argument is a MRR input RAW data file including full path,
! * second argument is OPTIONAL indicating the full path for the NetCDF output file.
!
! > ./MRR4ADMI /full_path_to_input_mrr/raw_datafile.mrr
! *  if second argument not specified, then the same path as the input file is used.
!
! > ./MRR4ADMI
! * When no argumed is given, then the software ask first for a full path input RAW file,
!   followed by a request of the full path for the processed output NetCDF file.
!
! In all cases the NetCDF output file has the same name as the input file but with the
! extension changed from .mrr to .nc
!
!   Copyright 2011-2014 Pablo Saavedra Garfias (pablosaa@uni-bonn.de)
!   See LICENSE.TXT
!

program mrrmain_tux
  use variables, only: MAXLINE, Nt, N_Dh, Nspc, NData, ANL, MRR_freq, sigma_B, Vnmodel
  use extras, only : calc_ref_Z, DSD2rainparam
  use psg_nc, only : save_db_netcdf
  implicit none
  namelist/varlist/Nt, ANL, MRR_freq, sigma_B, Vnmodel
  integer :: i
  character(len=MAXLINE) :: fileiname
  character(len=MAXLINE) :: output_path
  integer :: Ntot, CC, status
  real, dimension(NData) :: t_hr
  integer, dimension(N_Dh) :: height
  integer, dimension(6,2) :: dates
  real, dimension(N_Dh) :: TF
  integer, dimension(N_Dh,Nspc,NData) :: Fnn
  real, dimension(Nspc) :: Dn
  real, dimension(N_Dh, Nspc) :: VDn
  real, allocatable, dimension(:,:,:) :: DSD
  real, allocatable, dimension(:,:) :: Etha_NL, SNR, ZdBz
  real, allocatable, dimension(:,:) :: rr, lwc, Dm, Nw
  
  external read_raw_cfile

  ! Cheking command line input parameters:
  if(iargc()>2) then
     print *,"Command Line input argument not valid!"
     print *,"USAGE: >> ./MRR4raw mrr_input.raw [/output_path/]"
     print *,"The NetCDF output file will be: [/output_path/]mrr_input.nc"
     print *,"where [/output_path/] is optional."
     stop ": Try again! ;)"
  end if
  if(iargc()==1.OR.iargc()==2) then
     call getarg(1,fileiname)
  else
     write(*,'(A)',ADVANCE='NO') 'Please introduce full path to MRR RAW input data file: '
     read(*,*) fileiname
  end if

  ! checking for alternative output PATH as second argument input:
  if(iargc().eq.2) then 
     call getarg(2,output_path)
  elseif(iargc().eq.1) then
     i = index(fileiname,'/',.true.)
     output_path = fileiname(:i)
  else
     write(*,'(A)',ADVANCE='NO') 'Introduce full path for NetCDF output file: '
     read(*,*) output_path
  end if

  ! Checking if 'input' file contains alternative parameters:
  ! 'input' file should include the NAMELIST varlist as
  ! defined by: namelist/varlist/Nt, ANL, MRR_freq, sigma_B, Vnmodel
  ! where: Nt max number of profiles in a RAW data file (default: 9000)
  !        ANL profile of noise level estimation (instrument dependent)
  !        MRR_freq Operational frequency (Normally 24.1GHz)
  !        sigma_B Backscattering factor (default: spherical liquid particles)
  !        Vnmodel Terminal velocity speed model (default: METEK manual)
  open(100,FILE='input',STATUS='OLD',IOSTAT=status)
  if(status==0) then
     read(100,nml=varlist)
     close(100)
  else
     print*,'No alternative parameters have been read from input'
  end if

  ! Calling the subroutine to read the RAW data file:
  call read_raw_cfile(trim(fileiname)//char(0), Ntot, t_hr, height, Fnn, TF, CC, dates, status)
  if(status.ne.0) then
     stop "ERROR reading input file! Are you sure it is a RAW MRR data file?"
  end if

  ! Assigning dimentions to the precipitation variables:
  allocate(ZdBz(N_Dh,Ntot))   ! Radar reflactivity
  allocate(Etha_NL(N_Dh,Ntot))   ! Noise Level Spectral Reflectivity
  allocate(SNR(N_Dh,Ntot))       ! Signal-to-noise-ratio
  allocate(DSD(N_Dh,Nspc,Ntot))  ! Particle size distribution

  ! Computing Reflectivities by means of the algorithm provieded by METEK
  call calc_ref_Z(Ntot,Fnn,TF,height,CC,ZdBz,Etha_NL,SNR,Dn,DSD,VDn)

  ! Computing Rain parameters based on the DSD estimated from spectrum
  allocate(rr(Ntot,N_Dh),lwc(Ntot,N_Dh),Dm(Ntot,N_Dh),Nw(Ntot,N_Dh))
  
  do i=1,N_Dh  ! loop over the range bins
     call DSD2rainparam(Dn(:50),transpose(DSD(i,:50,:Ntot)),&
          &spread(VDn(i,:50),1,Ntot),&
          &LWC=lwc(:,i),RR=rr(:,i),Dm=Dm(:,i),Nw=Nw(:,i))
  end do

  ! Storing parameters in a NetCDF file
  call save_DB_netcdf(fileiname,output_path,dates,t_hr(:Ntot),real(height),Dn,&
       &ZdBz(:,:Ntot),transpose(rr),transpose(lwc),transpose(Dm),&
       &transpose(Nw),DSD)
  
  ! Free memory allocations:
  deallocate(DSD)
  deallocate(SNR)
  deallocate(Etha_NL)
  deallocate(ZdBz)
  deallocate(rr, lwc, Dm, Nw)
  stop
end program mrrmain_tux

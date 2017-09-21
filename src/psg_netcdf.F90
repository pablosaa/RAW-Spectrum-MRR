! ------------------------------------------------------------------------------------
! This is part of the program to process the RAW data files from a MicroRainRadar (MRR) by METEK
! This are subroutines to manage NetCDF archiving for the program named MRR4ADMI which
! is hosted at the github.com/pablosaa/RAW-MRR-spectrum repository.
!
!   Copyright 2011-2014 Pablo Saavedra Garfias (pablosaa@uni-bonn.de)
!   See LICENSE.TXT

! -------------------------------------------------------------
! subroutines to create and make a netcdf (general purposes)
!
! Pablo Saavedra, 2014. MIUB - University of Bonn
! ------------------------------------------------------------
module psg_nc
  use netcdf
  implicit none
contains
subroutine create_nc(fname,DIM,DIMVAR,nc,dim_id)

  implicit none
  character(len=*), intent(in) :: fname
  integer, dimension(:), intent(in) :: DIM
  character(len=*), dimension(:), intent(in) :: DIMVAR
  integer, intent(out) :: nc
  integer, dimension(:), intent(out) :: dim_id
  
  ! local variables
  integer :: Ndim, i

  Ndim = size(DIM)
  print*, 'Creating NetCDF file: '//trim(fname)
  call check(nf90_create(fname,nf90_clobber, nc))
  ! define dimensions:
  do i=1,Ndim
     call check(nf90_def_dim(nc,DIMVAR(i),DIM(i),dim_id(i)))
  end do
  return
  ! define variables:

end subroutine create_nc

! Subroutine for variables definition
subroutine def_vars(nc,dim_id,VAR,ATTNAME,ATTVAL,var_id)
  implicit none
  integer, intent(in) :: nc
  integer, dimension(:), intent(in) :: dim_id
  character(len=*), intent(in) :: VAR
  character(len=*), dimension(:), intent(in) :: ATTNAME
  character(len=*), dimension(:), intent(in) :: ATTVAL
  integer, intent(out) :: var_id

  integer :: i, Nval

  Nval = size(ATTNAME)

  ! In case the VARIABLE string is 'globals', then the global attributes
  ! of the netcdf will be set up, 
  if(VAR.ne.'globals') then
     call check(nf90_def_var(nc,VAR,NF90_REAL,dim_id,var_id))
  else
     var_id = NF90_GLOBAL
  end if
  do i=1,Nval
     if(ATTNAME(i).eq.'_FillValue') then
        call check(nf90_put_att(nc,var_id,ATTNAME(i),NF90_FILL_REAL))
     else
        call check(nf90_put_att(nc,var_id,ATTNAME(i),ATTVAL(i)))
     end if
  end do

  if(VAR.eq.'globals') call check(nf90_enddef(nc))
  return

end subroutine def_vars

! Subroutine to put the values of dimensions and variables:
subroutine put_vars(nc,var_id,VAR1D,shapes,LAST) !VAR2D,VAR3D,VAR4D,LAST)
  implicit none
  integer, intent(in) :: nc, var_id
  integer, dimension(:) :: shapes
  real, dimension(shapes(1),shapes(2),shapes(3)), intent(in) :: VAR1D
!!$  real, optional, dimension(:,:), intent(in) :: VAR2D
!!$  real, optional, dimension(:,:,:), intent(in) :: VAR3D
!!$  real, optional, dimension(:,:,:,:), intent(in) :: VAR4D
  logical, optional, intent(in) :: LAST


  call check(nf90_put_var(nc,var_id,VAR1D))


  if (present(LAST)) then
     if (LAST) then
        call check(nf90_close(nc))
        print*, 'NetCDF file closed succesfully'
     end if
  end if
  return
end subroutine put_vars

! Subroutine to check status of netcdf procedure:
subroutine check(status)

  integer, intent ( in) :: status
    
  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop "Stopped"
  end if
end subroutine check

  ! ===================================================
  ! Subroutine to storage the Moments in a NetCDF file
  ! ===================================================
  subroutine save_DB_netcdf(filename,output_path,the_date,T,H,DEQ,Z,RR,LWC,Dm,Nw,DSD)
    !use psg_nc
    use variables, only: MRR_freq, Vnmodel
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: output_path
    integer, intent(in), dimension(:,:) :: the_date
    real, intent(in), dimension(:) :: T, H, DEQ
    real, intent(inout), dimension(:,:) :: Z
    real, intent(in), dimension(:,:) :: LWC, RR, Dm, Nw
    real, intent(in), dimension(:,:,:) :: DSD
    integer, parameter :: Nvar=6
    integer :: Nt, Nh,Nd
    integer :: i, j, nc, dim_id(3), DIM(3), var_id(Nvar), tmp_id
    character(len=300) :: ncfile
    character(len=6), dimension(3) :: DIMVAR     ! dimension names
    character(len=15), dimension(8) :: ATTNAME   ! attribute names
    character(len=45), dimension(8) :: ATTVAL    ! attribute values

    character(len=3), dimension(Nvar) :: MOMS=(/'Z  ','RR ','LWC','Dm ','Nw ','DSD'/)

    ! creating the NetCDF output file
    i = index(filename,'/',.true.)+1
    j = index(filename,'.',.true.)
    ncfile = trim(output_path)//trim(filename(i:j))//'nc'
    
    Nt = size(T)  ! this is time e.g. unlimitted
    Nh = size(H)
    Nd = size(DEQ)
    where(Z.EQ.-99.) Z = NF90_FILL_REAL
    
    ! end of defining variables
    DIM=(/Nd,Nh,Nt/) !NF90_UNLIMITED/)
    DIMVAR = (/'DEQ   ','height','time  '/)  ! inverse orden for NetCDF

    call create_nc(ncfile,DIM,DIMVAR,nc,dim_id)

    ! defining dimensions
    ATTNAME(1) = 'long_name'
    ATTNAME(2) = 'short_name'
    ATTNAME(3) = 'units'
    ATTNAME(4) = '_FillValue'

    do i=1,size(DIMVAR)
       select case (trim(DIMVAR(i)))
       case ('time')
          ATTVAL(1) = 'time of the day'
          ATTVAL(2) = 'T'
          ATTVAL(3) = 'hours'
       case ('DEQ')
          ATTVAL(1) = 'Equivalent diameter'
          ATTVAL(2) = 'DEQ'
          ATTVAL(3) = 'milimeter'
       case ('height')
          ATTVAL(1) = 'Height a.g.l.'
          ATTVAL(2) = 'H'
          ATTVAL(3) = 'meter'
       end select
       call def_vars(nc,dim_id(i:i),DIMVAR(i),ATTNAME(1:3),&
            &ATTVAL(1:3),tmp_id)
    end do
    ! defining variables
    do i=1,Nvar
       select case(trim(MOMS(i)))
       case ('LWC')
          ATTVAL(1) = 'Liquid Water Content'
          ATTVAL(2) = 'LWC'
          ATTVAL(3) = 'g/m^3'
          ATTVAL(4) = '-99.'
       case ('RR')
          ATTVAL(1) = 'Rainrate'
          ATTVAL(2) = 'RR'
          ATTVAL(3) = 'mm/hr'
          ATTVAL(4) = '-99.'
       case ('Dm')
          ATTVAL(1) = 'Mean mass weithed diameter'
          ATTVAL(2) = 'Dm'
          ATTVAL(3) = 'mm'
          ATTVAL(4) = '-99.'
       case ('Nw')
          ATTVAL(1) = 'DSD intercept parameter?'
          ATTVAL(2) = 'Nw'
          ATTVAL(3) = 'mm^-1 m^-3'
          ATTVAL(4) = '-99.'
       case ('Z')
          ATTVAL(1) = 'Reflectivity'
          ATTVAL(2) = 'Z'
          ATTVAL(3) = 'dBZ'
          ATTVAL(4) = '-99'
       case ('DSD')
          ATTVAL(1) = 'Particle size distribution'
          ATTVAL(2) = 'DSD'
          ATTVAL(3) = 'mm^-1 m^-3'
          ATTVAL(4) = '-99.'
          call def_vars(nc,(/dim_id(2),dim_id(1),dim_id(3)/),&
               &trim(MOMS(i)),ATTNAME(1:4),&
               &ATTVAL(1:4),var_id(i))
       case default
          print*,'Problem???  '//MOMS(i)
          stop 'NetCDF variable not assigned at definiton!'        
       end select
       if(i/=6)  call def_vars(nc,dim_id(2:3),trim(MOMS(i)),ATTNAME(1:4),&
            &ATTVAL(1:4),var_id(i))
    end do

    ! Global attributes:
    ATTNAME(1) = 'Disdrometer'
    ATTNAME(2) = 'Hydrometeor'
    ATTNAME(3) = 'Frequency_GHz'  
    ATTNAME(4) = 'Source_RAWFile'
    ATTNAME(5) = 'Fall_Vel_model'
    ATTNAME(6) = 'Starting_Date'
    ATTNAME(7) = 'Ending_Date'
    ATTNAME(8) = 'Contact'
    ATTVAL(1) = 'METEK MRR'
    ATTVAL(2) = 'Water'
    write(ATTVAL(3),'(F6.2)') MRR_freq*1E-9
    ATTVAL(4) = trim(filename)
    ATTVAL(5) = Vnmodel
    write(ATTVAL(6),'(I2.2"."I2.2".20"I2.2"_"I2.2":"I2.2":"I2.2"UTC")') the_date(:,1)
    write(ATTVAL(7),'(I2.2"."I2.2".20"I2.2"_"I2.2":"I2.2":"I2.2"UTC")') the_date(:,2)
    ATTVAL(8) = 'P. Saavedra G.(pablosaa@uni-bonn.de)'
  
    call def_vars(nc,(/0/),'globals',ATTNAME,ATTVAL,tmp_id)
    ! END for variable definitions

    ! Putting the values to dimensions:
    call put_vars(nc,dim_id(1),DEQ,(/DIM(1),0,0/))
    call put_vars(nc,dim_id(2),H,(/DIM(2),0,0/))
    call put_vars(nc,dim_id(3),T,(/DIM(3),0,0/))
    
    ! Putting the values to variables:
    call put_vars(nc,var_id(1),Z,DIM(2:3))
    call put_vars(nc,var_id(2),RR,DIM(2:3))
    call put_vars(nc,var_id(3),LWC,DIM(2:3))
    call put_vars(nc,var_id(4),Dm,DIM(2:3))  ! correct dimemsions for number of data
    call put_vars(nc,var_id(5),Nw,DIM(2:3))
    call put_vars(nc,var_id(6),DSD,(/DIM(2),DIM(1),DIM(3)/),LAST=.true.)

  end subroutine save_DB_netcdf


end module psg_nc


module extras
  use variables, only : Nh, Nspec, Nt
  implicit none
contains

!     ****************************************************************
!     SUBROUTINE CALC_REF_Z
!     INPUT:
!     Fnn   = Raw Spectral power (range,spectrum,number of profiles=time).
!     TF    = Tranfer Function. (range)
!     HEIGHT= Range altitude, [m] (range)
!     CC    = Calibration Constant.
!     OUTPUT:
!     ZdBz = Attenuated reflectivity as a function of (ranges,time) [dBz]
!     Etha_NL= Spectral estimated Noise level (ranges,time).
!     SNR   = Signal to noise ratio (ranges,time).
  subroutine calc_ref_Z(Ndata,Fnn,TF,H,CC,ZdBz,Etha_NL,SNR,Dn,DSD,Vn_D)
    use variables, only : Nh, Nspec, PI, c_li, MRR_freq, Kappa2, sigma_B, Vnmodel
    implicit none
    integer, intent(in) :: Ndata, H(:), CC, Fnn(:,:,:)
    real, intent(in) :: TF(:)
    real, intent(out) :: ZdBz(:,:), Etha_NL(:,:), SNR(:,:), Vn_D(:,:)

    ! testing variables
    real, intent(out) :: Dn(:) ! new
    real, intent(out) :: DSD(:,:,:)  ! new DSD
    real :: del_u(Nh)
    real :: dVndD(Nh,Nspec)  ! -> Vn_D(N_Dh,Nspc)
    ! end of testing variables

    integer :: idx, idh, Sini, Send
    integer :: i_range(Nh), Fnoise(Nspec)
    real :: Etha_tot(Nh,Ndata), Etha(Nh, Nspec), F_corr(Nspec)
    real :: F_ave=-99, del_H, tmp(Nh), lambda
    real :: delDn(Nspec)

    DSD = 0. 
    !RR = 0. 
    select case (trim(Vnmodel))
       case('metek')
          ! The Dn are pre-calculated based on equation (2.5) METEK Physical Basis
          ! considering Vn = n*Delta_f*lambda/2= n*0.1887 [m/s]
          ! with n is the number of spectrum [0:63], 
          ! Following the relationship of terminal velocity and Dn for rain:
          ! with Vn_max ~ 9.64 m/s, therefore Dn_max ~ 9.95 mm
          ! Dn has Nspec elements:
          Dn = (/&
               0.1090,0.1416,0.1751,0.2094,0.2444,0.2801,0.3166,&
               0.3540,0.3922,0.4313,0.4713,0.5123,0.5544,0.5975,0.6418,&
               0.6873,0.7340,0.7821,0.8317,0.8828,0.9354,0.9898,1.0461,&
               1.1043,1.1646,1.2272,1.2922,1.3598,1.4303,1.5040,1.5810,&
               1.6618,1.7467,1.8361,1.9306,2.0308,2.1374,2.2513,2.3735,&
               2.5055,2.6488,2.8055,2.9786,3.1717,3.3902,3.6417,3.9380,&
               4.2987,4.7595,5.3987,6.4488,9.9506,9.9506,9.9506,9.9506,&
               9.9506,9.9506,9.9506,9.9506,9.9506,9.9506,9.9506,9.9506,&
               9.9506/) ! Drop diameters [mm]
       case DEFAULT
          stop "SORRY! No model to estimate particles size found!"
       end select

    
    delDn = (/(Dn(idx+1)-Dn(idx), idx=1,Nspec-1), 0.0/)  ! [mm]
    del_u = (/(1.+3.68E-5*H(idh)+1.71E-9*H(idh)**2,idh=1,Nh)/)  ! density altitude dependence factor
    Vn_D = matmul(reshape(del_u,(/Nh,1/)),&
         reshape(9.65-10.3*exp(-0.6*Dn),(/1,Nspec/))) ! Vn-Dn relationship

    dVndD = matmul(reshape(del_u,(/Nh,1/)),&
         reshape(6.18*exp(-0.6*Dn),(/1,Nspec/)))   !derivative of Vn resp. Dn
    Etha_tot = 0.
    SNR = 0.
    ZdBz = 0.
    lambda = 1E3*(c_li/MRR_freq)       ! MRR wavelenght [mm]
    del_H = real(H(2))-real(H(1))      ! MRR always keeps del_H constant
    i_range = (/(idh,idh=1,Nh)/)
    tmp = real(CC)*real(i_range)**2*del_H/(10.**20*TF)   ! [m]
    ! Calculating the Etha total reflectivity
    do idx=1,Ndata
       do idh=1,Nh
          if (idh>2.AND.idh<Nh-3) then
             Sini = 1
             Send = Nspec
          else
             Sini = 5
             Send = 60
          end if
          Fnoise = 0.
          Etha(idh,:) = 0.
          Fnoise(Sini:Send) = Fnn(idh,Sini:Send,idx)
          call Fnoise_corr(Fnoise,idh,Sini,Send,F_corr,F_ave)
          ! Integrating over the spectrum
          Etha_tot(idh,idx) = sum(F_corr(Sini:Send))
          SNR(idh,idx) = Etha_tot(idh,idx)/F_ave
          Etha_NL(idh,idx) = F_ave
          Etha(idh,Sini:Send) = tmp(idh)*F_corr(Sini:Send)
          ! Estimating Drop Size Distribution
          DSD(idh,:,idx) = -99.
          if(idh.LT.4) then
             DSD(idh,:,idx) = 1.E6*Etha(idh,:)*dVndD(idh,:)/sigma_B/0.1887
          else
             DSD(idh,:,idx) = 1.E6*Etha(idh,:)*dVndD(idh,:)*lambda**4/PI**5/Kappa2/0.1887/(Dn**6)
          end if
       end do ! end loop over range
       Etha_tot(:,idx) = tmp*Etha_tot(:,idx)
       ZdBz(:,idx) = 10.*log10(matmul(DSD(:,:,idx),delDn*Dn**6))
    end do    ! end loop over number of data
    ZdBz = -99.   ! for NaNs
    where (Etha_tot>0.0) ZdBz = 10.*log10(1.E6*Etha_tot*lambda**4/PI**5/Kappa2)
    where (SNR <= 0.) SNR = -99.
    print *, 'Reflectivity calculated'
    return
  end subroutine calc_ref_Z

  ! -----------------------------------------------------------------
  ! Subroutine Fnoise_corr
  ! Noise level estimation according the method specified on METEK's
  ! manual: Physical Basis.
  ! F_noise: Single raw spectral signal.
  ! range  : the range index
  ! AvNL   : Average Noise level, characteristic of every instrument.
  ! Sini,Send : initial and final spectral index.
  ! Fcorr  : Signal corrected by noise level.
  ! F_ave  : Estimated noise level for F_noise.
  subroutine Fnoise_corr(F_noise,range,Sini,Send,Fcorr,F_ave)
    use variables, only : Nspec, ANL
    implicit none
    integer, intent(in) :: range, Sini, Send
    integer, intent(inout) :: F_noise(:)
    real, intent(out) :: Fcorr(:), F_ave
    ! internal variables.
    real :: AvNL, Ave_last, tmp, F_max
    integer :: k, Ntotal, k1, k2
    integer, dimension(1) :: imax
    
    Fcorr = real(F_noise)
    AvNL = ANL(range)    ! Average Noise level 
    imax = 1             ! initial index of the spectral maximum
    ! computing the maximum value of the strectrum
    imax = maxloc(F_noise)
    F_max = real(maxval(F_noise))
    Ntotal = Send-Sini+1     ! Total elements accounted
    Ave_last = real(sum(F_noise))/real(Ntotal)
    ! Starting the noise level estimation
    do k = Sini,Send
       F_noise(imax) = 0d0       ! take out the spectral highest value
       Ntotal = 0
       F_ave = 0.
       Ntotal = count(F_noise(Sini:Send).NE.0)
       if (Ntotal.EQ.0) then
          F_ave = Ave_last
       else
          F_ave = sum(real(F_noise(Sini:Send)))/real(Ntotal)  ! average over all non-zero
       end if
       if (F_ave>=Ave_last.OR.F_ave<0.05*F_max) exit
       ! Finding out the nearest max value diff. to zero
       k1 = imax(1)
       k2 = imax(1)
       do while (F_noise(k1).EQ.0d0.AND.k1.GT.Sini)
          k1 = k1 - 1
       end do
       do while (F_noise(k2).EQ.0d0.AND.k2.LT.Send)
          k2 = k2 + 1
       end do
       if (F_noise(k1).GT.F_noise(k2)) then
          imax = k1
       else
          imax = k2
       end if
       Ave_last = F_ave
    end do
    ! Testing the quality of calculated Noise Level versus the Average 
    ! constant Noise Level extracted from different clear days.
    tmp = (F_ave-AvNL)/AvNL  !!!abs((F_ave-AvNL)/AvNL)
    if (tmp.lt.0.0.or.tmp.gt.0.2) F_ave = 1.25*AvNL
    !!!!where(F_noise.eq.0d0) F_noise = int(F_ave)
    !!!!Fcorr = Fcorr - real(F_noise)
    !!!!where(Fcorr.lt.1E-6) Fcorr = 0.0
    !!! if (AvNL.NE.0.0.AND.tmp.GT.0.2) F_ave = 1.2*AvNL
    ! checking the standard deviation of the left noise spectrum:
    !Ntotal = count(F_noise.NE.0)
    !tmp = sqrt(sum((real(F_noise)-F_ave)**2, MASK=F_noise.NE.0)/real(Ntotal))
    !where (F_noise.LT.(F_ave+2.*tmp).AND.F_noise.GT.0) Fcorr = 0.
    where (Fcorr>0.0) Fcorr = Fcorr - F_ave  ! correcting noise level
    where (Fcorr<0.0) Fcorr = 0.0  ! cleaning negative values
    Fcorr(1:Sini-1) = 0.0
    Fcorr(Send+1:Nspec) = 0.0
    return
  end subroutine Fnoise_corr

  ! ----------------------------------------------------------------------
  ! Soubroutine for estimation of the range of bottom and top rain layer.
  ! The output is the indexes for the range
  subroutine TopBottomRange(ZdBz,Ndata,idx1,idx2,FL)

    use variables, only : Nh
    implicit none
    real, intent(in) :: ZdBz(:,:)
    integer, intent(in) :: Ndata
    integer, intent(out) :: idx1(:), idx2(:), FL(:)
    ! local variables
    integer :: i,j,k
    integer, dimension(Nh,Ndata) :: Hmat
    integer, dimension(Ndata) :: tmp

    forall (i=1:Nh, j=1:NData) Hmat(i,j) = i

    ! Warning: here the threshould is <= 15 dBz:
    where (ZdBz.LT.16.OR.Hmat.GT.Nh-6) Hmat = -99

    ! Alising the bot/top range in order to avoid spikes.
    ! For the top range:
    tmp = maxval(Hmat,DIM=1,MASK=Hmat>0)
    where (tmp.LT.1.OR.tmp.GT.Nh) tmp = -99
    idx1 = tmp
    forall (k=2:Ndata-1,tmp(k)>0) idx1(k) = &
         sum(tmp(k-1:k+1),MASK=tmp(k-1:k+1).GT.0)/count(tmp(k-1:k+1).GT.0)

    ! For the bottom range:
    tmp = minval(Hmat,DIM=1,MASK=Hmat>0)
    where (tmp.LT.1.OR.tmp.GE.idx1) tmp = -99
    idx2 = tmp
    forall (k=2:Ndata-1,tmp(k)>0) idx2(k) = &
         sum(tmp(k-1:k+1),MASK=tmp(k-1:k+1).GT.0)/count(tmp(k-1:k+1).GT.0)

    FL = maxval(Hmat,DIM=1,MASK=ZdBz.GE.25.)+1
    where (FL.LE.1.OR.FL.GE.idx1) FL = -99

    return
  end subroutine TopBottomRange



  subroutine DSD2rainparam(Dn,NDn,Vn,diffDn,LWC,RR,Dm,Nw)
    ! =================================================================
    ! Subroutine to calculate rain parameters from a given DSD:
    ! * INPUT:
    ! Mandatory inputs are Drop diameters (Dn [mm]),
    ! Concentration (NDn [m^-3 mm^-1]),
    ! Falling velocitiy (Vn [m/s]).
    ! Optional input is the Drop diamters interval (delta_Dn [mm]) for integration
    ! * OUPUT (optionals):
    ! Liquid water content (LWC [g m^-3]), Rainrate (RR [mm/hr]),
    ! Mass-weight mean diamter (Dm [mm]) and
    ! DSD Scale paramter (Nw [m^-3 mm^-1]).
    ! * NUMBER OF ELEMENTS:
    ! Dn must be a vector containing the drop diameter spectrum
    ! NDn can be a matrix with number of obs. in rows and spectrum in columns
    ! Vn  same as NDn,
    ! The outputs are vectors with the same number of elements as NDn's rows.
    ! 
    ! Pablo Saavedra, (2010)
    ! pablosaa@uni-bonn.de
    ! ===================================================================

    use variables, only: RHOw, PI
    implicit none
    real, intent(in), dimension(:) :: Dn
    real, intent(in), dimension(:,:) :: NDn, Vn
    real, intent(in), optional, dimension(:) :: diffDn
    real, intent(out), optional, dimension(:) :: LWC, RR, Dm, Nw
    
    integer Ndata, Nd
    real, dimension(size(Dn)) :: dDn
    real, allocatable, dimension(:) ::  TMP
    ! constants
    !!real, parameter :: RHOw = 1.0  ! gr/cm^3
    !!real, parameter :: PI = 3.141659265359
  
    Nd = size(Dn)
    Ndata = size(NDn,1)

    if(.not.present(diffDn)) then   ! calculating delta Dn in case not given
       dDn(2:Nd) = Dn(2:Nd)-Dn(1:Nd-1)
       dDn(1) = dDn(2)
    else
       dDn = diffDn
    end if

    allocate(TMP(Ndata))
    TMP = matmul(NDn,Dn**3*dDn)
    ! Liquid water content:
    if(present(LWC)) LWC = 1E-3*RHOw*PI/6*TMP  ! [gr/m^3]
    ! Rain rate:
    if(present(RR)) RR = 6.0*PI*1E-4*matmul(NDn*Vn,Dn**3*dDn)   ! [mm/hr]
    ! Mass-weighted mean drop diameter:
    if(present(Dm)) Dm = matmul(NDn,Dn**4*dDn)/TMP  ! [mm]
    ! Scale parameter for drop concentration:
    if(present(Nw)) Nw = (256./6./Dm**4)*TMP          ! [m^-3 mm^-1]
  
    deallocate(TMP)
    return
  end subroutine DSD2rainparam

end module extras


!! OLD CODE
!!$  subroutine read_RAW_file(fname,idat,t_hr,H,Fnn,TF,CC)
!!$    implicit none
!!$    character(len=*), intent(in) :: fname
!!$    integer, intent(out) :: H(:), idat, CC, Fnn(:,:,:)
!!$    real, intent(out) :: t_hr(:), TF(:)
!!$    ! Internal variables
!!$    integer :: i, j, jspec, status !, len_data
!!$    integer :: yy, mm, dd, hh, mn, ss
!!$    logical :: NEWDATA
!!$    character(len=500) :: garbage
!!$
!!$    idat = 0
!!$    jspec = 0
!!$    NEWDATA = .FALSE.
!!$    ! Opening the MRR raw file
!!$    open(unit=18,file=fname,form='formatted',IOSTAT=status)
!!$    if (status /= 0) then
!!$       return
!!$    end if
!!$
!!$    do i=1,NData*(Nspc+3)
!!$       read(18,fmt='(A)',IOSTAT=status) garbage
!!$       if (status /= 0) then
!!$          exit
!!$       end if
!!$       if (garbage(1:3).EQ.'MRR'.AND.(.NOT.NEWDATA)) then
!!$
!!$          jspec = 0
!!$          idat = idat + 1
!!$          if (idat.GT.NData) then
!!$             write(6,*) 'Jump to end...'
!!$             exit
!!$          end if
!!$          read(garbage,FMT=20,ERR=15) yy, mm, dd, hh, mn, ss
!!$          j=index(garbage,'CC',.TRUE.)
!!$          read(garbage(j+2:j+10),'(I8)') CC
!!$                !CC = 2426000
!!$          t_hr(idat) = real(hh)+real(mn)/60.+real(ss)/3600.
!!$          NEWDATA = .TRUE.
!!$       else if (garbage(1:1).EQ.'H'.AND.NEWDATA) then
!!$          read(garbage,FMT=30,ERR=15) H
!!$       else if (garbage(1:2).EQ.'TF'.AND.NEWDATA) then
!!$          read(garbage,FMT=40,ERR=15) TF
!!$          NEWDATA = .FALSE.
!!$       else if (garbage(1:1).EQ.'F'.AND.(.NOT.NEWDATA)) then
!!$          ! Cleaning up the spectrum from rubish characters
!!$          forall (j=7:500,iachar(garbage(j:j))<48.OR.iachar(garbage(j:j))>57) &
!!$               garbage(j:j)=' '
!!$          ! print *, len_data(garbage),N_Dh
!!$          if (len_data(garbage).EQ.N_Dh.AND.idat.GT.0) then
!!$             read(garbage,'(1X,I2)') jspec
!!$             ! print *, jspec
!!$             jspec = jspec+1   ! in the file the spectrum starts from 0
!!$             read(garbage,FMT=50,IOSTAT=status,ERR=10) Fnn(:,jspec,idat)
!!$10           if (status>0) then
!!$                print *, 'load status > 0', status
!!$                print *, garbage
!!$                cycle
!!$             end if
!!$          else
!!$             print *, 'length not 32 or the file start with spectrum'
!!$          end if
!!$          garbage(:) = ' '
!!$       end if
!!$15     continue
!!$    end do
!!$
!!$    close(18)
!!$    write(6,*) 'End of File found...', idat
!!$
!!$    !     HERE STARTS THE DATA FORMATS
!!$20  FORMAT(4X,6I2) !,27X,I8)  ! Format for date + CC.
!!$30  FORMAT(3X,32I9)        ! Format for height.
!!$40  FORMAT(3X,32F9.6)      ! Format for TF and Z.
!!$50  FORMAT(3X,32I9)        ! Format for the spectrum
!!$    return
!!$  end subroutine read_RAW_file
!!$
!!$  integer function len_data(string)
!!$    implicit none
!!$    character(len=*), intent(in) :: string
!!$    integer :: bef
!!$    do bef=len(string),3,-1
!!$       if (string(bef:bef).NE.' ') then
!!$          len_data = (bef-3)/9
!!$          return
!!$       end if
!!$    end do
!!$    len_data = 0
!!$    return
!!$  end function len_data

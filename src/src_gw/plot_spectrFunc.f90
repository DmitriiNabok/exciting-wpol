
subroutine plot_spectrFunc()

  use modinput
  use modgw
  use mod_mpi_gw
  use m_getunit
  use mod_vxc,        only: vxcnn, read_vxcnn
  use mod_selfenergy, only: selfex, selfec, evalks, evalqp, eferqp, &
  &                         read_selfexnn, read_selfecnn, read_evalqp
  implicit none
  integer :: fid1, fid2, fid3, fid4, fid5, fid6, fid7
  integer :: n, iom, ik, ib
  character(30) :: frmt
  real(8),    allocatable :: enk(:), sf(:), sfunc(:,:,:), sftot(:)
  real(8),    allocatable :: ynk(:), sfd(:)
  complex(8), allocatable :: sx(:), sc(:), vxc(:)
  real(8)    :: scr, sci, om, om1, dE, eta, t1
  complex(8) :: sigma, dsigma

  input%gw%skipgnd = .true.
  call init_gw

if (rank == 0) then  

  if (allocated(vxcnn)) deallocate(vxcnn)
  allocate(vxcnn(ibgw:nbgw,kset%nkpt))
  vxcnn(:,:) = 0.d0
  call read_vxcnn()

  if (allocated(selfex)) deallocate(selfex)
  allocate(selfex(ibgw:nbgw,nkpt))
  selfex(:,:) = 0.d0
  call read_selfexnn()
  
  if (allocated(selfec)) deallocate(selfec)
  allocate(selfec(ibgw:nbgw,freq%nomeg,kset%nkpt))
  selfec(:,:,:) = 0.d0
  call read_selfecnn()

  n = nbgw-ibgw+1
  write(frmt,'("(",i8,"f14.6)")') 1+n
  call getunit(fid1)
  open(fid1,file='SF.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  call getunit(fid2)
  open(fid2,file='SE-Re.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  call getunit(fid3)
  open(fid3,file='SE-Im.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  call getunit(fid4)
  open(fid4,file='DF1.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  call getunit(fid5)
  open(fid5,file='DF2.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  call getunit(fid6)
  open(fid6,file='DF.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  call getunit(fid7)
  open(fid7,file='SF-D.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')

  allocate(sfunc(ibgw:nbgw,freq%nomeg,kset%nkpt))
  allocate(enk(ibgw:nbgw), sf(ibgw:nbgw))
  allocate(sx(ibgw:nbgw), sc(ibgw:nbgw), vxc(ibgw:nbgw))
  allocate(ynk(ibgw:nbgw), sfd(ibgw:nbgw))

  ! if (allocated(evalks)) deallocate(evalks)
  ! allocate(evalks(ibgw:nbgw,kset%nkpt))
  ! if (allocated(evalqp)) deallocate(evalqp)
  ! allocate(evalqp(ibgw:nbgw,kset%nkpt))
  ! call read_evalqp()

  call getSelfc(freq%nomeg, freq%freqs, selfec(nomax,:,ikvbm), evalsv(nomax,ikvbm), sigma, dsigma)
  dE = dble(selfex(nomax,ikvbm) + sigma - vxcnn(nomax,ikvbm))
  write(*,*) 'nomax, ikvbm, dE = ', nomax, ikvbm, dE

  do ik = 1, kset%nkpt

    write(fid1,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
    do iom = 1, freq%nomeg
      om  = freq%freqs(iom)
      enk = evalsv(ibgw:nbgw,ik)
      sx  = selfex(ibgw:nbgw,ik)
      sc  = selfec(ibgw:nbgw,iom,ik)
      ! do ib = ibgw, nbgw
      !   call getSelfc(freq%nomeg, freq%freqs, selfec(ib,:,ik), om-dE, sc(ib), dsigma)
      ! end do
      vxc = vxcnn(ibgw:nbgw,ik)
      sf  = 1.d0/pi*abs(aimag(sc)) / ( (om-enk-dble(sx+sc-vxc))**2 + aimag(sc)**2 )
      ! output
      om1 = om-efermi
      write(fid1,trim(frmt)) om1, sf
      write(fid2,trim(frmt)) om1, dble(sx+sc)
      write(fid3,trim(frmt)) om1, aimag(sx+sc)
      write(fid4,trim(frmt)) om1, om-enk
      write(fid5,trim(frmt)) om1, dble(sx+sc-vxc)
      write(fid6,trim(frmt)) om1, om-enk-dble(sx+sc-vxc)
      sfunc(:,iom,ik) = sf(:)

      !-------------------------------------
      ! Solution of the complex Dyson QPE
      !-------------------------------------
      do ib = ibgw, nbgw
        call getSelfc(freq%nomeg, freq%freqs, selfec(ib,:,ik), om, sigma, dsigma)
        ynk(ib) = -aimag(sigma) / (1.d0+dble(dsigma))
        t1 = om - enk(ib) - dble(sx(ib)) &
        &  - dble(sigma) + 2.d0*dble(sigma)*aimag(sigma)*ynk(ib) &
        &  + vxc(ib)
        t1 = t1*t1 + ynk(ib)*ynk(ib)
        sfd(ib) = 1.d0/pi*abs(ynk(ib)) / t1
      end do
      write(fid7,trim(frmt)) om1, sfd

    end do ! iom
    write(fid1,*); write(fid1,*)
    write(fid2,*); write(fid2,*)
    write(fid3,*); write(fid3,*)
    write(fid4,*); write(fid4,*)
    write(fid5,*); write(fid5,*)
    write(fid6,*); write(fid6,*)
    write(fid7,*); write(fid7,*)
    
  end do ! ik

  deallocate(enk,sf)
  deallocate(sx,sc,vxc)
  deallocate(ynk,sfd)

  close(fid1)
  close(fid2)
  close(fid3)
  close(fid4)
  close(fid5)
  close(fid6)
  close(fid7)

  !--------------------
  ! Gaussian smearing
  !--------------------
  ! if (.not.associated(input%gw%selfenergy%SpectralFunctionPlot)) &
  ! &  input%gw%selfenergy%SpectralFunctionPlot => getstructspectralfunctionplot(emptynode)
  ! eta = input%gw%selfenergy%SpectralFunctionPlot%eta
  ! write(*,*) 'eta = ', eta
  ! open(fid1,file='SF-G.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  ! do ik = 1, kset%nkpt
  !   do ib = ibgw, nbgw
  !     call apply_Gaussian_smearing(freq%nomeg, freq%freqs, sfunc(ib,:,ik), eta)
  !   end do
  !   do iom = 1, freq%nomeg
  !     write(fid1,trim(frmt)) freq%freqs(iom)-efermi, sfunc(:,iom,ik)
  !   end do
  !   write(fid1,*); write(fid1,*)
  ! end do
  ! close(fid1)


  !-------------------------
  ! Total spectral function
  !-------------------------
  allocate(sftot(freq%nomeg))
  sftot(:) = 0.d0
  do iom = 1, freq%nomeg
    do ik = 1, kset%nkpt
      sftot(iom) = sftot(iom) + kset%wkpt(ik)*sum(sfunc(ibgw:nbgw,iom,ik))
    end do
  end do
  open(fid1,file='SF-tot.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  do iom = 1, freq%nomeg
    write(fid1,'(2f18.8)') freq%freqs(iom)-efermi, sftot(iom)
  end do
  close(fid1)
  deallocate(sftot)

  deallocate(sfunc)

end if ! rank == 0

  return

contains

  subroutine apply_Gaussian_smearing(n,x,y,eta)
    implicit none
    integer, intent(in)    :: n
    real(8), intent(in)    :: x(n)
    real(8), intent(inout) :: y(n)
    real(8), intent(in)    :: eta
    integer :: i, j
    real(8) :: s
    real(8) :: z(n)

    ! eta = FWHM = 2 sqrt(2ln2) s
    s = eta / 2.35482

    do i = 1, n
      z(i) = 0.d0
      do j = 1, n
        z(i) = z(i) + y(j) * gaussian(x(i), x(j), s)
      end do
      z(i) = z(i) / dble(n)
    end do
    y = z

  end subroutine

  real(8) function gaussian(x,m,s)
    real(8), intent(in) :: x
    real(8), intent(in) :: m
    real(8), intent(in) :: s
    
    gaussian = exp(-(x-m)**2/(2.d0*s**2)) / sqrt(2.d0*pi*s**2)
    ! gaussian = exp(-(x-m)**2/(2.d0*s**2))
    
    return
  end function

    subroutine getSelfc(n, x, y, x0, f, df)
    integer,    intent(in)  :: n
    real(8),    intent(in)  :: x(n)
    complex(8), intent(in)  :: y(n)
    real(8),    intent(in)  :: x0
    complex(8), intent(out) :: f
    complex(8), intent(out) :: df
    ! local
    integer    :: i, i0, ii, j
    integer    :: np, np2
    real(8)    :: dx, t1, t2, xx
    complex(8) :: dy
    real(8),    allocatable :: c(:)
    complex(8), allocatable :: ya(:)
    real(8),    external    :: polynom

    ! only 1 point
    if (n <= 1) then
      write(*,*) 'ERROR(calc_evalqp_wpol::sigma): Not enough data points!'
      stop
    end if

    xx = max(x0, x(1))
    xx = min(x0, x(n))

if (.false.) then
    !-----------------------
    ! linear interpolation
    !-----------------------
    do i = 1, n
      if (x(i) >= xx) then
        if (i == 1) then
          df = (y(2)-y(1)) / (x(2)-x(1))
          f  = y(1) + df * (xx-x(1))
        else if (i == n) then
          df = (y(n)-y(n-1)) / (x(n)-x(n-1))
          f  = y(n-1) + df * (xx-x(n-1))
        else
          df = 0.5*(y(i+1)-y(i-1)) / (x(i)-x(i-1))
          f  = y(i-1) + df * (xx-x(i-1))
        end if
        exit
      end if
    end do
else
    !-----------------------
    ! polynomial fitting
    !-----------------------
    np = 4
    np2 = np/2
    allocate(ya(np), c(np))
    do i = 1, n
      if (x(i) >= xx) then
        if (i <= np2) then
          i0 = 1
        else if (i > n-np2) then
          i0 = n-np+1
        else
          i0 = i-np2
        end if
        do j = 1, np
          ii = i0+j-1
          ya(j) = y(ii)
        end do
        t1 = polynom( 0, np, x(i0), dble(ya), c, xx)
        t2 = polynom( 0, np, x(i0), aimag(ya), c, xx)
        f  = cmplx(t1, t2, 8)
        t1 = polynom( 1, np, x(i0), dble(ya), c, xx)
        t2 = polynom( 1, np, x(i0), aimag(ya), c, xx)
        df = cmplx(t1, t2, 8)
        exit
      end if
    end do
    deallocate(ya, c)
end if

    return
  end subroutine


end subroutine
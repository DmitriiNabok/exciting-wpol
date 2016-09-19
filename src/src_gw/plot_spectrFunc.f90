
subroutine plot_spectrFunc()

  use modinput
  use modgw
  use mod_mpi_gw
  use m_getunit
  use mod_vxc,        only: vxcnn, read_vxcnn
  use mod_selfenergy, only: selfex, selfec, read_selfexnn, read_selfecnn
  implicit none
  integer :: fid1, fid2, fid3, fid4, n, iom, ik, ie
  character(30) :: frmt
  real(8),    allocatable :: sfunc(:,:,:)
  real(8),    allocatable :: tvec(:)
  complex(8), allocatable :: zvec(:)
  complex(8) :: sigma
  real(8)    :: scr, sci, div, enk, vxc, om, eta

  input%gw%skipgnd = .true.
  call init_gw

  if (allocated(vxcnn)) deallocate(vxcnn)
  allocate(vxcnn(ibgw:nbgw,kset%nkpt))
  vxcnn(:,:) = 0.d0
  call read_vxcnn(kset,ibgw,nbgw,vxcnn)

  if (allocated(selfex)) deallocate(selfex)
  allocate(selfex(ibgw:nbgw,nkpt))
  selfex(:,:) = 0.d0
  call read_selfexnn(kset,ibgw,nbgw,selfex)
  
  if (allocated(selfec)) deallocate(selfec)
  allocate(selfec(ibgw:nbgw,freq%nomeg,kset%nkpt))
  selfec(:,:,:) = 0.d0
  call read_selfecnn(kset,freq,ibgw,nbgw,selfec)

  allocate(sfunc(ibgw:nbgw,freq%nomeg,kset%nkpt))
  sfunc(:,:,:) = 0.d0
  
  n = nbgw-ibgw+1
  write(frmt,'("(",i8,"f14.6)")') 1+n
  call getunit(fid1)
  open(fid1,file='SF.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  call getunit(fid2)
  open(fid2,file='SE-Re.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  call getunit(fid3)
  open(fid3,file='SE-Im.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  call getunit(fid4)
  open(fid4,file='DF.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')

  eta = input%gw%selfenergy%swidth

  allocate(tvec(ibgw:nbgw),zvec(ibgw:nbgw))
  do ik = 1, kset%nkpt
    write(fid1,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
    do iom = 1, freq%nomeg
      om  = freq%freqs(iom)
      do ie = ibgw, nbgw
        enk = evalsv(ie,ik)
        vxc = dble(vxcnn(ie,ik))
        zvec(ie) = selfex(ie,ik)+selfec(ie,iom,ik)
        scr = dble(zvec(ie))
        sci = imag(zvec(ie)) + eta
        tvec(ie) = om-enk-scr+vxc
        sfunc(ie,iom,ik) = 1.d0/pi*abs(sci)/(tvec(ie)**2+sci**2)
      end do
      write(fid1,trim(frmt)) om, sfunc(:,iom,ik)
      write(fid2,trim(frmt)) om, dble(zvec)
      write(fid3,trim(frmt)) om, imag(zvec)
      write(fid4,trim(frmt)) om, tvec
    end do
    write(fid1,*); write(fid1,*)
    write(fid2,*); write(fid2,*)
    write(fid3,*); write(fid3,*)
    write(fid4,*); write(fid4,*)
  end do
  deallocate(tvec,zvec)

  close(fid1)
  close(fid2)
  close(fid3)
  close(fid4)

  !--------------------
  ! Gaussian smearing
  !--------------------
  open(fid1,file='SF-G.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  do ik = 1, kset%nkpt
    do ie = ibgw, nbgw
      call apply_Gaussian_smearing(freq%nomeg,freq%freqs,sfunc(ie,:,ik))
    end do
    do iom = 1, freq%nomeg
      write(fid1,trim(frmt)) freq%freqs(iom), sfunc(:,iom,ik)
    end do
    write(fid1,*); write(fid1,*)
  end do
  close(fid1)

  return

contains

  subroutine apply_Gaussian_smearing(n,x,y)
    implicit none
    integer, intent(in)    :: n
    real(8), intent(in)    :: x(n)
    real(8), intent(inout) :: y(n)
    integer :: i, j
    real(8) :: s
    real(8) :: z(n)

    s = input%gw%selfenergy%SpectralFunctionPlot%eta

    z = 0.d0
    do i = 1, n
      ! z(i) = gaussian(x(i), 0.d0, s, 1.d0)
      if (y(i) > 1.d-6) then
        do j = 1, n
          z(j) = z(j) + gaussian(x(j), x(i), s, y(i))
        end do
      end if
    end do

    y = z

  end subroutine

  real(8) function gaussian(x,m,s,a)
    real(8), intent(in) :: x
    real(8), intent(in) :: m
    real(8), intent(in) :: s
    real(8), intent(in) :: a
    ! gaussian = a/(2.d0*pi*s**2)**0.5*exp(-(x-m)**2/(2.d0*s**2))
    gaussian = a*exp(-(x-m)**2/(2.d0*s**2))
    return
  end function


end subroutine

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
  real(8),    allocatable :: enk(:), sf(:)
  complex(8), allocatable :: sx(:), sc(:), vxc(:)
  real(8)    :: scr, sci, om

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
  open(fid4,file='DF.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')

  allocate(enk(ibgw:nbgw), sf(ibgw:nbgw))
  allocate(sx(ibgw:nbgw), sc(ibgw:nbgw), vxc(ibgw:nbgw))

  do ik = 1, kset%nkpt
    write(fid1,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
    do iom = 1, freq%nomeg
      om  = freq%freqs(iom)
      enk = evalsv(ibgw:nbgw,ik)
      sx  = selfex(ibgw:nbgw,ik)
      sc  = selfec(ibgw:nbgw,iom,ik)
      vxc = vxcnn(ibgw:nbgw,ik)
      sf  = 1.d0/pi*abs(aimag(sc)) / ( (om-enk-dble(sx+sc)+vxc)**2 + aimag(sc)**2 )
      ! output
      write(fid1,trim(frmt)) om, sf
      write(fid2,trim(frmt)) om, dble(sx+sc)
      write(fid3,trim(frmt)) om, aimag(sx+sc)
      write(fid4,trim(frmt)) om, om - ( enk + dble(sx+sc-vxc) )
    end do
    write(fid1,*); write(fid1,*)
    write(fid2,*); write(fid2,*)
    write(fid3,*); write(fid3,*)
    write(fid4,*); write(fid4,*)
  end do

  deallocate(enk,sf)
  deallocate(sx,sc,vxc)

  close(fid1)
  close(fid2)
  close(fid3)
  close(fid4)

  !--------------------
  ! Gaussian smearing
  !--------------------
  ! open(fid1,file='SF-G.DAT',form='FORMATTED',status='UNKNOWN',action='WRITE')
  ! do ik = 1, kset%nkpt
  !   do ie = ibgw, nbgw
  !     call apply_Gaussian_smearing(freq%nomeg,freq%freqs,sfunc(ie,:,ik))
  !   end do
  !   do iom = 1, freq%nomeg
  !     write(fid1,trim(frmt)) freq%freqs(iom)-efermi, sfunc(:,iom,ik)
  !   end do
  !   write(fid1,*); write(fid1,*)
  ! end do
  ! close(fid1)

end if ! rank == 0

  return

contains

  subroutine apply_Gaussian_smearing(n,x,y)
    implicit none
    integer, intent(in)    :: n
    real(8), intent(in)    :: x(n)
    real(8), intent(inout) :: y(n)
    integer :: i, j
    real(8) :: s, d
    real(8) :: z(n)

    ! eta = FWHM = 2 sqrt(2ln2) s
    s = input%gw%selfenergy%SpectralFunctionPlot%eta / 2.35482

    z = 0.d0
    do i = 1, n
      d = x(j) - x(i)
      if ( (y(i) > 1.d-6).and.(d < 100.d0) ) then
        do j = 1, n
          z(j) = z(j) + y(i)*gaussian(x(j), x(i), s)
        end do
      end if
    end do

    y = z

  end subroutine

  real(8) function gaussian(x,m,s)
    real(8), intent(in) :: x
    real(8), intent(in) :: m
    real(8), intent(in) :: s
    ! gaussian = exp(-(x-m)**2/(2.d0*s**2)) / sqrt(2.d0*pi*s**2)
    gaussian = exp(-(x-m)**2/(2.d0*s**2))
    return
  end function


end subroutine
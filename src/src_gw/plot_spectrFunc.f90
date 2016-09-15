
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
  real(8)    :: scr, sci, div, enk, vxc, om

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

  allocate(tvec(ibgw:nbgw),zvec(ibgw:nbgw))
  do ik = 1, kset%nkpt
    write(fid1,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
    do iom = 1, freq%nomeg
      om  = freq%freqs(iom)-efermi
      do ie = ibgw, nbgw
        enk = evalsv(ie,ik)-efermi
        vxc = dble(vxcnn(ie,ik))
        zvec(ie) = selfex(ie,ik)+selfec(ie,iom,ik)
        scr = dble(zvec(ie))
        sci = imag(zvec(ie))
        tvec(ie) = om-enk-scr+vxc
        ! sfunc(ie,iom,ik) = 1.d0/pi*abs(sci)/(tvec(ie)**2+sci**2)
        sfunc(ie,iom,ik) = abs( 1.d0/pi * 1.d0/(om-enk-zvec(ie)+vxc) )
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

  return
end subroutine
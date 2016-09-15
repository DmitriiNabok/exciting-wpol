
module mod_vxc

  ! APW-APW exchange-correlation  integrals
  real(8), allocatable :: vxcraa(:,:,:,:,:,:)
      
  ! local-orbital-APW exchange-correlation  integrals
  real(8), allocatable :: vxcrloa(:,:,:,:,:)
      
  ! local-orbital-local-orbital exchange-correlation  integrals
  real(8), allocatable :: vxcrlolo(:,:,:,:)

  ! G-space interstitial exchange-correlation potential
  complex(8), allocatable :: vxcig(:)
      
  ! diagonal matrix elements of the exchange-correlation potential
  complex(8), allocatable :: vxcnn(:,:)

contains

!------------------------------------------------------------------
  subroutine write_vxcnn(kset,ibgw,nbgw,zmat)
    use m_getunit
    use mod_kpointset
    implicit none
    type(k_set), intent(in) :: kset
    integer,     intent(in) :: ibgw, nbgw
    complex(8),  intent(in) :: zmat(ibgw:nbgw,kset%nkpt)
    integer :: fid, ik, ie

    call getunit(fid)

    ! text format
    open(fid,file='VXCNN.DAT',form='FORMATTED',status='UNKNOWN')
    do ik = 1, kset%nkpt
      write(fid,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
      do ie = ibgw, nbgw
        write(fid,'(i6,f18.6)') ie, dble(vxcnn(ie,ik))
      end do
      write(fid,*); write(fid,*)
    end do
    close(fid)

    ! fortran binary format
    open(fid,file='VXCNN.OUT',form='UNFORMATTED',status='UNKNOWN')
    write(fid) ibgw, nbgw, kset%nkpt, vxcnn
    close(fid)

    return
  end subroutine

!------------------------------------------------------------------
  subroutine read_vxcnn(kset,ibgw,nbgw,zmat)
    use m_getunit
    use mod_kpointset
    implicit none
    type(k_set), intent(in) :: kset
    integer,     intent(in) :: ibgw, nbgw
    complex(8),  intent(in) :: zmat(ibgw:nbgw,kset%nkpt)
    integer :: ib, nb, nk
    integer :: fid

    call getunit(fid)
    open(fid,file='VXCNN.OUT',form='UNFORMATTED',status='UNKNOWN')
    read(fid) ib, nb, nk
    close(fid)

    if (nk.ne.kset%nkpt) then
      write(*,*)'ERROR(mod_vxc::read_vxcnn): Wrong number of k-points'
      write(*,*)'    nk=', nk, '    nkpt=', kset%nkpt
      stop
    end if

    if ((ib.ne.ibgw).or.(nb.ne.nbgw)) then
      write(*,*)'ERROR(mod_vxc::read_vxcnn): Different number of bands'
      write(*,*)'    ib=',   ib, '    nb=', nb
      write(*,*)'  ibgw=', ibgw, '  nbgw=', nbgw
      stop
    end if
    
    open(fid,file='VXCNN.OUT',form='UNFORMATTED',status='UNKNOWN')
    read(fid) ib, nb, nk, vxcnn
    close(fid)
      
    return
  end subroutine
      
end module

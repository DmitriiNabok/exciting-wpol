
module mod_selfenergy

    ! Analytical continuation type
    integer :: iopac ! 0/1

    !--------------!
    ! self-energy  ! 
    !--------------!
    
    ! The exchange self-energy
    complex(8), allocatable :: selfex(:,:)
    
    ! Sum_ij{M^i*W^c_{ij}*conjg(M^j)}
    complex(8), allocatable :: mwm(:,:,:)
    target mwm
    
    ! The correlation sel-fenergy
    complex(8), allocatable :: selfec(:,:,:)
    complex(8), allocatable :: selfecSR(:,:,:)
    complex(8), allocatable :: selfecw2(:,:,:)
    
    ! Correction factors for (q^-1) and (q^-2) singularities
    real(8) :: singc1
    real(8) :: singc2
      
    !-------------!
    ! QP Energy   !
    !-------------!
    
    ! Original KS energies (evalsv will updated via self-consistent cycle)
    real(8), allocatable :: evalks(:,:)
    
    ! QP energies
    real(8) :: eferqp
    real(8), allocatable :: evalqp(:,:)
    
    ! Linearization (renormalization) factor
    real(8),    allocatable :: znorm(:,:)
    
    ! AC to the real axis of the correlation self-energy (selfec)
    complex(8), allocatable :: sigc(:,:)
    
    ! COHSEX approximation
    complex(8), allocatable :: sigsx(:,:) ! Screened exchange
    complex(8), allocatable :: sigch(:,:) ! Coulomb hole
    
contains

    !---------------------------------------------------------------------------
    subroutine init_selfenergy(ibgw,nbgw,nkpt,nomeg)
        use modinput
        implicit none
        integer, intent(in) :: ibgw, nbgw
        integer, intent(in) :: nkpt
        integer, intent(in) :: nomeg
       
        ! KS eigenvalues
        if (allocated(evalks)) deallocate(evalks)
        allocate(evalks(ibgw:nbgw,nkpt))
        evalks(:,:) = 0.d0
        
        ! Quasi-Particle energy
        if (allocated(evalqp)) deallocate(evalqp)
        allocate(evalqp(ibgw:nbgw,nkpt))
        evalqp(:,:) = 0.d0
        
        if (allocated(selfex)) deallocate(selfex)
        allocate(selfex(ibgw:nbgw,nkpt))
        selfex(:,:) = 0.d0
        
        if (input%gw%taskname.ne.'g0w0_x') then
          if (allocated(selfec)) deallocate(selfec)
          allocate(selfec(ibgw:nbgw,nomeg,nkpt))
          selfec(:,:,:) = 0.d0
          if (input%gw%selfenergy%secordw) then
            if (allocated(selfecSR)) deallocate(selfecSR)
            allocate(selfecSR(ibgw:nbgw,nomeg,nkpt))
            selfecSR(:,:,:) = 0.d0
            if (allocated(selfecw2)) deallocate(selfecw2)
            allocate(selfecw2(ibgw:nbgw,nomeg,nkpt))
            selfecw2(:,:,:) = 0.d0
          end if
          if (input%gw%taskname.ne.'cohsex') then
            ! Correlation self-energy at real frequencies after AC procedure
            if (allocated(sigc)) deallocate(sigc)
            allocate(sigc(ibgw:nbgw,nkpt))
            sigc(:,:) = 0.d0
            ! Renormalization (linearization) factors
            if (allocated(znorm)) deallocate(znorm)
            allocate(znorm(ibgw:nbgw,nkpt))
            znorm(:,:) = 0.d0
          else
            ! COHSEX approximation
            if (allocated(sigsx)) deallocate(sigsx)
            allocate(sigsx(ibgw:nbgw,nkpt))
            sigsx(:,:) = 0.d0
            if (allocated(sigch)) deallocate(sigch)
            allocate(sigch(ibgw:nbgw,nkpt))
            sigch(:,:) = 0.d0
          end if ! cohsex
        end if
        
    end subroutine
        
    !---------------------------------------------------------------------------
    subroutine delete_selfenergy()
      if (allocated(evalks)) deallocate(evalks)
      if (allocated(evalqp)) deallocate(evalqp)
      if (allocated(selfex)) deallocate(selfex)
      if (allocated(selfec)) deallocate(selfec)
      if (allocated(selfecSR)) deallocate(selfecSR)
      if (allocated(selfecw2)) deallocate(selfecw2)
      if (allocated(znorm))  deallocate(znorm)
      if (allocated(sigc))   deallocate(sigc)
      if (allocated(sigsx))  deallocate(sigsx)
      if (allocated(sigch))  deallocate(sigch)
    end subroutine
        
    !---------------------------------------------------------------------------
    subroutine write_selfenergy(ibgw,nbgw,nkpt,nomeg)
      use modinput
      implicit none
      integer, intent(in) :: ibgw, nbgw
      integer, intent(in) :: nkpt
      integer, intent(in) :: nomeg
      ! local variables
      integer :: fid, ie, ik, iom
      fid = 777
      ! exchange
      open(fid,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
      write(fid) ibgw, nbgw, nkpt, selfex
      close(fid)
      ! correlation
      if (input%gw%taskname.ne.'g0w0_x') then
        open(fid,file='SELFC.OUT',form='UNFORMATTED',status='UNKNOWN')
        write(fid) ibgw, nbgw, nomeg, nkpt, selfec
        close(fid)
        if (input%gw%taskname=='cohsex') then
          open(fid,file='COHSEX.OUT',form='UNFORMATTED',status='UNKNOWN')
          write(fid) ibgw, nbgw, nkpt, selfec, sigsx, sigch
          close(fid)
        end if
      end if
      !-----------------
      ! Text version
      !-----------------
      open(fid,file='SELFX.DAT',form='FORMATTED',status='UNKNOWN')
      write(fid,*) '# ik    ie    selfex'
      do ik = 1, nkpt
      do ie = ibgw, nbgw
        write(fid,'(2i6,2f18.6)') ik, ie, selfex(ie,ik)
      end do
      write(fid,*)
      end do
      close(fid)
      if (input%gw%taskname.ne.'g0w0_x') then
        open(fid,file='SELFC.DAT',form='FORMATTED',status='UNKNOWN')
        write(fid,*) '# iom    ik    ie    selfec'
        do iom = 1, nomeg
        do ik = 1, nkpt
        do ie = ibgw, nbgw
          write(fid,'(3i6,2f18.6)') iom, ik, ie, selfec(ie,iom,ik)
        end do
        write(fid,*)
        end do
        end do
        close(fid)
        if (input%gw%selfenergy%secordw) then
          !open(fid,file='SELFCW2.DAT',form='FORMATTED',status='UNKNOWN')
          open(fid,file='SELFCSR.DAT',form='FORMATTED',status='UNKNOWN')
          write(fid,*) '# iom    ik    ie    selfec'
          do iom = 1, nomeg
          do ik = 1, nkpt
          do ie = ibgw, nbgw
            write(fid,'(3i6,2f18.6)') iom, ik, ie, selfecSR(ie,iom,ik)
          end do
          write(fid,*)
          end do
          end do
          close(fid)
        end if
      end if
    end subroutine

!--------------------------------------------------------------------
! exchange selfenergy
!--------------------------------------------------------------------

    subroutine write_selfexnn()
      use modgw, only : kset, ibgw, nbgw
      use m_getunit
      implicit none
      integer :: fid, ik, ie

      call getunit(fid)
      ! text format
      open(fid,file='SELFX.DAT',form='FORMATTED',status='UNKNOWN')
      do ik = 1, kset%nkpt
        write(fid,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
        do ie = ibgw, nbgw
          write(fid,'(i6,f18.6)') ie, dble(selfex(ie,ik))
        end do
        write(fid,*); write(fid,*)
      end do
      close(fid)
      ! fortran binary format
      open(fid,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
      write(fid) ibgw, nbgw, kset%nkpt, selfex(ibgw:nbgw,1:kset%nkpt)
      close(fid)
      return
    end subroutine

    subroutine read_selfexnn()
      use modgw, only : kset, ibgw, nbgw
      use m_getunit
      implicit none
      integer :: ib, nb, nk
      integer :: fid

      call getunit(fid)
      open(fid,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
      read(fid) ib, nb, nk
      close(fid)

      if (nk /= kset%nkpt) then
        write(*,*)'ERROR(mod_selfenergy::read_selfexnn): Wrong number of k-points'
        write(*,*)'    nk=', nk, '    nkpt=', kset%nkpt
        stop
      end if

      if ((ib /= ibgw).or.(nb /= nbgw)) then
        write(*,*)'ERROR(mod_selfenergy::read_selfexnn): Different number of bands'
        write(*,*)'    ib=',   ib, '    nb=', nb
        write(*,*)'  ibgw=', ibgw, '  nbgw=', nbgw
        stop
      end if
    
      open(fid,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
      read(fid) ib, nb, nk, selfex(ibgw:nbgw,1:kset%nkpt)
      close(fid)

      return
    end subroutine

!--------------------------------------------------------------------
! correlation selfenergy
!--------------------------------------------------------------------

    subroutine write_selfecnn()
      use modgw, only : kset, freq, ibgw, nbgw
      use m_getunit
      implicit none
      integer :: iom, n, ik
      integer :: fid1, fid2
      character(80) :: frmt

      ! text format
      n = nbgw-ibgw+1
      write(frmt,'("(",i8,"f14.6)")') 1+n
      ! write(*,*) trim(frmt)

      call getunit(fid1)
      open(fid1, File='SELFC-WPOL-Re.DAT', Action='WRITE')
      call getunit(fid2)
      open(fid2, File='SELFC-WPOL-Im.DAT', Action='WRITE')
      do ik = 1, kset%nkpt
        write(fid1,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
        write(fid2,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
        do iom = 1, freq%nomeg
          write(fid1,trim(frmt)) freq%freqs(iom), dble(selfec(:,iom,ik))
          write(fid2,trim(frmt)) freq%freqs(iom), imag(selfec(:,iom,ik))
        end do
        write(fid1,*); write(fid1,*)
        write(fid2,*); write(fid2,*)
      end do
      close(fid1)
      close(fid2)

      ! binary format
      open(fid1,file='SELFC.OUT',form='UNFORMATTED',status='UNKNOWN')
      write(fid1) ibgw, nbgw, freq%nomeg, kset%nkpt, selfec(ibgw:nbgw,1:freq%nomeg,1:kset%nkpt)
      close(fid1)

      return
    end subroutine
    
    subroutine read_selfecnn()
      use modgw, only : kset, freq, ibgw, nbgw
      use m_getunit
      implicit none
      integer :: ib, nb, nk, no
      integer :: fid

      call getunit(fid)

      open(fid,file='SELFC.OUT',form='UNFORMATTED',status='UNKNOWN')
      read(fid) ib, nb, no, nk
      close(fid)

      if (nk /= kset%nkpt) then
        write(*,*)'ERROR(mod_selfenergy::read_selfecnn)): Wrong number of k-points'
        write(*,*)'    nk=', nk, '    nkpt=', kset%nkpt
        stop
      end if

      if ((ib /= ibgw).or.(nb /= nbgw)) then
        write(*,*)'WARNING(mod_selfenergy::read_selfecnn)): Different number of bands'
        write(*,*)'    ib=',   ib, '    nb=', nb
        write(*,*)'  ibgw=', ibgw, '  nbgw=', nbgw
        stop
      end if

      if (no /= freq%nomeg) then
        write(*,*)'ERROR(mod_selfenergy::read_selfecnn)): Wrong number of frequencies'
        write(*,*)'    no=', no, '    freq%nomeg=', freq%nomeg
        stop
      end if
    
      open(fid,file='SELFC.OUT',form='UNFORMATTED',status='UNKNOWN')
      read(fid) ib, nb, no, nk, selfec(ibgw:nbgw,1:freq%nomeg,1:kset%nkpt)
      close(fid)
     
      return
    end subroutine
end module

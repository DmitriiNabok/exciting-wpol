
module mod_wpol

  use modinput
  use modgw
  use mod_mpi_gw

  implicit none
  
  integer, private :: ndim, mdim
  integer, private :: mbdim
  integer, public  :: nvck

  ! index mapping array
  integer,    allocatable :: a2vck(:,:)
  integer,    allocatable :: vck2a(:,:,:)
  private :: a2vck, vck2a

  ! data arrays
  real(8),    allocatable :: tvck(:)
  complex(8), allocatable :: wvck(:,:)
  public :: tvck, wvck

  complex(8), allocatable :: md(:,:)
  complex(8), allocatable :: dmmd(:,:)
  complex(8), allocatable :: wij(:,:,:)
  private :: md, dmmd, wij

contains

!--------------------------------------------------------------------------------
  subroutine test_wpol()
    implicit none
    integer :: iq

    !=================
    ! Initialization
    !=================
    call init_gw
    call clean_gndstate

#ifdef MPI
    call set_mpi_group(kqset%nkpt)
    call mpi_set_range(nproc_row, &
    &                  myrank_row, &
    &                  kset%nkpt, 1, &
    &                  iqstart, iqend)
#else
    iqstart = 1
    iqend = kset%nkpt
#endif

    !===========================
    ! Momentum matrix elements
    !===========================
    if (.not.input%gw%rpmat) then
      call calcpmatgw
    end if
    
    ! Setup working array dimensions and index mappings
    call set_wpol_indices()

    !=========================
    ! Main loop over q-points
    !=========================

    if (myrank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if

    ! each process does a subset
    do iq = iqstart, iqend

      write(*,*)
      write(*,*) '(mod_wpol::test_wpol) q-point cycle, iq = ', iq
    
      Gamma = gammapoint(kqset%vqc(:,iq))
          
      !================================================
      ! Calculate interstitial product basis functions
      !================================================
      matsiz = locmatsiz+Gqset%ngk(1,iq)
      call diagsgi(iq)
      call calcmpwipw(iq)
    
      !======================================
      ! Calculate the bare Coulomb potential
      !======================================
      call calcbarcmb(iq)
      ! set v-diagonal MB
      call setbarcev(input%gw%barecoul%barcevtol)


      ! step 1: calculate M*D^{1/2} and D^2 + D^{1/2}*M^{+}*v*M*D^{1/2} matrices
      call calc_md_dmmd(iq)

      ! step 2: diagonalize DMMD
      call diagonalize_dmmd()

      ! step 3: calculate w_{vck} Eq. (2.20)
      call calc_wvck(iq)

      ! step 4: calculate W_{ij} Eq. (2.22)
      call calc_wmat()

      ! output results for visualization
      call print_wmat(iq,mbsiz,freq%nomeg,wij)

      call delete_coulomb_potential
      call clear_wpol()
      deallocate(tvck)
      deallocate(wvck)

    end do ! q-points

    ! delete index mapping arrays
    call del_wpol_indices()

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine set_wpol_indices()
    implicit none
    integer :: i, ik, n, m

    ! total number of states including the core ones
    if (input%gw%coreflag=='all') then
      ndim = nomax+ncg
    else
      ndim = nomax
    end if
    mdim  = nstdf-numin+1

    nvck = kqset%nkpt*nomax*(nstdf-numin+1)

    ! map a -> {vck}
    if (allocated(a2vck)) deallocate(a2vck)
    allocate(a2vck(3,nvck))

    ! map {vck} -> a
    if (allocated(vck2a)) deallocate(vck2a)
    allocate(vck2a(1:ndim,numin:nstdf,kqset%nkpt))

    i = 0
    do ik = 1, kqset%nkpt
      do m = numin, nstdf
      do n = 1, ndim
        i = i+1
        a2vck(1,i) = n
        a2vck(2,i) = m
        a2vck(3,i) = ik
        vck2a(n,m,ik) = i
      end do
      end do
    end do

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine del_wpol_indices()
    if (allocated(a2vck)) deallocate(a2vck)
    if (allocated(vck2a)) deallocate(vck2a)
  end subroutine

!--------------------------------------------------------------------------------
  subroutine calc_md_dmmd(iq)
    implicit none
    integer, intent(in) :: iq

    integer    :: ispn, ik, jk, ikp, jkp
    integer    :: ic, icg, ia, is, ias
    integer    :: i, j, n, m
    integer(8) :: recl
    real(8)    :: de
    real(8)    :: q0eps(3), modq0
    complex(8) :: zt1

    real(8),    allocatable :: d(:)
    complex(8), allocatable :: evecsv(:,:,:)

    mbdim = mbsiz
    if (Gamma) then
      mbdim = mbsiz+1
      ! q->0 direction
      q0eps(:) = input%gw%scrcoul%q0eps(:)
      modq0    = sqrt(q0eps(1)**2+q0eps(2)**2+q0eps(3)**2)
      q0eps(:) = q0eps(:)/modq0
      ! val-val
      if (allocated(pmatvv)) deallocate(pmatvv)
      allocate(pmatvv(nomax,numin:nstdf,3))
      inquire(iolength=recl) pmatvv
      open(fid_pmatvv,File=fname_pmatvv, &
      &    Action='READ',Form='UNFORMATTED',&
      &    Access='DIRECT',Status='OLD',Recl=recl)
      ! core-val 
      if (input%gw%coreflag=='all') then
        if (allocated(pmatcv)) deallocate(pmatcv)
        allocate(pmatcv(ncg,numin:nstdf,3))
        inquire(iolength=recl) pmatcv
        open(fid_pmatcv,File=fname_pmatcv, &
        &    Action='READ',Form='UNFORMATTED', &
        &    Access='DIRECT',Status='OLD',Recl=recl)
      end if
    end if

    ! global arrays
    allocate(md(mbdim,nvck))
    md(:,:) = 0.d0  
    allocate(dmmd(nvck,nvck))
    dmmd(:,:) = 0.d0

    ! local data
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))    
    allocate(evecsv(nmatmax,nstsv,nspinor))
    allocate(minmmat(mbsiz,ndim,numin:nstdf))
    allocate(d(nvck))

    ! Loop over k-points
    do ispn = 1, nspinor
    do ik = 1, kqset%nkpt
        
      write(*,*)
      write(*,*) '(mod_wpol::calc_md_dmmd): rank, (iq, ik):', myrank, iq, ik
    
      ! k-q point
      jk = kqset%kqid(ik,iq)

      ! irreducible k-point index
      ikp = kset%ik2ikp(ik)
      jkp = kset%ik2ikp(jk)
      
      ! get KS eigenvectors
      call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
      eveckp = conjg(evecsv(:,:,ispn))
      call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
      eveck = evecsv(:,:,ispn)
        
      ! Calculate M^i_{nm}+M^i_{cm}
      call expand_evec(ik,'t')
      call expand_evec(jk,'c')
      call expand_products(ik,iq,1,ndim,nomax,numin,nstdf,-1,minmmat)

      ! read the momentum matrix elements
      if (Gamma) call getpmatkgw(ik)
        
      ! Setup \tilde{M}_{ab}(q) * 2D^{1/2}
      do m = numin, nstdf
      do n = 1, ndim
        i = vck2a(n,m,ik)
        if (n <= nomax) then
          de = evalsv(m,jkp)-evalsv(n,ikp)
          if (Gamma) then
            ! p \cdot q
            zt1 =  pmatvv(n,m,1)*q0eps(1)+ &
            &      pmatvv(n,m,2)*q0eps(2)+ &
            &      pmatvv(n,m,3)*q0eps(3)
          end if
        else
          icg = n-nomax
          is  = corind(icg,1)
          ia  = corind(icg,2)
          ic  = corind(icg,3)
          ias = idxas(ia,is)
          de  = evalsv(m,jkp)-evalcr(ic,ias)
          if (Gamma) then
            ! p \cdot q
            zt1 =  pmatcv(icg,m,1)*q0eps(1)+ &
            &      pmatcv(icg,m,2)*q0eps(2)+ &
            &      pmatcv(icg,m,3)*q0eps(3)
          end if
        end if
        if (abs(de) < 1.d-6) then
          write(*,*) "(mod_wpol::calc_md_dmmd) Problem with eigenvalues! ", n, m
          stop
        end if

        d(i) = de
        md(1:mbsiz,i) = sqrt(de)*minmmat(1:mbsiz,n,m)

        ! singular term
        if (Gamma) md(mbsiz+1,i) = -sqrt(4.d0*pi/omega)*zt1/sqrt(de)

      end do
      end do

    end do ! ik
    end do ! nspinor

    deallocate(minmmat)
    deallocate(eveckalm)
    deallocate(eveckpalm)
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(evecsv)
    if (Gamma) then
      close(fid_pmatvv)
      deallocate(pmatvv)
      if (input%gw%coreflag=='all') then
        deallocate(pmatcv)
        close(fid_pmatcv)
      end if
    end if

    ! D^{1/2} M^{+} v M D^{1/2} = D^{1/2} \tilde{M}^{+} \tilde{M} D^{1/2}

    ! regular term
    call zgemm( 'c', 'n', nvck, nvck, mbdim, &
    &           zone, md, mbdim, md, mbdim, &
    &           zzero, dmmd, nvck)
    
    ! account for prefactor 4/Nk
    dmmd(:,:) = dmmd(:,:)*4.d0/dble(kqset%nkpt)

    ! D^2 + D^{1/2} M^{+} v M D^{1/2}
    do i = 1, nvck
      dmmd(i,i) = d(i)*d(i) + dmmd(i,i)
    end do
    deallocate(d)

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine diagonalize_dmmd()
    implicit none
    integer                 :: i, info, lwork
    ! RWORK dimension should be at least MAX(1,3*N-2)
    real(8),    allocatable :: rwork(:)
    complex(8), allocatable :: work(:)

    ! matrix eigenvalues
    allocate(tvck(nvck))

    write(*,*)
    write(*,*) 'Diagonalization with zheev'
    write(*,*)
    lwork = 2*nvck
    allocate(work(lwork),rwork(3*nvck))
    call zheev( 'v', 'u', nvck, dmmd, nvck, &
    &           tvck, work, lwork, rwork, info )
    deallocate(work,rwork)

    ! check for convergence
    if ( info > 0 ) then
      write(*,*)'Error(task_wpol): the algorithm failed to compute eigenvalues'
      write(*,*)'info = ', info
      stop
    end if

    ! following definition of Eq.(2.20) tau = sqrt(lambda)
    do i = 1, nvck
      if (abs(tvck(i)) > 1.d-8) then
        tvck(i) = sqrt(tvck(i))
      else
        write(*,*) "Warning(mod_wpol::diagonalize_dmmd) Negative or zero eigenvalues!"
        write(*,*) "i, tvck: ", i, tvck(i)
      end if
    end do

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine calc_wvck(iq)
    implicit none
    integer, intent(in) :: iq
    integer :: im, i
    real(8) :: t1

    ! v * M * D^{1/2} = v^{1/2} * \tilde{M} * 2D^{1/2}
    do im = 1, mbsiz
      md(im,:) = sqrt(barcev(im))*md(im,:)
    end do
    ! singular term
    if (Gamma) md(mbsiz+1,:) = sqrt(4.d0*pi)*md(mbsiz+1,:)

    ! w_{vck} (2.20)
    allocate(wvck(mbdim,nvck))
    call zgemm( 'n', 'n', mbdim, nvck, nvck, &
    &           zone, md, mbdim, dmmd, nvck, &
    &           zzero, wvck, mbdim)
    
    ! think of storing w_{vck} to be reused for calculating \Sigma_c

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine calc_wmat()
    implicit none
    integer    :: iom, i
    complex(8) :: zt1, zif, zieta

    ! W_{ij}(q,\omega) (2.22) (regular part only)
    allocate(wij(mbsiz,mbsiz,freq%nomeg))

    zieta = zi*input%gw%scrcoul%swidth
    select case (input%gw%freqgrid%fconv)
      case('refreq')
        ! real axis
        zif = zone
      case('imfreq')
        ! imaginary axis
        zif = zi
      case default
        write(*,*) "Accepted options: refreq or imfreq"
        stop
    end select

    do iom = 1, freq%nomeg

      do i = 1, nvck
        ! skip "bad" eigenvalues if they are
        if (abs(tvck(i)) > 1.d-8) then
          zt1 =  0.5d0 / tvck(i) *                     &
          &    ( 1.d0/(zif*freq%freqs(iom)-tvck(i)+zieta) &
          &     -1.d0/(zif*freq%freqs(iom)+tvck(i)-zieta) )
          md(:,i) = zt1*wvck(:,i) ! reuse the array md
        else
          md(:,i) = zzero
        end if
      end do

      ! Important: correct dimensions for q=0 case
      call zgemm( 'n', 'c', mbsiz, mbsiz, nvck, &
      &           zone, md(1:mbsiz,:), mbsiz, wvck(1:mbsiz,:), mbsiz, &
      &           zzero, wij(:,:,iom), mbsiz)

      ! W = v + wij
      ! do im = 1, mbsiz
      !   wij(im,im,iom) = barcev(im)+wij(im,im,iom)
      ! end do

    end do ! iom

    ! account for prefactor 4/Nk
    wij(:,:,:) = wij(:,:,:)*4.d0/dble(kqset%nkpt)

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine clear_wpol()
    if (allocated(md))   deallocate(md)
    if (allocated(dmmd)) deallocate(dmmd)
    if (allocated(wij))  deallocate(wij)
  end subroutine

!--------------------------------------------------------------------------------
  subroutine print_wmat(iq,mbsiz,nomeg,zmat)
    use m_getunit
    implicit none
    integer,    intent(in) :: iq
    integer,    intent(in) :: mbsiz
    integer,    intent(in) :: nomeg
    complex(8), intent(in) :: zmat(mbsiz,mbsiz,nomeg)
    integer :: iom, im
    integer :: fid
    character(80) :: fname

    ! store q-dependent Wij
    call getunit(fid)

    write(fname,'("WMAT-mat-q",I4.4,".OUT")') iq
    open(fid, File=trim(fname), Action='WRITE')
    do iom = 1, nomeg, 100
      do im = 1, mbsiz
        write(fid,'(i8, 2f16.6)') im, zmat(im,im,iom)
      end do
      write(fid,*)
    end do
    close(fid)

    write(fname,'("WMAT-iom-q",I4.4,".OUT")') iq
    open(fid, File=trim(fname), Action='WRITE')
    do im = mbsiz-100, mbsiz, 10
      do iom = 1, nomeg
        write(fid,'(3f16.6)') freq%freqs(iom), zmat(im,im,iom)
      end do
      write(fid,*)
    end do
    close(fid)

    return
  end subroutine

end module

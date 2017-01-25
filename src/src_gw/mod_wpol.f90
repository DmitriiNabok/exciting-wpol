
module mod_wpol

  use modinput
  use modgw
  use mod_mpi_gw

  implicit none
  
  integer, private :: ndim, mdim
  integer, public  :: nvck, nvck0, mbdim

  ! index mapping array
  integer,    allocatable :: a2vck(:,:)
  integer,    allocatable :: vck2a(:,:,:)
  private :: a2vck, vck2a

  ! data arrays
  real(8),    allocatable :: tvck(:)
  complex(8), allocatable :: wvck(:,:)
  complex(8), allocatable :: wvck0(:,:)
  public :: tvck, wvck

  real(8),    allocatable :: d(:)
  complex(8), allocatable :: md(:,:)
  complex(8), allocatable :: dmmd(:,:)
  complex(8), allocatable :: wij(:,:,:)
  private :: d, md, dmmd, wij

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

      ! step 2: diagonalize DMMD and calculate w_{vck} Eq. (2.20)
      call diagonalize_dmmd(iq)

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
  subroutine task_wpol()
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
    &                  kqset%nkpt, 1, &
    &                  iqstart, iqend)
    ! write(*,*) "myrank_row, iqstart, iqend =", myrank_row, iqstart, iqend
#else
    iqstart = 1
    iqend = kqset%nkpt
#endif

    ! Setup working array dimensions and index mappings
    call set_wpol_indices()

    !===========================
    ! Momentum matrix elements
    !===========================
    if (.not.input%gw%rpmat) call calcpmatgw()

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
      write(*,*) '(mod_wpol::task_wpol) q-point cycle, iq = ', iq
    
      Gamma = gammapoint(kqset%vqc(:,iq))

      ! Calculate interstitial product basis functions
      matsiz = locmatsiz+Gqset%ngk(1,iq)
      call diagsgi(iq)
      call calcmpwipw(iq)
    
      ! Calculate the bare Coulomb potential
      call calcbarcmb(iq)

      ! set v-diagonal MB
      call setbarcev(input%gw%barecoul%barcevtol)

      ! Calculate W_{ij} in pole representation
      call calc_md_dmmd(iq)
      call diagonalize_dmmd(iq)

      ! some clean up
      call delete_coulomb_potential
      call clear_wpol()

      ! output poles data
      call put_wpol(iq)

      ! clean unused data
      deallocate(tvck)
      deallocate(wvck)
      deallocate(mpwipw)
      deallocate(barc)

    end do ! q-points

    ! delete index mapping arrays
    call del_wpol_indices()

    return
  end subroutine 

!--------------------------------------------------------------------------------
  subroutine put_wpol(iq)
    use m_getunit
    implicit none
    integer,    intent(in) :: iq
    integer       :: fid, i
    integer(8)    :: recl
    character(80) :: fname

    ! store q-dependent data
    call getunit(fid)

    write(fname,'("TVCK-q",I4.4,".OUT")') iq
    open(fid, File=trim(fname), Action='WRITE', Err=10)
    write(fid,'("# ", i8)', Err=11) nvck
    do i = 1, nvck
      write(fid,'(i8, f16.6)', Err=11) i, tvck(i)
    end do
    close(fid)

    write(fname,'("WVCK-q",I4.4,".OUT")') iq
    inquire(IoLength=recl) mbdim, nvck, wvck(:,:)
    open(fid, File=trim(fname), Form='UNFORMATTED', Access='DIRECT', Recl=recl, Err=20)
    write(fid, Rec=1, Err=21) mbdim, nvck, wvck(:,:)
    close(fid)

    return

10  write(*,*) 'ERROR(mod_wpol::put_wpol) Error opening file TVCK-q.OUT for writing!'
    stop
11  write(*,*) 'ERROR(mod_wpol::put_wpol) Error writing data into TVCK-q.OUT!'
    stop
20  write(*,*) 'ERROR(mod_wpol::put_wpol) Error opening file WVCK-q.OUT for writing!'
    stop
21  write(*,*) 'ERROR(mod_wpol::put_wpol) Error writing data into WVCK-q.OUT!'
    stop

  end subroutine

!--------------------------------------------------------------------------------
  subroutine get_wpol(iq)
    use m_getunit
    implicit none
    integer,    intent(in) :: iq
    integer       :: fid, i, j
    integer(8)    :: recl
    character(80) :: fname, string

    ! store q-dependent data
    call getunit(fid)

    write(fname,'("TVCK-q",I4.4,".OUT")') iq
    open(fid, File=trim(fname), Action='READ', Err=10)
    read(fid,*) string, nvck
    
    if (allocated(tvck)) deallocate(tvck)
    allocate(tvck(nvck))
    do i = 1, nvck
      read(fid,*) j, tvck(i)
    end do
    close(fid)

    write(fname,'("WVCK-q",I4.4,".OUT")') iq
    inquire(IoLength=recl) mbdim, nvck
    open(fid, File=trim(fname), Form='UNFORMATTED', Access='DIRECT', Recl=recl, Err=20)
    read(fid, Rec=1) mbdim, nvck
    close(fid)  

    if (allocated(wvck)) deallocate(wvck)
    allocate(wvck(mbdim,nvck))

    inquire(IoLength=recl) mbdim, nvck, wvck(:,:)
    open(fid, File=trim(fname), Form='UNFORMATTED', Access='DIRECT', Recl=recl, Err=20)
    read(fid, Rec=1) mbdim, nvck, wvck(:,:)
    close(fid)

    return

10  write(*,*) 'ERROR(mod_wpol::put_wpol) Error opening file TVCK-q.OUT for reading!'
    stop
20  write(*,*) 'ERROR(mod_wpol::put_wpol) Error opening file WVCK-q.OUT for reading!'
    stop

  end subroutine

!--------------------------------------------------------------------------------
  subroutine calc_md_dmmd(iq)
    implicit none
    integer, intent(in) :: iq

    integer    :: ispn, ik, jk, ikp, jkp
    integer    :: ic, icg, ia, is, ias
    integer    :: i, n, m
    integer(8) :: recl
    real(8)    :: de
    real(8)    :: prefac
    complex(8) :: zt1, wkp
    complex(8), allocatable :: evecsv(:,:,:)

    mbdim = mbsiz
    if (Gamma) then      
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
      !
      select case(input%gw%scrcoul%sciavtype)
      case('isotropic')
        mbdim = mbsiz+1
      case('sphavrg')
        mbdim = mbsiz+3
      case default
        write(*,*) "ERROR(mod_wpol::calc_md_dmmd): Unknown averaging type!"
        stop
      end select  
    end if

    ! global arrays
    allocate(d(nvck0))
    d(:) = 0.d0
    allocate(md(mbdim,nvck0))
    md(:,:) = 0.d0  

    ! local data
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))    
    allocate(evecsv(nmatmax,nstsv,nspinor))
    allocate(minmmat(mbsiz,ndim,numin:nstdf))

    wkp = 2.d0/sqrt(dble(kqset%nkpt))
    prefac = sqrt(4.d0*pi/omega)*wkp

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
          if (Gamma) call q0_treatment_setup(i,prefac/sqrt(de),pmatvv(n,m,:))
        else
          icg = n-nomax
          is  = corind(icg,1)
          ia  = corind(icg,2)
          ic  = corind(icg,3)
          ias = idxas(ia,is)
          de  = evalsv(m,jkp)-evalcr(ic,ias)
          if (Gamma) call q0_treatment_setup(i,prefac/sqrt(de),pmatcv(icg,m,:))
        end if
        if (abs(de) < 1.d-6) then
          write(*,*) "WARNING(mod_wpol::calc_md_dmmd) Near degenerate eigenvalues ", n, m
          ! stop
        end if

        ! D matrix
        d(i) = de
        ! v^{1/2} M 2D^{1/2} matrix
        md(1:mbsiz,i) = minmmat(1:mbsiz,n,m)*sqrt(de)*wkp

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

    return
  end subroutine


!--------------------------------------------------------------------------------
  subroutine q0_treatment_setup(i,prefac,p)
    use modinput
    use modmain, only : pi, zi
    implicit none
    ! inpout
    integer,    intent(in) :: i
    real(8),    intent(in) :: prefac
    complex(8), intent(in) :: p(3)
    ! local
    real(8) :: c1, c2
    real(8) :: q0eps(3), modq0

    ! singular term: v^{1/2} M^{0} D^{1/2}

    select case(input%gw%scrcoul%sciavtype)
    case('isotropic')

      ! q->0 direction
      q0eps(:) = input%gw%scrcoul%q0eps(:)
      modq0    = sqrt(q0eps(1)**2+q0eps(2)**2+q0eps(3)**2)
      q0eps(:) = q0eps(:)/modq0

      md(mbsiz+1,i) = -prefac * ( p(1)*q0eps(1)+ &
      &                           p(2)*q0eps(2)+ &
      &                           p(3)*q0eps(3) )

    case('sphavrg')

      ! c1 = sqrt(4.d0*pi/3.d0)
      ! c2 = sqrt(2.d0*pi/3.d0)
      c1 = 1.d0
      c2 = sqrt(0.5d0)
      ! c1 = sqrt(1.d0/3.d0)
      ! c2 = sqrt(1.d0/6.d0)

      md(mbsiz+1,i) = -prefac * c2*(zi*p(2)+p(1))
      md(mbsiz+2,i) = -prefac * c1*p(3)
      md(mbsiz+3,i) = -prefac * c2*(zi*p(2)-p(1))

    case default
      write(*,*) "ERROR(mod_wpol::calc_md_dmmd): Unknown averaging type!"
      stop

    end select

    return
  end subroutine


!--------------------------------------------------------------------------------
  subroutine diagonalize_dmmd(iq)
    use mod_wpol_pert
    use mod_lanczos
    use mod_rank1_update, only : sevr, generate_update_vectors
    implicit none
    integer, intent(in) :: iq
    integer :: i, j, n, niter, blks
    real(8) :: rho
    real(8),    allocatable :: eval(:)
    complex(8), allocatable :: zmat(:,:), evec(:,:), z(:)
    complex(8), external    :: zdotc

    write(*,*)
    write(*,*) 'Info(mod_wpol::diagonalize_dmmd)'
    write(*,*) '    Matrix size: ', nvck0

    if (.not.associated(input%gw%eigensolver)) &
    &  input%gw%eigensolver => getstructeigensolver(emptynode)
   
    select case (input%gw%eigensolver%method)

      case ('pert')
        !----------------------
        ! Perturbation theory
        !----------------------
        write(*,*) '    (test) First-order perturbation theory'
        nvck = nvck0
        allocate(dmmd(nvck,nvck))
        call zgemm( 'c', 'n', nvck, nvck, mbdim, &
        &           zone, md, mbdim, md, mbdim,  &
        &           zzero, dmmd, nvck)
        if (Gamma) call add_q0_contribution()
        ! D*D + D^{1/2}*M^{+}*M*D^{1/2}
        do i = 1, nvck
          dmmd(i,i) = d(i)*d(i) + dmmd(i,i)
        end do
        !cutoff = input%gw%eigensolver%cutoff
		    cutoff=maxval(abs(dmmd))
        allocate(tvck(nvck))
        call wpol_pert(nvck,d,dmmd,tvck,cutoff)
        ! w_{vck} (2.20)
        allocate(wvck(mbdim,nvck))
        call zgemm( 'n', 'n', mbdim, nvck, nvck, &
        &           zone, md, mbdim, dmmd, nvck, &
        &           zzero, wvck, mbdim)
        deallocate(dmmd)
         
      case ('lanczos')
        !----------------------
        ! Band Lanczos method 
        !----------------------
        write(*,*) '    Band Lanczos algorithm'
        niter = min(input%gw%eigensolver%niter,nvck0)
        blks  = input%gw%eigensolver%blkSize
        nvck  = blks*niter

        allocate(evec(nvck0,nvck))
        evec(:,:) = zzero

        ! Starting basis estimate
        select case (input%gw%eigensolver%basis)

          case ('unit')
            do i = 1, min(nvck0, blks)
              evec(i,i) = zone
            end do

          case ('rand')
            do i = 1, blks
              evec(:,i) = md(i,:)
            end do
            call mkl_qr(nvck0, blks, evec(:,1:blks))

          case ('svd')
            allocate(zmat(mbdim,nvck0))
            zmat(:,:) = md(:,:)
            call orthogonalize(mbdim, nvck0, zmat)
            do i = 1, blks
              evec(:,i) = conjg(zmat(i,:))
            end do
            deallocate(zmat)

          case default
            write(*,*) 'ERROR(mod_wpol::diagonalize_dmmd): Unknown basis generator!'
            stop

        end select

        allocate(tvck(nvck))
        call lanczos_band_matvec(matrix_vector, nvck0, niter, blks, evec, tvck)
        
        ! w_{vck} (2.20)
        allocate(wvck(mbdim,nvck))
        call zgemm( 'n', 'n', mbdim, nvck, nvck0, &
        &           zone, md, mbdim, evec, nvck0, &
        &           zzero, wvck, mbdim)

        ! clear memory
        deallocate(evec)

      case ('lapack')
        !-------------------------
        ! Direct diagonalization
        !-------------------------
        write(*,*) '    LAPACK diagonalization'

        ! D*D + D^{1/2}*M^{+}*M*D^{1/2}
        nvck = nvck0
        allocate(dmmd(nvck,nvck))
        call zgemm( 'c', 'n', nvck, nvck, mbdim, &
        &           zone, md, mbdim, md, mbdim,  &
        &           zzero, dmmd, nvck)
        do i = 1, nvck
          dmmd(i,i) = d(i)*d(i) + dmmd(i,i)
        end do

        allocate(tvck(nvck))
        ! call mkl_zheev(nvck,dmmd,tvck)
        call mkl_zheevr(nvck,dmmd,tvck)

        ! w_{vck} (2.20)
        allocate(wvck(mbdim,nvck))
        call zgemm( 'n', 'n', mbdim, nvck, nvck, &
        &           zone, md, mbdim, dmmd, nvck, &
        &           zzero, wvck, mbdim)
        deallocate(dmmd)

      case ('rank1')
        !-------------------------
        ! Rank-1 update
        !-------------------------
        write(*,*) '    Rank-1 update diagonalization'

        ! SVD of M*D^{1/2}
        n = min(mbdim, nvck0)
        allocate(eval(n))
        allocate(evec(mbdim,nvck0))
        evec(:,:) = md(:,:)
        call generate_update_vectors(mbdim, nvck0, evec, eval)

        ! Rank-1 update iterations
        nvck = nvck0

        ! D*D
        allocate(tvck(nvck))
        do i = 1, nvck
          tvck(i) = d(i)*d(i)
        end do
        deallocate(d)

        allocate(dmmd(nvck,nvck))
        dmmd(:,:) = zzero
        do i = 1, nvck
          dmmd(i,i) = zone
        end do

        niter = min(input%gw%eigensolver%niter, n)
        allocate(z(nvck))
        do i = 1, niter
          rho = eval(i)*eval(i)
          call zgemv('c', nvck, nvck, zone, dmmd, nvck, &
          &           conjg(evec(i,:)), 1, zzero, z, 1)
          call sevr(nvck, tvck, z, rho, dmmd)
        end do
        deallocate(z, evec, eval)

        ! w_{vck} (2.20)
        allocate(wvck(mbdim,nvck))
        call zgemm( 'n', 'n', mbdim, nvck, nvck, &
        &           zone, md, mbdim, dmmd, nvck, &
        &           zzero, wvck, mbdim)
        
        deallocate(dmmd)

      case default
        write(*,*) "ERROR(mod_wpol::diagonalize_dmmd): Unknown eigensolver!"
        stop

    end select

    ! following definition of Eq.(2.20) tau = sqrt(lambda)
    do i = 1, nvck
      if (tvck(i) > 1.d-8) then
        tvck(i) = sqrt(tvck(i))
      else
        ! set to (numerical) zero
        tvck(i) = 0.d0
      end if
    end do

    return
  end subroutine

!--------------------------------------------------------------------
! Calculate product r = H*q = D*D*q + M^{+}*M*q + P^{+}*P*q (if q=0)
!--------------------------------------------------------------------
  subroutine matrix_vector(n, zq, zr)
    use modmain, only : zzero, zone    
    implicit none
    integer,    intent(in)  :: n
    complex(8), intent(in)  :: zq(n)
    complex(8), intent(out) :: zr(n)
    ! local
    integer :: i, j
    complex(8), allocatable :: zv(:), zp(:)
    real(8)    :: q0eps(3), modq0, c1, c2
    complex(8) :: p1(3), p2(3)
    complex(8), external    :: zdotu

    ! D*D*q
    do i = 1, n
      zr(i) = d(i)*d(i)*zq(i)
    end do

    ! M*q
    allocate(zv(mbdim))
    call zgemv('n', mbdim, n, zone, md, mbdim, zq, 1, zzero, zv, 1)
    ! M^{+}*M*q
    call zgemv('c', mbdim, n, zone, md, mbdim, zv, 1, zone,  zr, 1)
    deallocate(zv)

    return
  end subroutine


!--------------------------------------------------------------------------------
  subroutine calc_wmat()
    implicit none
    integer    :: iom, i, im
    real(8)    :: om, eta
    complex(8) :: zt1

    ! plot smearing
    eta = input%gw%scrcoul%swidth

    ! W_{ij}(q,\omega) (2.22) (regular part only)
    allocate(wij(mbsiz,mbsiz,freq%nomeg))

    do iom = 1, freq%nomeg

      om = freq%freqs(iom)

      do i = 1, nvck
        zt1 =  0.5d0 / tvck(i) *        &
        &    ( 1.d0/(om-tvck(i)+zi*eta) &
        &     -1.d0/(om+tvck(i)-zi*eta) )
        md(:,i) = zt1*wvck(:,i) ! reuse the array md
      end do

      ! v w_{vck}
      do im = 1, mbsiz
        md(im,:) = barcev(im)*md(im,:)
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

    nvck0 = kqset%nkpt*ndim*mdim

    ! map a -> {vck}
    if (allocated(a2vck)) deallocate(a2vck)
    allocate(a2vck(3,nvck0))

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
  subroutine clear_wpol()
    if (allocated(d))    deallocate(d)
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
    do im = 1, mbsiz, 10
      do iom = 1, nomeg
        write(fid,'(3f16.6)') freq%freqs(iom), zmat(im,im,iom)
      end do
      write(fid,*)
    end do
    close(fid)

    return
  end subroutine

end module

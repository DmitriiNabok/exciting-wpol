
module mod_wpol_selfc

  use modinput
  use modgw
  use mod_selfenergy, only : singc1, singc2, selfec, &
  &                          write_selfecnn
  use mod_mpi_gw
  use mod_wpol

  implicit none
  integer,    private :: mdim
  real(8),    private :: eta

  public  :: task_wpol_selfc
  private :: calc_selfc_wpol_q, calc_minmkq

  real(8), parameter :: eps = 1.d-8

contains

!--------------------------------------------------------------------------------
  subroutine task_wpol_selfc()
    use m_getunit
    implicit none
    integer :: iq, fid

    eta = input%gw%selfenergy%swidth

    !=================
    ! Initialization
    !=================
    call init_gw
    call clean_gndstate

    !------------------------
    ! total number of states
    !------------------------
    if (input%gw%coreflag=='all') then
      mdim = nstse+ncg
    else
      mdim = nstse
    end if

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

    ! q->0 singularity treatment scheme
    select case (trim(input%gw%selfenergy%singularity))
      case('none')
        singc1 = 0.d0
        singc2 = 0.d0
      case('mpb')
        ! Auxiliary function method
        call setsingc
      case('crg')  
        ! Auxiliary function method
        call calc_q0_singularities
      case default
        write(*,*) 'ERROR(mod_wpol::test_wpol) Unknown singularity treatment scheme!'
        stop
    end select

    !===========================
    ! Momentum matrix elements
    !===========================
    if (.not.input%gw%rpmat) call calcpmatgw

    !=========================
    ! Main loop over q-points
    !=========================

    if (myrank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if

    allocate(selfec(ibgw:nbgw,freq%nomeg,kset%nkpt))
    selfec(:,:,:) = zzero

    ! each process does a subset
    do iq = iqstart, iqend

      write(*,*)
      write(*,*) '(mod_selfc_wpol::test_selfc_wpol) q-point cycle, iq = ', iq
    
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
      ! call calc_wmat()
      ! call print_wmat(iq)
      call delete_coulomb_potential
      call clear_wpol()

      ! Calculate q-dependent \Sigma^c_{nn}(k,q;\omega)
      call calc_selfc_wpol_q(iq)

      ! clean unused data
      deallocate(tvck)
      deallocate(wvck)
      deallocate(mpwipw)
      deallocate(barc)

    end do ! q-points

#ifdef MPI
      call mpi_sum_array(0,selfec,nbandsgw,freq%nomeg,kset%nkpt,mycomm_row)
#endif

    ! print to file the results
    if (myrank==0) call write_selfecnn()

    ! delete index mapping arrays
    call del_wpol_indices()
    deallocate(selfec)
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine calc_selfc_wpol_q(iq)
    implicit none
    integer, intent(in) :: iq
    integer    :: ispn, ik, jk, jkp, n, iom, m, i
    integer(8) :: recl
    real(8)    :: wkq
    complex(8) :: zt1, zwt, coefs1, coefs2, sigma1, sigma2
    complex(8), allocatable :: mw(:,:)

    wkq    = 1.d0 / dble(kqset%nkpt)
    coefs1 = singc1 * sqrt(4.d0*pi/omega)
    coefs2 = singc2 * 4.d0*pi/omega

    allocate(mw(mdim,nvck))
    allocate(minmmat(1:mbsiz,ibgw:nbgw,1:mdim))

    ! loop over k-points
    do ispn = 1, nspinor
    do ik = 1, kset%nkpt

      write(*,*)
      write(*,'(a,3i8)') '(mod_selfc_wpol::calc_selfc_wpol_q): rank, (iq, ik):', myrank, iq, ik

      ! k-q point
      jk = kqset%kqid(ik,iq)
      jkp = kset%ik2ikp(jk)

      ! Product basis coefficients M^i_{nm}(k,q)
      call calc_minmkq(ik,iq,ispn)

      ! loop over states (diagonal matrix elements only)
      do n = ibgw, nbgw

        ! precalculate \sum_{i} conjg(M^i_{nm})*w_{vck} in Eq.(2.28)
        call zgemm( 'c', 'n', mdim, nvck, mbsiz, &
        &           zone, minmmat(1:mbsiz,n,:), mbsiz, wvck(1:mbsiz,:), mbsiz, &
        &           zzero, mw, mdim)

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iom,zt1)
!$OMP DO
#endif
        ! loop over frequencies
        do iom = 1, freq%nomeg
          zt1 = 0.d0
          ! sum over states
          zt1 = zt1 + sum_valence(jkp,iom,mw)
          zt1 = zt1 + sum_core(jkp,iom,mw)
          selfec(n,iom,ik) = selfec(n,iom,ik) + wkq*zt1
        end do ! iom
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

        if (Gamma) then
          ! loop over frequencies
          do iom = 1, freq%nomeg
            call sum_singular(n,jkp,iom,mw,sigma1,sigma2)
            selfec(n,iom,ik) = selfec(n,iom,ik)+ &
            &                  coefs1*sigma1   + &
            &                  coefs2*sigma2
          end do
        end if

      end do ! n

    end do ! ik
    end do ! nspinor
     
    deallocate(minmmat)
    deallocate(mw)

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine sum_singular(n,jkp,iom,mw,sigma1,sigma2)
    implicit none
    integer,    intent(in) :: n, jkp, iom
    complex(8), intent(in) :: mw(mdim,nvck)
    complex(8)             :: sigma1, sigma2
    ! local
    integer    :: m, i
    real(8)    :: enk, om
    real(8)    :: x, y, t1, t2
    complex(8) :: zt1, zt2, p(3)
    real(8)    :: q0eps(3), modq0
    real(8)    :: c1, c2

    enk = evalsv(n,jkp)

    om = freq%freqs(iom)

    sigma1 = zzero
    sigma2 = zzero

    select case(input%gw%scrcoul%sciavtype)
    
    case('isotropic')

      ! q->0 direction
      q0eps(:) = input%gw%scrcoul%q0eps(:)
      modq0    = sqrt(q0eps(1)**2+q0eps(2)**2+q0eps(3)**2)
      q0eps(:) = q0eps(:)/modq0

#ifdef USEOMP
!$omp parallel &
!$omp default(shared) &
!$omp private(i,zt1,zt2)
!$omp do reduction (+:sigma1,sigma2)
#endif
      do i = 1, nvck
        !------------------
        zt1 = tvck(i) * ( om - enk + sign(1,nomax-n)*(tvck(i)-zi*eta) )
        zt1 = 0.5d0 / zt1
        !------------------
        ! t1 = tvck(i) * ( om - enk + sign(1,nomax-n)*tvck(i) )
        ! t2 = tvck(i) * sign(1,nomax-n) * eta
        ! zt1 = 0.5d0*cmplx( t1 / ( t1*t1 + t2*t2 ), &
        ! &                  t2 / ( t1*t1 + t2*t2 ), 8)
        !------------------
        zt2 = wvck(mbsiz+1,i)
        !---------------------------------------------------
        ! contribution from the first term: 1/q^2
        sigma2 = sigma2 + zt1 * zt2*conjg(zt2)
        !---------------------------------------------------
        ! contribution from the second+third term: 1/q
        sigma1 = sigma1 + zt1 * ( zt2*conjg(mw(n,i)) + conjg(zt2)*mw(n,i) )
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif      

    case('sphavrg')

      c1 = sqrt(1.d0/3.d0)
      c2 = sqrt(1.d0/6.d0)
#ifdef USEOMP
!$omp parallel &
!$omp default(shared) &
!$omp private(i,zt1,p,zt2)
!$omp do reduction (+:sigma1,sigma2)
#endif      
      do i = 1, nvck
        zt1 = tvck(i) * ( om - enk + sign(1,nomax-n)*(tvck(i)-zi*eta) )
        zt1 = 0.5d0 / zt1
        p(1) = c2*(zi*wvck(mbsiz+2,i)+wvck(mbsiz+1,i))
        p(2) = c1*wvck(mbsiz+3,i)
        p(3) = c2*(zi*wvck(mbsiz+2,i)-wvck(mbsiz+1,i))
        zt2 = ( p(1)*conjg(p(1))+ &
        &       p(2)*conjg(p(2))+ &
        &       p(3)*conjg(p(3)) )
        !---------------------------------------------------
        ! contribution from the first term: 1/q^2
        sigma2 = sigma2 + zt1 * zt2
        !---------------------------------------------------
        ! There is no contribution from the second+third term: 1/q
        sigma1 = zzero
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif

    case default
      write(*,*) "ERROR(mod_wpol_selfc::sum_singular): Unknown averaging type!"
      stop
    end select

    return
  end subroutine

!--------------------------------------------------------------------------------
  function sum_valence(jkp,iom,mw) result(zsum)
    implicit none
    integer,    intent(in) :: jkp, iom
    complex(8), intent(in) :: mw(mdim,nvck)
    complex(8)             :: zsum
    ! local
    integer    :: m, i
    real(8)    :: enk, om
    real(8)    :: x, y, t1, t2
    complex(8) :: zt1
    complex(8), allocatable :: mwt(:)
    complex(8), external :: zdotc

    zsum = 0.d0
    allocate(mwt(nvck))
    mwt(:) = 0.d0

    om = freq%freqs(iom)

    ! sum over states
    do m = 1, nstse
      enk = evalsv(m,jkp)
      ! apply frequency/state dependent prefactor
      do i = 1, nvck
        !------------------
        zt1 = tvck(i) * ( om - enk + sign(1,nomax-m)*(tvck(i)-zi*eta) )
        zt1 = 0.5d0 / zt1
        !------------------
        ! t1 = tvck(i) * ( om - enk + sign(1,nomax-m)*tvck(i) )
        ! t2 = tvck(i) * sign(1,nomax-m) * eta
        ! zt1 = 0.5d0*cmplx( t1 / ( t1*t1 + t2*t2 ), &
        ! &                  t2 / ( t1*t1 + t2*t2 ), 8)
        !------------------
        mwt(i) = zt1*mw(m,i)
      end do
      ! sum over vck
      zsum = zsum + zdotc(nvck,mw(m,:),1,mwt,1)
    end do ! m
    deallocate(mwt)

    return
  end function

!--------------------------------------------------------------------------------
  function sum_core(jkp,iom,mw) result(zsum)
    implicit none
    integer,    intent(in) :: jkp, iom
    complex(8), intent(in) :: mw(mdim,nvck)
    complex(8)             :: zsum
    ! local
    integer    :: m, i
    integer    :: icg, is, ia, ias, ic
    real(8)    :: enk, om
    real(8)    :: x, y, t1, t2
    complex(8) :: zt1
    complex(8), allocatable :: mwt(:)
    complex(8), external    :: zdotc, zdotu

    om = freq%freqs(iom)

    zsum = 0.d0
    allocate(mwt(nvck))
    ! sum over states
    do m = nstse+1, mdim
      icg = m-nstse
      is  = corind(icg,1)
      ia  = corind(icg,2)
      ias = idxas(ia,is)
      ic  = corind(icg,3)
      enk = evalcr(ic,ias)
      ! apply frequency/state dependent prefactor
      do i = 1, nvck
        zt1 = tvck(i) * ( om - enk + (tvck(i)-zi*eta) )
        zt1 = 0.5d0 / zt1
        mwt(i) = zt1*mw(m,i)
      end do
      ! sum over vck
      zsum = zsum + zdotc(nvck,mw(m,:),1,mwt,1)
    end do ! m
    deallocate(mwt)

    return
  end function

!--------------------------------------------------------------------------------
  subroutine calc_minmkq(ik,iq,ispn)
    implicit none
    integer, intent(in) :: ik, iq, ispn
    integer :: jk, im
    complex(8), allocatable :: evecsv(:,:,:)

    ! k-q point
    jk = kqset%kqid(ik,iq)

    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))    
      
    ! get KS eigenvectors
    allocate(evecsv(nmatmax,nstsv,nspinor))
    call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
    eveckp = conjg(evecsv(:,:,ispn))
    call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
    eveck = evecsv(:,:,ispn)
    deallocate(evecsv)
        
    ! Calculate M^i_{nm}+M^i_{nc}
    call expand_evec(ik,'t')
    call expand_evec(jk,'c')

    call expand_products(ik,iq,ibgw,nbgw,-1,1,mdim,nstse,minmmat)
    
    deallocate(eveckalm)
    deallocate(eveckpalm)
    deallocate(eveck)
    deallocate(eveckp)

    return
  end subroutine

end module
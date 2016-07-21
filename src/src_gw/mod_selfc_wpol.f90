
module mod_selfc_wpol

  use modinput
  use modgw
  use mod_mpi_gw
  use mod_wpol

  implicit none
  integer,    private :: mdim
  complex(8), private :: zieta, zif

  public  :: test_selfc_wpol, print_selfc_wpol
  private :: calc_selfc_wpol_q, calc_minmkq, omegaWeight

contains


!--------------------------------------------------------------------------------
  subroutine test_selfc_wpol()
    implicit none
    integer :: iq

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
    write(*,*) "myrank_row, iqstart, iqend =", myrank_row, iqstart, iqend
#else
    iqstart = 1
    iqend = kqset%nkpt
#endif

    ! Setup working array dimensions and index mappings
    call set_wpol_indices()

    zieta = zi*input%gw%selfenergy%swidth
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
    if (.not.input%gw%rpmat) then
      call calcpmatgw
    end if

    !=========================
    ! Main loop over q-points
    !=========================

    if (myrank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if

    allocate(selfec(ibgw:nbgw,freq%nomeg,kset%nkpt))
    selfec(:,:,:) = 0.d0

    ! each process does a subset
    do iq = iqstart, iqend

      write(*,*)
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
      call calc_wvck(iq)
      ! call calc_wmat()
      ! call print_wmat(iq)
      call delete_coulomb_potential
      call clear_wpol()

      ! Calculate q-dependent \Sigma^c_{nn}(k,q;\omega)
      call calc_selfc_wpol_q(iq)
      ! write(*,*) 'selfec_q=', sum(selfec)

      ! clean unused data
      deallocate(tvck)
      deallocate(wvck)
      deallocate(mpwipw)
      deallocate(barc)

    end do ! q-points

#ifdef MPI
    if ((nproc_row>1).and.(myrank_col==0)) then
      call mpi_sum_array(0,selfec,nbandsgw,freq%nomeg,kset%nkpt,mycomm_row)
      write(*,*) "sum selfec done"
    end if
#endif

    ! print to file the results
    if (myrank==0) call print_selfc_wpol(selfec,1,ibgw,nbgw)

    ! delete index mapping arrays
    call del_wpol_indices()
    deallocate(selfec)

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine calc_selfc_wpol_q(iq)
    implicit none
    integer, intent(in) :: iq
    integer    :: ispn, ik, jk, jkp, n, iom, m, i
    integer(8) :: recl
    real(8)    :: wkq
    complex(8) :: zt1, zwt, coefs1, coefs2
    ! complex(8), allocatable :: wght(:,:,:)
    complex(8), allocatable :: mw(:,:), mwt(:), zwt0(:)

    complex(8), external    :: zdotc, zdotu

    wkq    = 1.d0 / dble(kqset%nkpt)
    coefs1 = singc1 * sqrt(4.d0*pi/omega)
    coefs2 = singc2 * 4.d0*pi/omega

    ! allocate(wght(nvck,mdim,freq%nomeg))
    allocate(mw(mdim,nvck),mwt(nvck))
    if (Gamma) allocate(zwt0(nvck))
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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iom,zt1,m,i,zwt,mwt,zwt0)
!$OMP DO
#endif
        ! loop over frequencies
        do iom = 1, freq%nomeg

          ! sum over states
          zt1 = 0.d0
          do m = 1, mdim
            ! apply frequency/state dependent prefactor
            do i = 1, nvck
              zwt = omegaWeight(jkp,i,m,iom)
              mwt(i) = zwt*mw(m,i)
              if ((Gamma) .and. (m==n)) zwt0(i) = zwt
            end do
            ! sum over vck
            zt1 = zt1 + zdotc(nvck,mw(m,:),1,mwt,1)
          end do ! m
          selfec(n,iom,ik) = selfec(n,iom,ik)+wkq*zt1

          if (Gamma) then
            !---------------------------------------------------
            ! contribution from the first term: 1/q^2
            do i = 1, nvck
              mwt(i) = zwt0(i)*wvck(mbsiz+1,i)
            end do
            zt1 = coefs2*zdotc(nvck,wvck(mbsiz+1,:),1,mwt,1)
            selfec(n,iom,ik) = selfec(n,iom,ik)+zt1
            !---------------------------------------------------
            ! contribution from the second term: 1/q
            do i = 1, nvck
              mwt(i) = zwt0(i)*wvck(mbsiz+1,i)
            end do
            zt1 = coefs1*zdotc(nvck,mw(n,:),1,mwt,1)
            selfec(n,iom,ik) = selfec(n,iom,ik)+zt1
            !---------------------------------------------------
            ! contribution from the third term: 1/q
            do i = 1, nvck
              mwt(i) = zwt0(i)*conjg(wvck(mbsiz+1,i))
            end do
            zt1 = coefs1*zdotu(nvck,mw(n,:),1,mwt,1)
            selfec(n,iom,ik) = selfec(n,iom,ik)+zt1
          end if
          
        end do ! iom
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

      end do ! n

    end do ! ik
    end do ! nspinor
     
    deallocate(minmmat)
    deallocate(mw,mwt)
    if (Gamma) deallocate(zwt0)

    return
  end subroutine


!--------------------------------------------------------------------------------
  complex(8) function omegaWeight(jkp,i,m,iom)
    implicit none
    integer,    intent(in)  :: jkp, i, m, iom
    ! local
    integer    :: icg, is, ia, ias, ic
    real(8)    :: enk, sgn
    complex(8) :: zt1

    if (m <= nstse) then
      enk = evalsv(m,jkp)
    else
      icg = m-nstse
      is  = corind(icg,1)
      ia  = corind(icg,2)
      ias = idxas(ia,is)
      ic  = corind(icg,3)
      enk = evalcr(ic,ias)
    end if

    if ((m <= nomax).or.(m > nstse)) then
      ! valence or core state
      sgn = 1.d0
    else
      ! conduction state
      sgn = -1.d0
    end if

    zt1 = tvck(i) * ( zif*freq%freqs(iom) -(enk-efermi) + sgn*(tvck(i)-zieta) )
    omegaWeight = 0.5d0 / zt1

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


!--------------------------------------------------------------------------------
  subroutine print_selfc_wpol(zmat,ik,ib,nb)
    use m_getunit
    implicit none
    ! complex(8), intent(in) :: zmat(:,:,:)
    complex(8), intent(in) :: zmat(ibgw:nbgw,freq%nomeg,kset%nkpt)
    integer,    intent(in) :: ik
    integer,    intent(in) :: ib, nb

    integer :: iom, n
    integer :: fid1, fid2
    character(80) :: frmt

    ! write(*,*) 'SIZE=', size(zmat,dim=1), size(zmat,dim=2), size(zmat,dim=3)

    call getunit(fid1)
    open(fid1, File='SELFC-WPOL-Re.OUT', Action='WRITE')
    call getunit(fid2)
    open(fid2, File='SELFC-WPOL-Im.OUT', Action='WRITE')

    n = nb-ib+1
    write(frmt,'("(",i8,"f14.6)")') 1+n
    ! write(*,*) trim(frmt)

    write(fid1,*) '# ik = ', ik
    write(fid2,*) '# ik = ', ik
    do iom = 1, freq%nomeg
      write(fid1,trim(frmt)) freq%freqs(iom), dble(zmat(ib:nb,iom,ik))
      write(fid2,trim(frmt)) freq%freqs(iom), imag(zmat(ib:nb,iom,ik))
    end do

    close(fid1)
    close(fid2)

    return
  end subroutine

end module
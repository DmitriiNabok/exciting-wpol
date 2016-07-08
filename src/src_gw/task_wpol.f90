      
subroutine task_wpol()

    use modinput
    use modmain,               only : zzero, zone, zi, evalsv, efermi
    use modgw
    use mod_mpi_gw
    use m_getunit
            
    implicit none
    integer :: fid
    integer :: iq, ik, jk, ikp, jkp
    integer :: i, n, m, iom
    real(8) :: de, t1
    complex(8) :: zt1, zwkp
    integer :: ndim, mdim, nmdim
    integer :: ispn
    integer :: im, nvck
    integer :: ic, icg, ia, is, ias
    
    character(80) :: fname
    
    complex(8) :: zieta, zif
    
    ! mapping array 
    integer,    allocatable :: a2vck(:,:)
    integer,    allocatable :: vck2a(:,:,:)

    complex(8), allocatable :: evecsv(:,:,:)

    real(8),    allocatable :: d(:)
    complex(8), allocatable :: md(:,:)
    real(8),    allocatable :: dmmdev(:)
    complex(8), allocatable :: dmmd(:,:)
    complex(8), allocatable :: wvck(:,:)
    complex(8), allocatable :: wij(:,:,:)

    ! .. Local Scalars ..
    integer :: info, lwork
    !
    ! .. Local Arrays ..
    !    RWORK dimension should be at least MAX(1,3*N-2)
    real(8),    allocatable :: rwork(:)
    complex(8), allocatable :: work(:)

!===========================================================================
! Initialization
!===========================================================================
    
    ! prepare GW global data
    call init_gw
    
    ! clean not used anymore global exciting variables
    call clean_gndstate
    
!===========================================================================
! Main loop: BZ integration
!===========================================================================    

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
    iomstart = 1
    iomend = freq%nomeg

    if (myrank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if

    ! total number of states including the core ones
    if (input%gw%coreflag=='all') then
      ndim = nomax+ncg
    else
      ndim = nomax
    end if
    mdim = nstdf-numin+1
    nmdim = ndim*mdim
    
    ! arrays to store products of KS eigenvectors with the matching coefficients
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))    

    nvck = kqset%nkpt*nomax*(nstdf-numin+1)
    call indexMapping()
    allocate(dmmdev(nvck))
    allocate(dmmd(nvck,nvck))

    ! smearing parameter
    zieta = zi*1.d-8

    select case (input%gw%freqgrid%fconv)
      case('refreq')
        ! real axis
        zif = zone
      case('imfreq')
        ! imaginary axis
        zif = zi
      case default
        write(*,*) "Correct options: refreq or imfreq"
        stop
    end select
    
    ! each process does a subset
    do iq = iqstart, iqend
    
      write(*,*)
      write(*,*) '(task_wpol): q-point cycle, iq = ', iq
    
      Gamma = gammapoint(kqset%vqc(:,iq))
          
      !========================================
      ! Calculate interstitial basis functions
      !========================================
      matsiz = locmatsiz+Gqset%ngk(1,iq)
      call diagsgi(iq)
      call calcmpwipw(iq)
    
      !======================================
      ! Calculate the bare Coulomb potential
      !======================================
      call calcbarcmb(iq)
      ! Set v-diagonal MB
      call setbarcev(input%gw%barecoul%barcevtol)
      
      allocate(minmmat(1:mbsiz,1:ndim,numin:nstdf))
      allocate(d(1:nvck))
      allocate(md(1:mbsiz,1:nvck))

      ! Loop over k-points
      do ispn = 1, nspinor
      do ik = 1, kqset%nkpt
        
        write(*,*)
        write(*,*) 'task_wpol: rank, (iq, ik):', myrank, iq, ik
    
        ! k-q point
        jk = kqset%kqid(ik,iq)

        ! irreducible k-point index
        ikp = kset%ik2ikp(ik)
        jkp = kset%ik2ikp(jk)
      
        ! get KS eigenvectors
        allocate(evecsv(nmatmax,nstsv,nspinor))
        call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
        eveckp = conjg(evecsv(:,:,ispn))
        call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
        eveck = evecsv(:,:,ispn)
        deallocate(evecsv)
        
        ! Calculate M^i_{nm}+M^i_{cm}
        call expand_evec(ik,'t')
        call expand_evec(jk,'c')
        call expand_products(ik,iq,1,ndim,nomax,numin,nstdf,-1,minmmat)
        
        ! Setup \tilde{M}_{ab}(q) * 2D^{1/2}
        do m = numin, nstdf
        do n = 1, ndim
          if (n <= nomax) then
            de = evalsv(m,jkp)-evalsv(n,ikp)
          else
            icg = n-nomax
            is  = corind(icg,1)
            ia  = corind(icg,2)
            ic  = corind(icg,3)
            ias = idxas(ia,is)
            de  = evalsv(m,jkp)-evalcr(ic,ias)
          end if
          if (abs(de) < 1.d-6) then
            write(*,*) "Suspicious degeneracy! ", n, m
            stop
          end if
          i = vck2a(n,m,ik)
          d(i) = de
          md(:,i) = sqrt(de)*minmmat(:,n,m)
        end do
        end do

      end do ! ik
      end do ! nspinor

      ! clean unused data
      deallocate(minmmat)
      deallocate(mpwipw)
      deallocate(barc)
      deallocate(vmat)

      ! 2D^{1/2} M^{+} v M 2D^{1/2} = 2D^{1/2} \tilde{M}^{+} \tilde{M} 2D^{1/2}
      call zgemm( 'c','n',nvck,nvck,mbsiz, &
      &           zone,md,mbsiz,md,mbsiz, &
      &           zzero,dmmd,nvck)

      ! account for prefactor 4/Nk
      dmmd(:,:) = dmmd(:,:)*4.d0/dble(kqset%nkpt)

      ! D^2 + D^{1/2} M^{+} v M D^{1/2}
      do i = 1, nvck
        dmmd(i,i) = d(i)*d(i) + dmmd(i,i)
      end do
      deallocate(d)

      write(*,*)
      write(*,*) 'Diagonalization with zheev'
      write(*,*)
      lwork = 2*nvck
      allocate(work(lwork),rwork(3*nvck))
      call zheev( 'v', 'u', nvck, dmmd, nvck, &
      &           dmmdev, work, lwork, rwork, info )
      deallocate(work,rwork)

      ! check for convergence
      if ( info > 0 ) then
         write(*,*)'Error(task_wpol): the algorithm failed to compute eigenvalues'
         write(*,*)'info = ', info
         stop
      end if

      ! v * M * D^{1/2} = v^{1/2} * \tilde{M} * 2D^{1/2}
      do im = 1, mbsiz
        md(im,:) = sqrt(barcev(im))*md(im,:)
      end do

      ! w_{vck} (2.20)
      allocate(wvck(mbsiz,nvck))
      call zgemm( 'n','n',mbsiz,nvck,nvck, &
      &           zone,md,mbsiz,dmmd,nvck, &
      &           zzero,wvck,mbsiz)

      !============================
      ! W_{ij}(q,\omega) (2.22)
      !============================
      allocate(wij(mbsiz,mbsiz,iomstart:iomend))
      do iom = iomstart, iomend

        do i = 1, nvck
          if (abs(dmmdev(i)) > 1.d-8) then
            t1 = sqrt(dmmdev(i))
            !-----------------------------------------------
            ! smearing method
            zt1 =  0.5d0 / t1 *                     &
            &    ( 1.d0/(zif*freq%freqs(iom)-t1+zieta) &
            &     -1.d0/(zif*freq%freqs(iom)+t1-zieta) )
            !-----------------------------------------------
            md(:,i) = zt1*wvck(:,i) ! reuse the array md
          else
            md(:,i) = 0.d0
          end if
        end do

        call zgemm( 'n','c',mbsiz,mbsiz,nvck, &
        &           zone,md,mbsiz,wvck,mbsiz, &
        &           zzero,wij(:,:,iom),mbsiz)

        ! W = v + wij
        ! do im = 1, mbsiz
        !   wij(im,im,iom) = barcev(im)+wij(im,im,iom)
        ! end do

      end do ! iom

      ! account for prefactor 4/Nk
      wij = wij*4.d0/dble(kqset%nkpt)

      ! store q-dependent Wij
      call getunit(fid)
      write(fname,'("WMAT-mat-q",I4.4,".OUT")') iq
      open(fid, File=trim(fname), Action='WRITE')
      do iom = 1, freq%nomeg, 100
      do im = 1, mbsiz
        write(fid,'(i8, 2f16.6)') im, wij(im,im,iom)
      end do
      write(fid,*)
      end do
      close(fid)

      write(fname,'("WMAT-iom-q",I4.4,".OUT")') iq
      open(fid, File=trim(fname), Action='WRITE')
      do im = mbsiz-100, mbsiz, 10
      do iom = iomstart, iomend
        write(fid,'(i8, 2f16.6)') iom, wij(im,im,iom)
      end do
      write(fid,*)
      end do
      close(fid)

      deallocate(md)
      deallocate(wvck)
      deallocate(wij)

    end do ! iq

    deallocate(dmmdev)
    deallocate(dmmd)

    ! clean up
    deallocate(eveckalm)
    deallocate(eveckpalm)
    deallocate(eveck)
    deallocate(eveckp)

    if (allocated(kiw)) deallocate(kiw)
    if (allocated(ciw)) deallocate(ciw)
    
    if (allocated(evalsv)) deallocate(evalsv)
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
    
    return

contains

  subroutine indexMapping()
    implicit none

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

end subroutine

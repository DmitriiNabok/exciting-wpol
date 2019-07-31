module mod_non_diagonal_wpol

  use modinput
  use modgw
  use modmain,        only : zzero, evalsv, efermi
  use mod_selfenergy, only : singc1, singc2
  use mod_selfc_wpol, only : mdim, calc_minmkq, calc_q_surr, norm2, cross, tvck_s, q_surr
  use mod_wpol
  use mod_mpi_gw
  use mod_hdf5
  
  implicit none

  real(8),    private :: eta, etasign, wkq, REEQP, IMEQP, f_ac, wkq_eta

  complex(8), allocatable, public :: selfec_nn(:,:),selfex_nn(:,:,:),vxc_nn(:,:,:)

  !valence-valence momentum matrix elements
  complex(8), allocatable :: pmatvv_k(:,:,:,:)
  ! core-valence momentum matrix elements
  complex(8), allocatable :: pmatcv_k(:,:,:,:)
  
  public  :: calc_selfc_nn_wpol, calc_selfx_nn, calc_vxc_nn
  private :: calc_selfc_nn_wpol_q, calc_selfx_nn_q
  
  real(8), parameter :: eps = 1.d-8

contains

subroutine calc_selfx_nn()

  use m_getunit
            
  implicit none
  integer(4) :: ikp, iq, fid, ik
  integer(4) :: recl
  !integer :: ie

  !===========================================================================
  ! Initialization
  !===========================================================================
  !call init_gw
  !call clean_gndstate
    
  ! occupancy dependent BZ integration weights
  call kintw

  !===========================================================================
  ! Main loop: BZ integration
  !===========================================================================    
	
#ifdef MPI
  call set_mpi_group(kqset%nkpt)
  call mpi_set_range(nproc_row, &
  &                  myrank_row, &
  &                  kqset%nkpt, 1, &
  &                  iqstart, iqend)
#else
  iqstart = 1
  iqend = kqset%nkpt
#endif

  !---------------------------------------
  ! treatment of singularities at G+q->0
  !---------------------------------------
  singc1 = 0.d0
  singc2 = 0.d0
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
      write(*,*) 'ERROR(task_gw): Unknown singularity treatment scheme!'
      stop
  end select

  ! each process does a subset
  do iq = iqstart, iqend
    
    write(*,*)
    write(*,*) '(mod_non_diagonal_selfx_nn): q-point cycle, iq = ', iq
    
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

    !===============================
    ! Calculate \Sigma^{x}_{kn}(q)
    !===============================
    call calc_selfx_nn_q(iq)
      
    ! clean unused data
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
      
  end do ! iq
    
  if (allocated(kiw)) deallocate(kiw)
  if (allocated(ciw)) deallocate(ciw)
    
#ifdef MPI
    call mpi_sum_array(0,selfex_nn,nbandsgw,nbandsgw,kset%nkpt,mycomm_row)
#endif
    
  !if (myrank==0) call write_selfexnn()
  
  ! clear memory  
  ! call delete_freqgrid(freq)
  ! call delete_k_vectors(kset)
  ! call delete_G_vectors(Gset)
  ! call delete_Gk_vectors(Gkset)
  ! call delete_kq_vectors(kqset)
  ! call delete_Gk_vectors(Gqset)
  ! call delete_Gk_vectors(Gqbarc)
  
  return
end subroutine

!--------------------------------------------------------------------------------
subroutine calc_selfx_nn_q(iq)

!!USES:
    use modmain,        only : nstfv, apwordmax, lmmaxapw, natmtot, nspnfv, &
    &                          pi, idxas, zzero, nmatmax, zone, occsv
      
!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq
    
!!LOCAL VARIABLES:            
    integer(4) :: ik, ikp, jk
    integer(4) :: mdim, nmdim
    real(8)    :: tstart, tend, t0, t1
    integer(4) :: ie1,ie1_, ie2, im
    integer(4) :: ia, is, ias, ic, icg
    real(8)    :: wkq, sxs2, de, de_ , q0eps(3), modq0
    complex(8) :: sx, vc
    complex(8) :: mvm     ! Sum_ij{M^i*V^c_{ij}*conjg(M^j)}
    complex(8), allocatable :: evecsv(:,:,:)
   
    integer :: k, l, ispn, ist, jst
    complex(8) :: zsum
    
    ! external routines 
    complex(8), external :: zdotc    
    
    call timesec(tstart)
   
    if (vccut) then
      sxs2 = 0.d0
      !----------------------------------------
      ! Set v-diagonal mixed product basis set
      !----------------------------------------
      mbsiz = matsiz
      if (allocated(barc)) deallocate(barc)
      allocate(barc(matsiz,mbsiz))
      barc(:,:) = zzero
      do im = 1, matsiz
        vc = cmplx(barcev(im),0.d0,8)
        barc(:,im) = vmat(:,im)*sqrt(vc)
      end do
    else
      ! singular term prefactor (q->0)
      sxs2 = -4.d0*pi*vi
      !----------------------------------------
      ! Set v-diagonal mixed product basis set
      !----------------------------------------
      call setbarcev(0.d0)
    end if
    
    !--------------------------------------------------
    ! total number of states (n->m + n->c transisions)
    !--------------------------------------------------
    if ((input%gw%coreflag=='all').or. &
    &   (input%gw%coreflag=='xal')) then
      mdim = nomax+ncg
    else
      mdim = nomax
    end if
    nmdim = (nbgw-ibgw+1)*mdim
    
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))
    
    allocate(minmmat(mbsiz,ibgw:nbgw,1:mdim))
    minmmat(:,:,:) = zzero
    msize = sizeof(minmmat)*b2mb
    ! write(*,'(" calcselfx: rank, size(minmmat) (Mb):",i4,f12.2)') myrank, msize
    
    !================================
    ! loop over irreducible k-points
    !================================
    do ispn = 1, nspinor
    do ikp = 1, kset%nkpt
    
      ! write(*,*)
      ! write(*,*) '(calcselfx): k-point loop ikp=', ikp
    
      ! k vector
      ik = kset%ikp2ik(ikp)
      ! k-q vector 
      jk = kqset%kqid(ik,iq)
      
      ! get KS eigenvectors
      call timesec(t0)
      allocate(evecsv(nmatmax,nstsv,nspinor))
      call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
      !call getevecfv(kqset%vkl(:,jk),Gkset%vgkl(:,:,:,jk),evecsv)
      eveckp = conjg(evecsv(:,:,ispn))
      
      call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
      !call getevecfv(kqset%vkl(:,ik),Gkset%vgkl(:,:,:,ik),evecsv)
      eveck = evecsv(:,:,ispn)
      deallocate(evecsv)
      call timesec(t1)
      time_io = time_io+t1-t0

      ! Calculate M^i_{nm}+M^i_{cm}
      call expand_evec(ik,'t')
      call expand_evec(jk,'c')
      call expand_products(ik,iq,ibgw,nbgw,-1,1,mdim,nomax,minmmat)
      
      ! q->0 direction
      q0eps(:) = input%gw%scrcoul%q0eps(:)
      modq0    = sqrt(q0eps(1)**2+q0eps(2)**2+q0eps(3)**2)
      q0eps(:) = q0eps(:)/modq0
      

      !========================================================
      ! Calculate the contribution to the exchange self-energy
      !========================================================
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ie1,ie1_,ie2,sx,mvm,icg,ia,is,ic,ias)
!$OMP DO
#endif    
      do ie1 = ibgw, nbgw
       do ie1_ = ibgw, nbgw   

        ! sum over occupied states
        sx = zzero 
        
        do ie2 = 1, mdim
          !======================= 
          ! Valence contribution
          !======================= 
          if (ie2 <= nomax) then
            mvm = zdotc(mbsiz,minmmat(:,ie1,ie2),1,minmmat(:,ie1_,ie2),1)
            sx = sx-kiw(ie2,jk)*mvm

           !momentum matrix term for n,n' and q->0 ; p_nv q * p_n'v q terms
		   ! de = evalsv(ie2,ikp)-evalsv(ie1,ikp)
           ! de_ = evalsv(ie2,ikp)-evalsv(ie1_,ikp)
           ! if (Gamma .and. (ie1/=ie2) .and. (ie1_/=ie2) .and. (abs(de)>1.d-3) .and. (abs(de_)>1.d-3) )  then     
           ! !isotropic 
           ! sx = sx + sxs2 * kiw(ie2,ik) * & 
              ! & zdotc(3,q0eps(:),1,pmatvv_k(ie2,ie1,:,ikp),1)/de * &
              ! & zdotc(3,pmatvv_k(ie2,ie1_,:,ikp),1,q0eps(:),1)/de_     
           ! end if 

          else
            !============================= 
            ! Core electron contribution
            !============================= 
            icg = ie2-nomax
            is = corind(icg,1)
            ia = corind(icg,2)
            ias = idxas(ia,is)
            ic = corind(icg,3)
            mvm = zdotc(mbsiz,minmmat(:,ie1,ie2),1,minmmat(:,ie1_,ie2),1)
            sx = sx-ciw(ic,ias)*mvm
          
          end if ! occupied states
        end do ! ie2
        
        ! add singular term (q->0)
        if (Gamma.and.(ie1==ie1_).and.(ie1<=nomax)) sx = sx+sxs2*singc2*kiw(ie1,ik)*kqset%nkpt
        !if (Gamma.and.(dabs(kiw(ie1,ik))>1.d-6)) sx = sx+sxs2*singc2*kiw(ie1,ik)*kqset%nkpt

        selfex_nn(ie1,ie1_,ikp) = selfex_nn(ie1,ie1_,ikp)+sx
        
       end do ! ie1_
      end do ! ie1
      
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
   

      
      
     ! debugging info
      if (input%gw%debug) then
        write(fdebug,*) 'EXCHANGE SELF-ENERGY: iq=', iq, ' ikp=', ikp
        write(fdebug,*) 'state   Sigma_x'
        do ie1_ = ibgw, nbgw
         do ie1 = ibgw, nbgw
          write(fdebug,*) ie1, selfex_nn(ie1,ie1_,ikp)
         end do
        end do
        write(fdebug,*)
      end if
      
    end do ! ikp
    end do ! ispn

    deallocate(minmmat)
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

    ! timing
    call timesec(tend)
    time_selfx = time_selfx+tend-tstart
    
    return

end subroutine


!--------------------------------------------------------------------------------
  subroutine calc_vxc_nn()
!
! This subroutine calculates the non-diagonal matrix elements of 
! the exchange correlation potential (only for valence states).

!!USES:
    use modinput
    use modmain, only : apwordmax, lmmaxapw, lmmaxvr, natmtot, nlomax, &
    &                   nstfv, nspinor, nstsv, nmatmax, nspecies, zzero, &
    &                   nmat, natoms, vxcmt, vxcir, zone
    use mod_vxc
    use modmpi
    use m_getunit

!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ik, ik0
    integer(4) :: ist, jst, i, j, k, l, ispn
    integer(4) :: ia, is
    integer(4) :: ngp, fid
    real(8)    :: tstart, tend
    complex(8) :: zsum, zt1
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfv(:,:), evecsv(:,:), evec(:,:,:)
    complex(8), allocatable :: h(:)
    real(8)    :: t0, t1

!!EXTERNAL ROUTINES: 
    complex(8), external :: zdotc

    call timesec(tstart)

    ! allocate exchange-correlation integral arrays
    if (allocated(vxcraa)) deallocate(vxcraa)
    allocate(vxcraa(apwordmax, &
    &               0:input%groundstate%lmaxmat, &
    &               apwordmax, &
    &               0:input%groundstate%lmaxapw, &
    &               0:lmmaxvr, &
    &               natmtot))
    if (allocated(vxcrloa)) deallocate(vxcrloa)
    allocate(vxcrloa(nlomax, &
    &                apwordmax, &
    &                0:input%groundstate%lmaxmat, &
    &                0:lmmaxvr, natmtot))
    if (allocated(vxcrlolo)) deallocate(vxcrlolo)
    allocate(vxcrlolo(nlomax, &
    &                 nlomax, &
    &                 0:lmmaxvr, &
    &                 natmtot))


    ! Calculate radial integrals
    call vxcrad
    
    ! Fourier transform the interstitial part of Vxc
    if (allocated(vxcig)) deallocate(vxcig)
    call genvxcig  

    allocate(apwalm(Gkset%ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evec(nmatmax,nstsv,nspinor))
    allocate(h(nmatmax))

	
    do ik = firstofset(rank,kset%nkpt), lastofset(rank,kset%nkpt)
    ! do ik = 1, kset%nkpt
        
      ik0 = kset%ikp2ik(ik)
      
      ngp = Gkset%ngk(1,ik0)

!#ifdef _HDF5_      
!      ! get the eigenvectors from file
!      write(cik,'(I4.4)') ik0
!      path = "/kpoints/"//trim(adjustl(cik))
!      call hdf5_read(fgwh5,path,"evecsv",evecsv(1,1),(/nmatmax,nstsv/))
!#else
      call getevecsvgw_new('GW_EVECSV.OUT',ik0,kqset%vkl(:,ik0),nmatmax,nstsv,nspinor,evec)
!#endif

      ! find the matching coefficients
      call match(ngp, &
      &          Gkset%gkc(:,1,ik0), &
      &          Gkset%tpgkc(:,:,1,ik0), &
      &          Gkset%sfacgk(:,:,1,ik0), &
      &          apwalm)

      do ispn = 1, nspinor
      
        
         do j = ibgw, nbgw
 
          h(:) = zzero
          ! muffin-tin contributions
          do is = 1, nspecies
          do ia = 1, natoms(is)
            call vxcaa(is,ia,ngp,apwalm,evec(:,j,ispn),h)
            call vxcalo(is,ia,ngp,apwalm,evec(:,j,ispn),h)
            call vxclolo(is,ia,ngp,evec(:,j,ispn),h)
          end do
          end do
          ! interstitial contribution
          call vxcistl(ngp,Gkset%igkig(:,1,ik0),evec(:,j,ispn),h)
          
          do i = ibgw, nbgw
           vxc_nn(i,j,ik) = vxc_nn(i,j,ik)+zdotc(nmat(1,ik),evec(:,i,ispn),1,h,1)
          end do ! i
          
         end do ! j 

      end do ! ispn
      
    end do ! ik
    
    
    deallocate(h)
    deallocate(apwalm)
    deallocate(evec)
    deallocate(vxcraa)
    deallocate(vxcrloa)
    deallocate(vxcrlolo)
    if (allocated(vxcig)) deallocate(vxcig)
    
#ifdef MPI
   call mpi_allgatherv_ifc(kset%nkpt,nbandsgw*nbandsgw,zbuf=vxc_nn)
   call barrier
#endif
! 
!     if (rank==0) then
! #ifdef _HDF5_
!       do ik = 1, kset%nkpt
!         write(cik,'(I4.4)') ik
!         path = "/kpoints/"//trim(adjustl(cik))
!         if (.not.hdf5_exist_group(fgwh5,"/kpoints",cik)) &
!         &  call hdf5_create_group(fgwh5,"/kpoints",cik)
!         call hdf5_write(fgwh5,path,"vxcnn",vxcnn(1,ik),(/nbandsgw/))
!       end do
! #else
!       call write_vxcnn()
! #endif
!     end if
      
    call timesec(tend)
    time_vxc = time_vxc + (tend-tstart)
    
  end subroutine

!--------------------------------------------------------------------------------
  subroutine calc_selfc_nn_wpol(ik,RE_EQP,IM_EQP,eta0,f_ac0)
    use m_getunit
    implicit none
    integer, intent(in) :: ik
    real(8), intent(in) :: RE_EQP,IM_EQP, eta0, f_ac0
    !local
    integer :: iq, fid
    integer :: i, nq(3)
        
    !eta = eta0 !input%gw%selfenergy%swidth
    etasign = eta0
    f_ac = f_ac0 !input%gw%selfenergy%f_ac
    IMEQP = IM_EQP
    REEQP = RE_EQP
    ! weights
    wkq = 1.d0 / dble(kqset%nkpt)
    ! reciprocal lattice vectors bvec(i,3)
    ! volume = bvec(1,:) * bvec(2,:) x bvec(3,:) or (2 pi)**3 / omega
    ! surface = 2 sum bvec(i,:) x bvec(j,:)
    nq(:)=input%gw%ngridq(:)
    wkq_eta = wkq / (2.d0*pi) * abs(dot_product(bvec(1,:),cross(bvec(2,:),bvec(3,:))))
    wkq_eta = wkq_eta / &
    &       ( norm2(cross(bvec(2,:),bvec(3,:)))/(nq(2)*nq(3)) + &
    &         norm2(cross(bvec(1,:),bvec(3,:)))/(nq(1)*nq(3)) + &
    &         norm2(cross(bvec(1,:),bvec(2,:)))/(nq(1)*nq(2)) )
        

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
        write(*,*) 'ERROR(mod_diag_selfe_wpol::task_calc_selfc_nn_wpol) Unknown singularity treatment scheme!'
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

    ! each process does a subset
    do iq = iqstart, iqend

      write(*,*)
      write(*,*) '(mod_non_diagonal_wpol::calc_selfc_nn_wpol) q-point cycle, iq = ', iq
    
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
      ! call calc_md_dmmd(iq)
      ! call diagonalize_dmmd(iq)
      ! call delete_coulomb_potential
      ! call clear_wpol()
      call get_wpol(iq)
      
      ! for eta/energy-derivative with respect to q
      call calc_q_surr(iq)
      if (allocated(tvck_s)) deallocate(tvck_s)
      allocate(tvck_s(3,nvck))
      call get_tvck(tvck_s(1,:),q_surr(1))
      call get_tvck(tvck_s(2,:),q_surr(2))
      call get_tvck(tvck_s(3,:),q_surr(3))
      

      ! Calculate q-dependent \Sigma^c_{nn'}(k,q,\omega)
      call calc_selfc_nn_wpol_q(ik,iq) !gets selfec_nn     

      ! clean unused data
      deallocate(tvck_s)
      deallocate(tvck)
      deallocate(wvck)
      deallocate(mpwipw)
      deallocate(barc)
      
    end do ! q-points
        

#ifdef MPI
    call mpi_sum_array(0,selfec_nn,nbandsgw,nbandsgw,mycomm_row)
#endif

    ! print self-energy to file
    !if (myrank==0) call write_selfecnn()

    ! delete index mapping arrays
    call del_wpol_indices()
!     call delete_freqgrid(freq)
!     call delete_k_vectors(kset)
!     call delete_G_vectors(Gset)
!     call delete_Gk_vectors(Gkset)
!     call delete_kq_vectors(kqset)
!     call delete_Gk_vectors(Gqset)
!     call delete_Gk_vectors(Gqbarc)

    return
  end subroutine


!--------------------------------------------------------------------------------
  subroutine calc_selfc_nn_wpol_q(ik,iq)
    implicit none
    integer, intent(in) :: iq, ik
    integer    :: ispn, jk, jkp, n, n_, i
    integer    :: d_jk(3), d_jkp(3)
    integer(8) :: recl
    real(8)    :: h(3)
    complex(8) :: zt1, zwt, coefs1, coefs2, sigma1, sigma2
    complex(8), allocatable :: mw(:,:),mw_(:,:)

    coefs1 = singc1 * sqrt(4.d0*pi/omega)
    coefs2 = singc2 * 4.d0*pi/omega

    allocate(mw(mdim,nvck))
    allocate(mw_(mdim,nvck))
    allocate(minmmat(1:mbsiz,ibgw:nbgw,1:mdim))

    ! loop over k-points -> is beeing done inside mod_diag_selfe_wpol
    do ispn = 1, nspinor

      write(*,*)
      write(*,'(a,3i8)') '(mod_non_diagonal::calc_selfc_nn_wpol_q): rank, (iq, ik):', myrank, iq, ik

      ! k-q point
      jk = kqset%kqid(ik,iq)
      jkp = kset%ik2ikp(jk)
                  
                  
      ! k-deltaq points for epsilon_mk derivative
      do i = 1, 3
        d_jk(i)  = kqset%kqid(ik,q_surr(i))
        d_jkp(i) = kset%ik2ikp(d_jk(i))
        ! h(i) = | q - q_i |
        h(i) = norm2( kqset%vkc(:,iq)-kqset%vkc(:,q_surr(i)) )
      end do
      
          
      ! Product basis coefficients M^i_{nm}(k,q)
      call calc_minmkq(ik,iq,ispn)

      ! loop over states 
      !additional loop over n for non diagonal elements!
      
      do n = ibgw, nbgw   
        
        ! precalculate \sum_{i} conjg(M^i_{nm})*w^i_{vck} in Eq.(28) in paper -> mw(n,s)
        call zgemm( 'c', 'n', mdim, nvck, mbsiz, &
        &           zone, minmmat(1:mbsiz,n,1:mdim), mbsiz, &
        &           wvck(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw, mdim)
      
      do n_ = ibgw, nbgw

        !for other n_  -> mw_
        call zgemm( 'c', 'n', mdim, nvck, mbsiz, &
        &           zone, minmmat(1:mbsiz,n_,1:mdim), mbsiz, &
        &           wvck(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw_, mdim)

        !calculate sigma matrix
        zt1 = 0.d0
        ! sum over states
        zt1 = zt1 + sum_valence(jkp,mw,mw_,d_jkp,h)
        zt1 = zt1 + sum_core(jkp,mw,mw_)

        
        selfec_nn(n,n_) = selfec_nn(n,n_) + wkq*zt1

        !treat the singularity at GAMMA 
        if (Gamma) then
         call sum_singular(n,n_,jkp,mw,mw_,sigma1,sigma2,d_jkp,h)
         selfec_nn(n,n_) = selfec_nn(n,n_) + coefs1*sigma1 + coefs2*sigma2
        end if

      end do ! n_     
      end do ! n

    end do ! nspinor
     
    deallocate(minmmat)
    deallocate(mw,mw_)

   return
  end subroutine
  
  
!--------------------------------------------------------------------------------
   subroutine sum_singular(n,n_,jkp,mw,mw_,sigma1,sigma2,d_jkp,h)
    implicit none
    integer, intent(in) :: n,n_, jkp, d_jkp(3)
    real(8), intent(in) :: h(3)
    complex(8), intent(in) :: mw(mdim,nvck), mw_(mdim,nvck)
    complex(8)             :: sigma1, sigma2
    ! local
    integer    :: i, j
    real(8)    :: enk,en_k , d_enk(3), d_en_k(3), d_tvck(3), d_tvck_(3)
    complex(8) :: zt1, zt2, p(3), zt1_, pq, zt3, zt3_
    real(8)    :: q0eps(3), modq0
    real(8)    :: c1, c2, de, eta_

	complex(8), external :: zdotc
	
    enk  = evalsv(n,jkp)
    en_k = evalsv(n_,jkp)

    ! d/dq enk
    do j = 1, 3
     d_enk(j)  = enk  - evalsv(n,d_jkp(j))
     d_en_k(j) = en_k - evalsv(n_,d_jkp(j))
    end do

    sigma1 = zzero
    sigma2 = zzero

    select case(input%gw%scrcoul%sciavtype)
    
    case('isotropic')

      ! q->0 direction
      q0eps(:) = input%gw%scrcoul%q0eps(:)
      modq0    = sqrt(q0eps(1)**2+q0eps(2)**2+q0eps(3)**2)
      q0eps(:) = q0eps(:)/modq0

	  !p_n'nk * vec(q)/q
	  pq = zdotc(3,q0eps(:),1,pmatvv_k(n_,n,:,jkp),1) 
	  de = enk-en_k
	  
#ifdef USEOMP
!$omp parallel &
!$omp default(shared) &
!$omp private(i,j,zt1,zt1_,zt2,zt3,zt3_,d_tvck,eta,d_tvck_,eta_)
!$omp do reduction (+:sigma1,sigma2)
#endif
      do i = 1, nvck
        
        !d/dq tvck
        do j = 1, 3
          d_tvck(j) = sign(1,nomax-n)*(tvck(i)-tvck_s(j,i))
          d_tvck_(j) = sign(1,nomax-n_)*(tvck(i)-tvck_s(j,i))
        end do
        
        eta  = etasign * wkq_eta * norm2( (d_enk-d_tvck)/h )
        eta_ = etasign * wkq_eta * norm2( (d_en_k-d_tvck_)/h )

        !---------------------------------------------------
        zt1 = tvck(i) * ( REEQP + zi*IMEQP - enk + sign(1,nomax-n)*tvck(i) )
        zt1 = 0.5d0 / zt1      
        zt1_ = tvck(i) * ( REEQP + zi*IMEQP - en_k  + sign(1,nomax-n_)*tvck(i) )
        zt1_ = 0.5d0 / zt1_ 

		zt3 = tvck(i) * ( REEQP + zi*eta - enk + sign(1,nomax-n)*tvck(i) )
        zt3 = 0.5d0 / zt3
		zt3_ = tvck(i) * ( REEQP + zi*eta_ - en_k + sign(1,nomax-n_)*tvck(i) )
        zt3_ = 0.5d0 / zt3_
        
        zt2 = wvck(mbsiz+1,i)
        
        
        !---------------------------------------------------
        ! contribution from the first term: 1/q^2
        if (n==n_) then
        sigma2 = sigma2 + zt1*zt2*conjg(zt2) + 2*f_ac*zi*AIMAG(zt3*zt2*conjg(zt2))
        end if
        !---------------------------------------------------
        ! contributions from the second+third term: 1/q
        !for all n,n_:
        sigma1 = sigma1 + zt1_*mw(n_,i)*conjg(zt2) + zt1*zt2*conjg(mw_(n,i))  + 2*f_ac*zi*AIMAG(zt3_*mw(n_,i)*conjg(zt2)) + 2*f_ac*zi*AIMAG(zt3*zt2*conjg(mw_(n,i)))
        !non diagonal contribution 
        if ((n/=n_) .and. (de>1.d-3)) then
        sigma1 = sigma1 - sqrt(4.d0*pi/omega) * zt1*zt2*conjg(zt2) * pq / de     - 2*f_ac*zi*AIMAG(sqrt(4.d0*pi/omega) * zt3*zt2*conjg(zt2) * pq / de)
        sigma1 = sigma1 - sqrt(4.d0*pi/omega) * zt1_*zt2*conjg(zt2) * pq / (-de) - 2*f_ac*zi*AIMAG(sqrt(4.d0*pi/omega) * zt3_*zt2*conjg(zt2) * pq / (-de)) !zdotc(3,pmatvv_k(n,n_,:,jkp),1,q0eps(:),1)/(en_k-enk)
        end if
        
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
  function sum_valence(jkp,mw,mw_,d_jkp,h) result(zsum)
    implicit none
    integer,    intent(in) :: jkp, d_jkp(3)
    real(8), intent(in) :: h(3)
    complex(8), intent(in) :: mw(mdim,nvck),mw_(mdim,nvck)
    complex(8)             :: zsum
    ! local
    integer    :: m, i, j
    real(8)    :: enk, d_enk(3), d_tvck(3)
    complex(8) :: zt1,zt2
    complex(8), allocatable :: mwt(:),mwt2(:)
    complex(8), external :: zdotc

    zsum = 0.d0
    allocate(mwt(nvck))
	allocate(mwt2(nvck))
    mwt(:) = 0.d0
    mwt2(:) = 0.d0
	
    ! sum over states
    do m = 1, nstse
      enk = evalsv(m,jkp)
      ! d/dq enk
      do j=1,3
       d_enk(j) = enk - evalsv(m,d_jkp(j))
      end do
#ifdef USEOMP
!$omp parallel &
!$omp default(shared) &
!$omp private(i,j,zt1,zt2,d_tvck,eta)
!$omp do
#endif      
      ! apply frequency/state dependent prefactor
      do i = 1, nvck 
      
        ! d/dq tvck
        do j = 1, 3
         d_tvck(j) = sign(1,nomax-m)*(tvck(i)-tvck_s(j,i))
        end do        
        ! eta = wkq^(1/3) 1/pi V_BZ/O_BZ * d/dq eps+-tau
        eta = etasign * wkq_eta * norm2( (d_enk-d_tvck) / h )

		zt1 = tvck(i) * ( REEQP + zi*IMEQP - enk + sign(1,nomax-m)*tvck(i) )
        zt1 = 0.5d0 / zt1
        mwt(i) = zt1*mw(m,i)
		
		!zt2 = tvck(i) * ( REEQP  - enk + sign(1,nomax-m)*(tvck(i)-zi*eta) )
		zt2 = tvck(i) * ( REEQP + zi*eta - enk + sign(1,nomax-m)*tvck(i) )
        zt2 = 0.5d0 / zt2
		mwt2(i) = zt2*mw(m,i)
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif      
      ! sum over vck
      zsum = zsum + zdotc(nvck,mw_(m,:),1,mwt,1) + 2*f_ac*zi*AIMAG(zdotc(nvck,mw_(m,:),1,mwt2,1))

    end do ! m
    deallocate(mwt,mwt2)

    return
  end function

!--------------------------------------------------------------------------------
  function sum_core(jkp,mw,mw_) result(zsum)
    implicit none
    integer,    intent(in) :: jkp
    complex(8), intent(in) :: mw(mdim,nvck),mw_(mdim,nvck)
    complex(8)             :: zsum
    ! local
    integer    :: m, i, j
    integer    :: icg, is, ia, ias, ic
    real(8)    :: enk
    real(8)    :: x, y, t1, t2
    complex(8) :: zt1, zt2
    complex(8), allocatable :: mwt(:),  mwt2(:)
    complex(8), external    :: zdotc, zdotu

    zsum = 0.d0
    allocate(mwt(nvck))
	allocate(mwt2(nvck))
    ! sum over states
    do m = nstse+1, mdim
      icg = m-nstse
      is  = corind(icg,1)
      ia  = corind(icg,2)
      ias = idxas(ia,is)
      ic  = corind(icg,3)
      enk = evalcr(ic,ias)
#ifdef USEOMP
!$omp parallel &
!$omp default(shared) &
!$omp private(i,zt1,zt2,eta)
!$omp do
#endif
      ! apply frequency/state dependent prefactor
      do i = 1, nvck
        
        !no core state summation for Si..
        eta = etasign * input%gw%selfenergy%swidth
      
        zt1 = tvck(i) * ( REEQP + zi*IMEQP - enk + tvck(i) )
        zt1 = 0.5d0 / zt1
        mwt(i) = zt1*mw(m,i)
				
		zt2 = tvck(i) * ( REEQP + zi*eta - enk + tvck(i) )
        zt2 = 0.5d0 / zt2
		mwt2(i) = zt2*mw(m,i)
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
      ! sum over vck
      zsum = zsum + zdotc(nvck,mw_(m,:),1,mwt,1) + 2*f_ac*zi*AIMAG(zdotc(nvck,mw_(m,:),1,mwt2,1))
    end do ! m
    deallocate(mwt,mwt2)

    return
  end function  

end module

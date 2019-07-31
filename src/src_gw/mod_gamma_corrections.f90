 module mod_gamma_correction
!================================
! calculates self-energy correction GWGAMMA^1 on frequency grid
! usage:
! input.xml 
!  -> add <task name="gwgamma_correction"/> to gwplan
!  -> set up freqgrid   
! 
! Output will be seperated into selfx_GvvGG, selfx_GvWGG, selfx_GWvGG, selfc_GWWGG
!
!================================
  use modinput
  use modgw
  use mod_mpi_gw  
  use mod_selfenergy, only : read_evalqp, evalks, evalqp, eferqp,singc1, singc2
  use mod_non_diagonal_wpol, only : calc_selfx_nn, calc_vxc_nn , calc_selfc_nn_wpol ,vxc_nn , selfex_nn, selfec_nn
  use mod_wpol_diagonalization, only : mkl_zheev, mkl_zheevr
  use modmain
  use mod_wpol !, only: tvck,wvck,get_wpol,nvck,mbdim !,del_wpol_indices,set_wpol_indices
  
  use modfvsystem

  implicit none
  !integer,    public :: mdim, nvdim

  real(8), public :: eta, wkq 
  
  complex(8), public :: pq_q
  
  complex(8), allocatable, public ::  selfex_GvvGG_nn(:,:),selfex_GvWGG_nn(:,:), selfex_GWvGG_nn(:,:) , selfec_GWWGG_nn(:,:) !,  selfex_GvvGG_nn_k(:,:,:), selfex_GvWGG_nn_k(:,:,:), selfex_GWvGG_nn_k(:,:,:)

  public :: calc_selfx_GvvGG_wpol_nn, calc_selfx_GvWGG_wpol_nn, calc_minmkq_gc, gc_q0_treatment_setup

contains

  subroutine task_gamma_correction()
    implicit none 

    integer :: ikp,ib,iq,jk,ik,jkp
    
    !Sigma = GW + GWWGG

    call init_gw   

    !================================
    !loop over irreducible k-points
    !================================
    
   do ikp = 1, kset%nkpt
    
    !if (ikp /= 1) cycle
      
    eta  = input%gw%selfenergy%swidth
  
    !(1) Sigma^(1)_exchange GvvGG
    allocate(selfex_GvvGG_nn(ibgw:nbgw,freq%nomeg))
    selfex_GvvGG_nn(:,:) = zzero
    call calc_selfx_GvvGG_wpol_nn(ikp)
    if (myrank==0) call write_selfecnn_gc(ikp,selfex_GvvGG_nn(:,:),'GvvGG')
    deallocate(selfex_GvvGG_nn)
!     
! 
    !(2) Sigma^(1)_exchange GvW^cGG
      allocate(selfex_GvWGG_nn(ibgw:nbgw,freq%nomeg))
      selfex_GvWGG_nn(:,:) = zzero
      call calc_selfx_GvWGG_wpol_nn(ikp)
      if (myrank==0)   call write_selfecnn_gc(ikp,selfex_GvWGG_nn(:,:),'GvWGG')
      deallocate(selfex_GvWGG_nn)
! 			   
! 			   
! !     !!(3) Sigma^(1)_exchange GW^cvGG
      allocate(selfex_GWvGG_nn(ibgw:nbgw,freq%nomeg))
      selfex_GWvGG_nn(:,:) = zzero
      call calc_selfx_GWvGG_wpol_nn(ikp)
      if (myrank==0)  call write_selfecnn_gc(ikp,selfex_GWvGG_nn(:,:),'GWvGG')
      deallocate(selfex_GWvGG_nn) 
  
!     !! Sigma^(1)_correclation
      allocate(selfec_GWWGG_nn(ibgw:nbgw,freq%nomeg))
      selfec_GWWGG_nn(:,:) = zzero
      call calc_selfc_GWWGG_wpol_nn(ikp)
      if (myrank==0)   call write_selfecnn_gc(ikp,selfec_GWWGG_nn(:,:),'GWWGG')
      deallocate(selfec_GWWGG_nn)
        
   
  end do !ikp
  
  end subroutine
  
!   !--------------------------------------------------------------------------------
  subroutine calc_selfx_GvvGG_wpol_nn(ikp)
 
  use m_getunit
            
  implicit none
  !real(8), intent(in) :: E_n_QP
  integer(4), intent(in) :: ikp
  
  integer(4) ::  iq, fid, ik_,iq2,iqmax
  integer(4) :: recl
  real(8) :: wk

  !===========================================================================
  ! Initialization
  !===========================================================================
  !call init_gw
  !call clean_gndstate
    
  ! occupancy dependent BZ integration weights
  !call kintw
   wkq = 1.d0 / dble(kqset%nkpt)
   

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
  
  !if (myrank==0) write(*,*) kset%nkpt  , kqset%nkpt  , wkq

   
   !===========================================================================
  ! Main loop: BZ integration
  !===========================================================================   
  

iqmax = INT(kqset%nkpt*kqset%nkpt)

#ifdef MPI
  call set_mpi_group(iqmax)
  call mpi_set_range(nproc_row, &
  &                  myrank_row, &
  &                  iqmax, 1, &
  &                  iqstart, iqend)
#else
  iqstart = 1
  iqend = iqmax
#endif

  ! each process does a subset
  !do iq = iqstart, iqend
  do iq2 = iqstart, iqend
  
  
    iq = mod(iq2,kqset%nkpt)
    if (iq==0) iq=kqset%nkpt

    ik_ = 1 + INT((iq2-1)/kqset%nkpt)
  
    write(*,*)
    write(*,*) '(mod_gamma_corrections/GvvGG_wpol): q-point cycle, iq = ', iq2,iq,ik_

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

	call setbarcev(input%gw%barecoul%barcevtol)
	 
    !===============================
    ! Calculate \Sigma^{x}_{kn}(q)
    !===============================
    
    call calc_selfx_GvvGG_wpol_nn_q(ikp,ik_,iq)  
  
    ! clean unused data
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
      
  end do ! iq

#ifdef MPI
    call mpi_sum_array(0,selfex_GvvGG_nn,nbandsgw,freq%nomeg,mycomm_row)
	!call mpi_sum_array(0,selfex_GvvGG_nn_k,nbandsgw,freq%nomeg,kqset%nkpt,mycomm_row) 
#endif
  
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


subroutine calc_selfx_GvWGG_wpol_nn(ikp)
use m_getunit
            
  implicit none
  integer(4), intent(in) :: ikp

  integer(4) ::  iq, fid,ik_,iq2,iqmax
  integer(4) :: recl

  !===========================================================================
  ! Main loop: BZ integration
  !===========================================================================    
   !weights
   wkq = 1.d0 / dble(kqset%nkpt)

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

    !===========================
    ! Momentum matrix elements
    !===========================
    if (.not.input%gw%rpmat) call calcpmatgw

    if (myrank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if

 ! allocate(selfex_GvWGG_nn_k(ibgw:nbgw,freq%nomeg,kqset%nkpt))


iqmax = INT(kqset%nkpt*kqset%nkpt)

#ifdef MPI
  call set_mpi_group(iqmax)
  call mpi_set_range(nproc_row, &
  &                  myrank_row, &
  &                  iqmax, 1, &
  &                  iqstart, iqend)
#else
  iqstart = 1
  iqend = iqmax
#endif

 call set_wpol_indices()

  ! each process does a subset
  !do iq = iqstart, iqend
  do iq2 = iqstart, iqend
  
  
    iq = mod(iq2,kqset%nkpt)
    if (iq==0) iq=kqset%nkpt

    ik_ = 1 + INT((iq2-1)/kqset%nkpt)
  
    write(*,*)
    write(*,*) '(mod_gamma_corrections/GvWGG_wpol): q-point cycle, iq = ', iq2,iq,ik_  
  
  
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
     
    call setbarcev(input%gw%barecoul%barcevtol) 
     
    call get_wpol(iq)
    
    call calc_selfx_GvWGG_wpol_nn_q(ikp,ik_,iq)
    
  
    ! clean unused data
    deallocate(tvck)
    deallocate(wvck)
    
        
    ! clean unused data
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
      
  end do ! iq2

  
#ifdef MPI
    call mpi_sum_array(0,selfex_GvWGG_nn,nbandsgw,freq%nomeg,mycomm_row)
#endif

  return
end subroutine



subroutine calc_selfx_GWvGG_wpol_nn(ikp)
use m_getunit
            
  implicit none
  integer(4), intent(in) :: ikp

  integer(4) ::  iq, fid,ik_,iq2,iqmax
  integer(4) :: recl

  !===========================================================================
  ! Main loop: BZ integration
  !===========================================================================    
   !weights
   wkq = 1.d0 / dble(kqset%nkpt)

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

    !===========================
    ! Momentum matrix elements
    !===========================
    if (.not.input%gw%rpmat) call calcpmatgw

    if (myrank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if

iqmax = INT(kqset%nkpt*kqset%nkpt)


#ifdef MPI
  call set_mpi_group(iqmax)
  call mpi_set_range(nproc_row, &
  &                  myrank_row, &
  &                  iqmax, 1, &
  &                  iqstart, iqend)
#else
  iqstart = 1
  iqend = iqmax
#endif

  call set_wpol_indices()

  do iq2 = iqstart, iqend
  
  
    iq = mod(iq2,kqset%nkpt)
    if (iq==0) iq=kqset%nkpt

    ik_ = 1 + INT((iq2-1)/kqset%nkpt)
  
    write(*,*)
    write(*,*) '(mod_gamma_corrections/GWvGG_wpol): q-point cycle, iq = ', iq2,iq,ik_  
    
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
     
    call setbarcev(input%gw%barecoul%barcevtol) 
         
    call calc_selfx_GWvGG_wpol_nn_q(ikp,ik_,iq)
    
    ! clean unused data
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
      
  end do ! iq
  
#ifdef MPI
    call mpi_sum_array(0,selfex_GWvGG_nn,nbandsgw,freq%nomeg,mycomm_row)
	!call mpi_sum_array(0,selfex_GWvGG_nn_k,nbandsgw,freq%nomeg,kqset%nkpt,mycomm_row) 
#endif

   
  return
end subroutine


subroutine calc_selfc_GWWGG_wpol_nn(ikp)
use m_getunit
            
  implicit none
  integer(4), intent(in) :: ikp

  integer(4) ::  iq, fid,ik_,iq2,iqmax
  integer(4) :: recl

  !===========================================================================
  ! Main loop: BZ integration
  !===========================================================================    
   !weights
   wkq = 1.d0 / dble(kqset%nkpt)

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

    !===========================
    ! Momentum matrix elements
    !===========================
    if (.not.input%gw%rpmat) call calcpmatgw

    if (myrank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if

iqmax = INT(kqset%nkpt*kqset%nkpt)


#ifdef MPI
  call set_mpi_group(iqmax)
  call mpi_set_range(nproc_row, &
  &                  myrank_row, &
  &                  iqmax, 1, &
  &                  iqstart, iqend)
#else
  iqstart = 1
  iqend = iqmax
#endif

  call set_wpol_indices()

  do iq2 = iqstart, iqend
  
  
    iq = mod(iq2,kqset%nkpt)
    if (iq==0) iq=kqset%nkpt

    ik_ = 1 + INT((iq2-1)/kqset%nkpt)
  
    write(*,*)
    write(*,*) '(mod_gamma_corrections/GWWGG_wpol): q-point cycle, iq = ', iq2,iq,ik_  
    
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
     
    call setbarcev(input%gw%barecoul%barcevtol) 
     
    call get_wpol(iq)
    
    call calc_selfc_GWWGG_wpol_nn_q(ikp,ik_,iq)
    
        
    ! clean unused data
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
      
  end do ! iq
  
#ifdef MPI
    call mpi_sum_array(0,selfec_GWWGG_nn,nbandsgw,freq%nomeg,mycomm_row)
#endif

  return
end subroutine



  subroutine gc_q0_treatment_setup(prefac,p)
    use modinput
    use modmain, only : pi, zi
    implicit none
    ! inpout
    real(8),    intent(in) :: prefac
    complex(8), intent(in) :: p(3)
    ! local
    real(8) :: c1, c2
    real(8) :: q0eps(3), modq0
    
    select case(input%gw%scrcoul%sciavtype)
        case('isotropic')

      ! q->0 direction
      q0eps(:) = input%gw%scrcoul%q0eps(:)
      modq0    = sqrt(q0eps(1)**2+q0eps(2)**2+q0eps(3)**2)
      q0eps(:) = q0eps(:)/modq0

      pq_q = prefac * ( p(1)*q0eps(1)+ &
      &                           p(2)*q0eps(2)+ &
      &                           p(3)*q0eps(3) )

        case default
      write(*,*) "ERROR(mod_wpol::calc_md_dmmd): Unknown averaging type!"
        stop

    end select
  end subroutine

   
  subroutine calc_minmkq_gc(ik,iq,ispn,minmmat,mf)
    implicit none
    integer, intent(in) :: ik, iq, ispn
    
    character(2), intent(in) :: mf
    
    complex(8), intent(inout) :: minmmat(:,:,:)
    
    integer :: jk, im
    complex(8), allocatable :: evecsv(:,:,:)

	!ik=kset%ikp2ik(ikp)
	
	minmmat(:,:,:) = zzero

    ! k-q point
    jk = kqset%kqid(ik,iq)

    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))    

   !M1
   if (mf=="nv") then 
   
	!Psi_n1_k^* Psi_n2_k-q chi_q^*
    allocate(evecsv(nmatmax,nstsv,nspinor))
    call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
    eveckp = evecsv(:,:,ispn)
    call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
    eveck = conjg(evecsv(:,:,ispn))
    deallocate(evecsv)
    call expand_evec(ik,'l')   !eveck
    call expand_evec(jk,'r')  !eveckp
    call expand_products(ik,iq,ibgw,nbgw,-1,1,nstse,nstse,minmmat)
     
     
   else if (mf=="vv") then    !else v,c
      
    allocate(evecsv(nmatmax,nstsv,nspinor))
    call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
    eveckp = conjg(evecsv(:,:,ispn))
    call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
    eveck = evecsv(:,:,ispn)
    deallocate(evecsv)
    call expand_evec(ik,'t')
    call expand_evec(jk,'c')
    call expand_products(ik,iq,1,nstse,nstse,1,nstse,-1,minmmat)
 
   end if 

     
    deallocate(eveckalm)
    deallocate(eveckpalm)
    deallocate(eveck)
    deallocate(eveckp)

    return
  end subroutine  

    
 
    subroutine write_selfecnn_gc(ikp,selfc,cha)
        use modgw, only : kset, freq, ibgw, nbgw
      use m_getunit
      implicit none
      
       complex(8) , intent(in) :: selfc(ibgw:nbgw,1:freq%nomeg)
      integer, intent(in) :: ikp
       character(5), intent(in) :: cha
      
      
      integer :: iom, n, ik
      integer :: fid1, fid2
      character(80) :: frmt, fname1, fname2


      ! text format
      n = nbgw-ibgw+1
      write(frmt,'("(",i8,"f14.6)")') 1+n
      ! write(*,*) trim(frmt)
      !write(frmt,'("(",i8,"f14.6)")') 1+1

      call getunit(fid1)
      write(fname1,'("selfx-RE_",A5,"-k",I4.4,".DAT")') cha,ikp
      open(fid1, File=trim(fname1), Action='WRITE')
      call getunit(fid2)
      write(fname2,'("selfx-IM_",A5,"-k",I4.4,".DAT")') cha,ikp
      open(fid2, File=trim(fname2), Action='WRITE')
      !do ik = 1, kset%nkpt
        !write(fid1,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
        !write(fid2,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
        do iom = 1, freq%nomeg
          write(fid1,trim(frmt)) freq%freqs(iom), dble(selfc(:,iom))
          write(fid2,trim(frmt)) freq%freqs(iom), imag(selfc(:,iom))
        end do
        !write(fid1,*); write(fid1,*)
        !write(fid2,*); write(fid2,*)
      !end do
      close(fid1)
      close(fid2)

      ! binary format
!       write(fname1,'("selfx-",A5,"-k",I4.4,".OUT")') cha,ikp
!       open(fid1, File=trim(fname1),form='UNFORMATTED',status='UNKNOWN')
!       write(fid1) ibgw, nbgw, freq%nomeg, kset%nkpt, selfc(ibgw:nbgw,1:freq%nomeg)
!       close(fid1)

      return
    end subroutine
    
    

end module  

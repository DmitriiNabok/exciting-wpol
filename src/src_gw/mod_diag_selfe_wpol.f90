module mod_diag_selfe_wpol

!--------------------------------------------------------------------------------
! NOTES ON input.xml 
!
! requires a wpol calculation, marked in the comment block <!-- --> below 
! actual self-energy Diagonalization is invoked by <task name="selfe_diag"/>
!
!    <gwplan>
!      <!-- <task name="vxc"/>
!       <task name="selfx"/>
!       <task name="wpol"/>
!       <task name="selfc_wpol"/>
!       <task name="evalqp_wpol"/>-->
!     <task name="selfe_diag"/>
!    </gwplan>
!
! Additonally requries a parameter f_ac, which controls the Riemann sheet; f_ac=1 corresponds to second sheet
!     <selfenergy 
!       singularity="mpb"
!       swidth="1.d-4"
!       f_ac="1.d0"
!      >
!
! Also check parameter nbgw_diag in code below, which controls the output number of states
!
! Results are written to eval-kXXX.OUT; order of complex energies is: minimization input -> output -> whole grid..
!
!--------------------------------------------------------------------------------


  use modinput
  use modgw
  use mod_selfenergy, only : read_evalqp, evalks, evalqp, eferqp
  use mod_non_diagonal_wpol
  use mod_wpol_diagonalization, only : mkl_zheev, mkl_zheevr
  use modmain
  use mod_wpol !, only: tvck,wvck,get_wpol,nvck,mbdim !,del_wpol_indices,set_wpol_indices
  
  use modfvsystem

  implicit none
  
  complex(8), allocatable :: evaldiag(:)
  complex(8), allocatable :: dSigma(:,:), zevec_left(:,:)
  complex(8), allocatable :: evec_qp_right(:,:,:,:),evec_qp_left(:,:,:,:),eval_qp_diag(:,:)
  complex(8), allocatable :: G_nmk(:,:,:)
  complex(8), allocatable :: evec_ks(:,:,:),evec_overlap(:,:) 
  

  integer, public  :: nbgw_diag , nmatp, ispn
  
  real(8), public :: omega_pl,kappa
  
  public :: task_diag_selfe_wpol

contains

!--------------------------------------------------------------------------------
  subroutine task_diag_selfe_wpol()
    implicit none 
    
    type(evsystem) :: system
    integer :: ik, ib, ik0
    integer :: i, j,  n, m,  n_qp, n_sat
    real(8) :: eval_gamma, eta0, gamma0, omega_sat
     !, evecfv(:,:,:)
    complex(8), allocatable :: apwalm_s(:,:,:,:,:)
	complex(8), external :: zdotc
	    
    call init_gw   
	
	!set number of states for output; usually nbgw, but needs a lot of computing time
	nbgw_diag = nbgw
	
    allocate(vxc_nn(ibgw:nbgw,ibgw:nbgw,kset%nkpt))
    vxc_nn(:,:,:) = zzero
    !calculate Vxc_nn'k -> vxc_nn(n,n',k)
    call calc_vxc_nn()
    	
	allocate(selfex_nn(ibgw:nbgw,ibgw:nbgw,kset%nkpt))
    selfex_nn(:,:,:) = zzero
    !calculate Sigma^x_nn'k -> selfex_nn(n,n',k)
    call calc_selfx_nn()

	! KS eigenvalues 
    if (allocated(evalks)) deallocate(evalks)
    allocate(evalks(ibgw:nbgw,kset%nkpt))
    evalks(:,:) = 0.d0
    ! Quasi-Particle energies (linearized solver)
    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(ibgw:nbgw,kset%nkpt))
    evalqp(:,:) = 0.d0
	!get initial QP energies -> evalqp(n,k)
    call read_evalqp()

	!momentum matrix: p_nn'k
	allocate(pmatvv_k(ibgw:nbgw,ibgw:nbgw,3,kset%nkpt))
	if ((input%gw%coreflag=='all') .or. (input%gw%coreflag=='xal')) allocate(pmatcv_k(1:ncg,ibgw:nbgw,3,kset%nkpt))
	pmatvv_k(:,:,:,:) = zzero
	pmatcv_k(:,:,:,:) = zzero
	!call calcpmat_non_diag()
        
    !new QP coefficients (nstsv->nbandsgw)
    allocate(evec_qp_right(nmatmax,nbgw_diag,2,1)) !last 1 = nspinor
    allocate(evec_qp_left(nbgw_diag,nmatmax,2,1))
    allocate(eval_qp_diag(nbgw_diag,7))
    eval_qp_diag(:,:)=zzero
	
	!can sum over spins, but for plotting G_ispn=1 is suffiecient if there is no spin polarization
    ispn = 1
    
    
	 !satellites, plasma freq. for nsat
 	
    

    !loop over all E_mk^QP -> for each diagonalize full self-energy matrix
    !----------------------------------------------------------------------------------------------------------------------------------------------------
     do ik = 1, kset%nkpt !Sigma_nk(E_mk)
!     
!     !only Gamma=1, L=3 , X=7 for 4x4x4 grid
!     !if ((ik/=1) .and. (ik/=3) .and. (ik/=7)) cycle
!     !GAMMA
!     !if (ik/=1) cycle
!     !L
!     !if (ik/=3) cycle
!     !X
!     !if (ik/=7) cycle
! 
     ik0 = kset%ikp2ik(ik)
    
    ! call get_satellites(ik0) 
     

      !KS coefficients
    allocate(evec_ks(nmatmax,nstsv,nspinor)) 

     call getevecsvgw_new('GW_EVECSV.OUT',ik0,kqset%vkl(:,ik0),nmatmax,nstsv,nspinor,evec_ks)
     !overlap matrix S
 	 nmatp = nmat(1,ik)
 	 allocate(apwalm_s(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
	 call newsystem(system,input%groundstate%solver%packedmatrixstorage,nmatp)
	 call hamiltonandoverlapsetup(system,ngk(1,ik),apwalm_s, &
     &                            igkig(:,1,ik),vgkc(:,:,1,ik))
     !conjg(c)*S
      allocate(evec_overlap(nstsv,nmatp)) 
     call zgemm('c','n',nstsv,nmatp,nmatp, &
     &          zone,evec_ks(1:nmatp,1:nstsv,ispn),nmatp, &
     &          system%overlap%za,nmatp, &
     &          zzero,evec_overlap(1:nstsv,1:nmatp),nstsv)
	 deallocate(apwalm_s)	
	 



      !start ib -----------------------------------------------------------------------------------------------------------------------------------------
     do ib = ibgw, nbgw_diag
      
     !new QP energies ->  eval_qp_diag(ib,2)
     !gamma->0 where E<E_F-E_g
     !Gamma~(E_QP-E_F)^2
     
     allocate(selfec_nn(ibgw:nbgw,ibgw:nbgw))
     allocate(dSigma(ibgw:nbgw,ibgw:nbgw))
     allocate(zevec_left(ibgw:nbgw,ibgw:nbgw))
     allocate(evaldiag(nbandsgw))
     
     !Fermi-liquid quadratic behaviour for Im E^QP
     gamma0 = 1.d0/26.2d0 * (evalqp(ib,ik)-eferqp)**2 / eferqp   !for silicon at -12 eV @ GAMMA ~ 1 eV  => empirical factor 1/26.2
 
     !sign of imag. part for occupied/unocc. respectively
	 if (sign(1,nomax-ib) == 1) then
      eval_gamma = gamma0
      eta0  = - 1.d0 !input%gw%selfenergy%swidth
	 else
	  eval_gamma = - gamma0
	  eta0  = 1.d0 !input%gw%selfenergy%swidth
	 end if
     
     kappa=1.d0
     
     !input QP energy
     eval_qp_diag(ib,1) = evalqp(ib,ik) + zi*eval_gamma 
     
     !energy grid over AIMAG, get min of | E^QP - E^QP^0 | 
     call get_eval_qp_diag(ib,ik,evalqp(ib,ik),0.50*eval_gamma,eta0,input%gw%selfenergy%f_ac,n_qp,3)
     call get_eval_qp_diag(ib,ik,evalqp(ib,ik),0.75*eval_gamma,eta0,input%gw%selfenergy%f_ac,n_qp,4)
     call get_eval_qp_diag(ib,ik,evalqp(ib,ik),     eval_gamma,eta0,input%gw%selfenergy%f_ac,n_qp,5)
     call get_eval_qp_diag(ib,ik,evalqp(ib,ik),1.25*eval_gamma,eta0,input%gw%selfenergy%f_ac,n_qp,6)
     call get_eval_qp_diag(ib,ik,evalqp(ib,ik),1.50*eval_gamma,eta0,input%gw%selfenergy%f_ac,n_qp,7)


! 	 !-------------------------------------------------------------------------------------------------------------------------------------------------------------
!      !satellites for valence states
!      !-------------------------------------------------------------------------------------------------------------------------------------------------------------
!      ! sign(1,nomax-n_) = 1 for occupied
!      if (sign(1,nomax-ib) == 1) then
!       !shift qp-energy by plasma frequency to obtain replica
!       omega_sat = evalqp(ib,ik)-1.3d0*omega_pl
!       eval_gamma = 1.d0/26.d0 * (omega_sat-eferqp)**2 / eferqp
!       eta0  = - 1.d0
! 
! 	  eval_qp_diag(ib,1) = omega_pl + zi*0.d0
!       eval_qp_diag(ib,2) = omega_sat + zi*eval_gamma
! 	  
!      !get the satellite values
!      call get_eval_qp_diag(ib,ik,omega_sat,eval_gamma,eta0,input%gw%selfenergy%f_ac,n_qp,3)
! 
!      !left EV: conjg(VL)*conjg(c_n)*S
!      call zgemm('c','n',1,nmatp,nbandsgw, &
!      &          zone,zevec_left(ibgw:nbgw,n_qp),nbandsgw, &
!      &          evec_overlap(ibgw:nbgw,1:nmatp),nbandsgw, &
!      &          zzero,evec_qp_left(ib,1:nmatp,2,ispn),1)
! 	 !right EV: c_n*VR
!      call zgemm('n','n',nmatp,1,nbandsgw, &
!      &          zone,evec_ks(1:nmatp,ibgw:nbgw,ispn),nmatp, &
!      &          dSigma(ibgw:nbgw,n_qp),nbandsgw, &
!      &          zzero,evec_qp_right(1:nmatp,ib,2,ispn),nmatp)
! 	end if
! 	!-------------------------------------------------------------------------------------------------------------------------------------------------------------

    deallocate(selfec_nn,dSigma,zevec_left,evaldiag)

   end do !ib-----------------------------------------------------------------------------------------------------------------------------------------
!      
	deallocate(evec_overlap,evec_ks)!from overlap matrix
    !nbandsgw -> nbgw_diag-ibgw
    
    allocate(G_nmk(nbgw_diag,nbgw_diag,2))
    G_nmk(:,:,:)=zzero
    
    !spectral function, i.e. Im G_nmk = 1/pi < psi_n^L | psi_m^R > < psi_m^L | psi_n^R >  
    !->sum over m
    
     !< psi_n^L | psi_m^R > = conjg(VL)*conjg(c)*S*c*VR 
    call zgemm('n','n',nbgw_diag,nbgw_diag,nmatp, &
     &          zone,evec_qp_left(ibgw:nbgw_diag,1:nmatp,1,ispn),nbgw_diag, &
     &          evec_qp_right(1:nmatp,ibgw:nbgw_diag,1,ispn),nmatp, &
     &          zzero,G_nmk(:,:,1),nbgw_diag)

	G_nmk(n,m,1) = G_nmk(n,m,1) * G_nmk(m,n,1)
	do n=1,nbgw_diag
	do m=1,nbgw_diag
	 G_nmk(n,m,1) = G_nmk(n,m,1) * G_nmk(m,n,1)
	end do
	end do
	!print G_nm
    if (myrank==0) then
     call print_emat(ik,G_nmk(:,:,1),'eqp')
    end if
    
!satellites
!     ! G_sat_nm(n,m)
!     call zgemm('n','n',nbgw_diag,nbgw_diag,nmatp, &
!      &          zone,evec_qp_left(ibgw:nbgw_diag,1:nmatp,2,ispn),nbgw_diag, &
!      &          evec_qp_right(1:nmatp,ibgw:nbgw_diag,2,ispn),nmatp, &
!      &          zzero,G_nmk(:,:,2),nbgw_diag)
! 		
! 	do n=1,nbgw_diag
! 	do m=1,nbgw_diag
! 	 G_nmk(n,m,2) = G_nmk(m,n,2)*G_nmk(n,m,2)
! 	end do
! 	end do 
! 	!print G_nm
!     if (myrank==0) then
!      call print_emat(ik,G_nmk(:,:,2),'sat')
!     end if
!          
    !print eigenvalues
    if (myrank==0) then
     call print_eval(ik,eval_qp_diag(:,:))
    end if
    !print Im G_nk = sum_m G_nm*Gamma_m / [ (omega-E_m)^2 + Gamma_m^2 ]

    deallocate(G_nmk)


     end do !ik -----------------------------------------------------------------------------------------------------------------------------------------     

	deallocate(pmatvv_k)
	if ((input%gw%coreflag=='all') .or. (input%gw%coreflag=='xal')) deallocate(pmatcv_k) !.or. (input%gw%coreflag=='xal')
    !
    deallocate(evalks,evalqp)    
    deallocate(selfex_nn,vxc_nn)
    
    deallocate(evec_qp_left,evec_qp_right,eval_qp_diag)
        
    !clean up
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)

  end subroutine
  
!------------------------------------------------------------------------------
     !get new eval QP into eval_qp_diag(ib,i) from old one x+iy
     subroutine get_eval_qp_diag(ib,ikp,x,y,eta,f_ac,n_qp,i)
      real(8),intent(in) :: x,y,eta,f_ac
      integer,intent(in) :: i,ikp,ib
      integer,intent(out) :: n_qp
	  !local
      integer :: j

      !calculate Sigma^c_nn'k (E_n^QP) -> selfec_nn(n,n',k)
      selfec_nn(:,:) = zzero
      call calc_selfc_nn_wpol(ikp,x,y,eta,f_ac)
      !dSigma=Sigma^x+Sigma^c-V^XC  
      dSigma(:,:) = zzero
      zevec_left(:,:) = zzero
      dSigma(:,:) =  selfec_nn(:,:) + selfex_nn(:,:,ikp) - vxc_nn(:,:,ikp) 
      do j=ibgw,nbgw
       dSigma(j,j) = dSigma(j,j) + evalks(j,ikp)
      end do

      evaldiag(:) = zzero
      !a non-hermitian solver   
      call mkl_zgeev(nbandsgw,dSigma,evaldiag,zevec_left)
      !eigenvalues -> evaldiag(j), (right) eigenvectors -> dSigma(:,j) ,(left) eigenvectors -> zevec_left(:,j) 

      !n_qp = minloc(abs(x-REAL(evaldiag(:))),1)
      n_qp = minloc(abs(x+zi*y-evaldiag(:)),1)
      
      eval_qp_diag(ib,i) = evaldiag(n_qp)

      !if solution is better -> eval_qp_diag(ib,2) and calc vectors 
      if (abs(x+zi*y-evaldiag(n_qp)) < kappa ) then     
      
       kappa=abs(x+zi*y-evaldiag(n_qp))
      
       eval_qp_diag(ib,1) = x+zi*y
       eval_qp_diag(ib,2) = evaldiag(n_qp)
      
	   ! left EV: conjg(VL)*conjg(c_n)*S
       call zgemm('c','n',1,nmatp,nbandsgw, &
       &          zone,zevec_left(ibgw:nbgw,n_qp),nbandsgw, &
       &          evec_overlap(ibgw:nbgw,1:nmatp),nbandsgw, &
       &          zzero,evec_qp_left(ib,1:nmatp,1,ispn),1)
 	   !right EV: c_n*VR
       call zgemm('n','n',nmatp,1,nbandsgw, &
       &          zone,evec_ks(1:nmatp,ibgw:nbgw,ispn),nmatp, &
       &          dSigma(ibgw:nbgw,n_qp),nbandsgw, &
       &          zzero,evec_qp_right(1:nmatp,ib,1,ispn),nmatp)   

      end if
	     
     end subroutine


 subroutine get_satellites(iq)
	 
	 implicit none
	 integer, intent(in) :: iq
	 real(8) ::  n_el
	 integer :: i,j,s
	 
	 complex(8) :: eps_av
	 complex(8),allocatable :: epsil(:,:)
	 
	 !get epsilon^-1 = 1-w_s w^s^dagger/tau_s^2
	 Gamma = gammapoint(kqset%vqc(:,iq))
     ! Calculate interstitial product basis functions
     matsiz = locmatsiz+Gqset%ngk(1,iq)
     call diagsgi(iq)
     call calcmpwipw(iq)
     ! Calculate the bare Coulomb potential
     call calcbarcmb(iq)
	 call setbarcev(input%gw%barecoul%barcevtol)
     call get_wpol(iq)
     
     
     allocate(epsil(mbdim,mbdim))
     epsil(:,:)=zzero
     
     ! epsil0 = epsil0 - REAL(sum(wvck(j,s)*conjg(wvck(:,s))/tvck(s)**2))
     do i=1,1
     !do j=1,mbdim
      !if (i==j) epsil(i,j)=1.d0
      epsil(i,i)=1.d0
      do s=1,nvck  
       !epsil(i,j) = epsil(i,j)  - 2*wvck(i,s)*conjg(wvck(j,s))/tvck(s)
       epsil(i,i) = epsil(i,i)  - 2*wvck(i,s)*conjg(wvck(i,s))/tvck(s)
      enddo
     !enddo
     enddo

	 !clean up
     if (allocated(tvck)) deallocate(tvck)
     if (allocated(wvck)) deallocate(wvck)
	 if (allocated(mpwipw)) deallocate(mpwipw)
     if (allocated(barc)) deallocate(barc)

     eps_av=0.d0
	 ! average
	 do i=1,1
	  eps_av = eps_av + epsil(i,i) 
     enddo
     
     deallocate(epsil)
          
     !density
     !omega_pl ~ 16.7 eV ~ 0.6 Ha
	 !omega_pl = sqrt( n_el*e^2 / epsilon m_e ) 
	 
	 !IMPORTANT HARDCODED for SI 
	 n_el = 8.d0/omega
	 omega_pl = sqrt(abs(eps_av)*n_el*4*pi)	 
	 
 end subroutine	 
	 
!--------------------------------------------------------------------------------
  subroutine print_emat(ik,emat,cha)
    use m_getunit
    implicit none
    integer,    intent(in) :: ik
    !real(8), intent(in)    :: eqp,eks
    complex(8), intent(in) :: emat(nbgw_diag,nbgw_diag)
    character(3), intent(in) :: cha
    
    !local
    integer :: ib,ia
    integer :: fid
    character(80) :: fname

    ! store k,state-dependent evaldiag
    call getunit(fid)

    write(fname,'("emat_",A3,"-k",I4.4,".OUT")') cha,ik
    open(fid, File=trim(fname), Action='WRITE')
     do ia=ibgw,nbgw_diag
      do ib = ibgw, nbgw_diag
        write(fid,'(2f16.8)')  emat(ia,ib)
      end do
     !do ib = ibgw, nbgw
     !  write(fid,'(2f16.8)')  emat(ia,ib,2)
     !end do
      write(fid,*) 
     end do 
    close(fid)

    return
  end subroutine  
  
  
!--------------------------------------------------------------------------------
  subroutine print_eval(ik,eval)
    use m_getunit
    implicit none
    integer,    intent(in) :: ik
    !real(8), intent(in)    :: eqp,eks
    complex(8), intent(in) :: eval(nbgw_diag,7)
    
    !local
    integer :: ib,ia
    integer :: fid
    character(80) :: fname

    ! store k,state-dependent evaldiag
    call getunit(fid)

    write(fname,'("eval-k",I4.4,".OUT")') ik
    open(fid, File=trim(fname), Action='WRITE')
     !write(fid,'(2f16.6)') eqp,eks
      write(fid,'(5X,A,8X,A,8X,A,12X,A,12X,A,12X,A)') 'E^QP^(0)', 'gamma^(0)', '', '','',''
      do ib = ibgw, nbgw_diag
        !write(fid,'(10f16.8)')  evalks(ib,ik), eval(ib,1), eval(ib,2) , eval(ib,3) , real(eval(ib,3))-omega_pl, eval(ib,4)
        write(fid,'(14f16.8)')  eval(ib,1), eval(ib,2) , eval(ib,3) , eval(ib,4), eval(ib,5) , eval(ib,6), eval(ib,7)
!         write(fid,'(A)') '_____________________________'
!         do ia = ibgw, nbgw
!         write(fid,'(2f16.8)')  eval(ib,ia)
!         end do
!         write(fid,'(A)') '_____________________________'
      end do  
      write(fid,*)
    close(fid)

    return
  end subroutine
  
  
  subroutine mkl_zgeev(ndim,zevec,deval,zevec_left)
    implicit none
    integer,    intent(in)    :: ndim
    complex(8), intent(inout) :: zevec(ndim,ndim)
    complex(8), intent(out)   :: deval(ndim), zevec_left(ndim,ndim)
    ! local
    integer :: info, lwork
    
    complex(8) :: VL(ndim,ndim), VR(ndim,ndim)
    real(8),    allocatable :: rwork(:)
    complex(8), allocatable :: work(:)
    
    !
    write(*,*) 'Info(mod_zmatrix_diagonalization::mkl_zgeev):'
    write(*,*) '   Diagonalization with ZGEEV'
    write(*,*)
    !
    ! Query the optimal workspace
    !
    lwork = -1
    allocate(work(1), rwork(1))
    call zgeev( 'V', 'V' , ndim, zevec, ndim, deval, VL, ndim, VR, ndim, &
    &            work, lwork, rwork, info )
    lwork  = int(work(1))
    !lrwork = 3*ndim-2
    deallocate(work, rwork)
    !
    ! solve eigenproblem
    !
    allocate(work(lwork), rwork(2*ndim))
    call zgeev( 'V', 'V' , ndim, zevec, ndim, deval, VL, ndim, VR, ndim, &
    &            work, lwork, rwork, info )
    deallocate(work, rwork)
    
    !return right eigenvectors  lambda_i,VR(:,i) are the eigenvalues and -vectors
    zevec(:,:) = VR(:,:) 
    
    !return left eigenvectors lambda_i^*, VL(:,i)
    zevec_left(:,:) = VL(:,:) 
    
    
    if (info > 0) then
      write(*,*)
      write(*,*) 'Error(mod_zmatrix_diagonalization::mkl_zgeev):'
      write(*,*) '    ZGEEV algorithm failed to compute eigenvalues'
      write(*,*) '    info = ', info
      stop
    end if
    return
  end subroutine
  

  
end module

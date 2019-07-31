subroutine calc_selfx_GvvGG_wpol_nn_q(ikp,ik_,iq)

!!USES:
    use modmain,        only : nstfv, apwordmax, lmmaxapw, natmtot, nspnfv, &
    &                          pi, idxas, zzero, nmatmax, zone, occsv, zi
    use mod_gamma_correction
    
!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq, ikp, ik_ 
     !  real(8), intent(in) :: E_n_QP
    
!!LOCAL VARIABLES:            
    integer(4) :: ik, jk, iq_ , jkp1 ,jkp4 ,ikp_
    real(8)    :: tstart, tend, t0, t1
    integer(4) :: ie,iom, ie1, ie3, ie4, im
    integer(4) :: ia, is, ias, ic, icg
    real(8)    :: q0eps(3), modq0, enk1,enk3,enk4 , prefac ,de, om
    complex(8) :: sx, vc, sigma1, sigma2
    complex(8) :: mvm     ! Sum_ij{M^i*V^c_{ij}*dconjg(M^j)}
    complex(8), allocatable :: evecsv(:,:,:)
    complex(8), allocatable :: minmmat1(:,:,:), minmmat2(:,:,:) , minmmat3(:,:,:) , minmmat4(:,:,:)
   
    integer :: recl
   
    logical :: Gamma_ 
   
    integer :: ispn !k, l, ist, jst
    !complex(8) :: zsum
    
    ! external routines 
    complex(8), external :: zdotc    
    
    prefac = sqrt(4.d0*pi/omega) 
    
    !set reducible ik
    ik = kset%ikp2ik(ikp)
    !ik_ = kset%ikp2ik(ikp_) 
    ikp_ = kset%ik2ikp(ik_) 
    
    Gamma_ = gammapoint(kqset%vqc(:,iq))

	 jk = kqset%kqid(ik_,iq)
     jkp1 = kset%ik2ikp(jk)
    	
		
	!spinor sum, unecessary since nspinor=1 for unpolarized
    do ispn = 1, nspinor

	 !n_,n_4,k,q
     allocate(minmmat2(mbsiz,ibgw:nbgw,1:nstse))
     call calc_minmkq_gc(ik,iq,ispn,minmmat2,'nv') 

	 !n3,n1
     allocate(minmmat3(mbsiz,1:nstse,1:nstse))   
     !n_3,n_1,k',q
     call calc_minmkq_gc(ik_,iq,ispn,minmmat3,'vv')	

	 
     !loop über q'
    do iq_ = 1, kqset%nkpt
	
    jk = kqset%kqid(ik_,iq_)
    jkp4 = kset%ik2ikp(jk)
	
	Gamma = gammapoint(kqset%vqc(:,iq_))
    !========================================
    ! Calculate interstitial basis functions
    !========================================
    matsiz = locmatsiz+Gqset%ngk(1,iq_)
    call diagsgi(iq_)
    call calcmpwipw(iq_)
    !======================================
    ! Calculate the bare Coulomb potential
    !======================================
    call calcbarcmb(iq_)
	call setbarcev(input%gw%barecoul%barcevtol)

	
     ! read the momentum matrix elements
     if ((Gamma).and.(Gamma_)) then     
      ! val-val
      if (allocated(pmatvv)) deallocate(pmatvv)
      allocate(pmatvv(nomax,numin:nstse,3))
      inquire(iolength=recl) pmatvv
      open(fid_pmatvv,File=fname_pmatvv, &
      &    Action='READ',Form='UNFORMATTED',&
      &    Access='DIRECT',Status='OLD',Recl=recl)
     end if

     write(*,*) "ik ", ik,"ik_ ", ik_, "iq ", iq,"iq_ ", iq_

     !n,n1,k,q'
     allocate(minmmat1(mbsiz,ibgw:nbgw,1:nstse))
     call calc_minmkq_gc(ik,iq_,ispn,minmmat1,'nv')


     !n3,n4,k',q'  with n3=v, n4=c
     allocate(minmmat4(mbsiz,1:nstse,1:nstse))
     call calc_minmkq_gc(ik_,iq_,ispn,minmmat4,'vv')

!      ! read the momentum matrix elements
!      !if ((Gamma).or.(Gamma_)) call getpmatkgw(ik_)


do iom = 1, freq%nomeg   
 om = freq%freqs(iom)
 
     do ie = ibgw, nbgw
    
      sx = zzero
      sigma1 = zzero
      sigma2 = zzero
      
!sum over c1,v3,c4  __________________________________________________________________________________________       
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ie1,ie3,ie4,enk1,enk3,enk4,mvm)
!$OMP DO
#endif    
        do ie1 = numin, nstse ! unocc
        do ie3 = 1, nomax ! occ
        do ie4 = numin, nstse! unocc

        !KS energies   
        enk1 = evalsv(ie1,jkp1) 
        enk3 = evalsv(ie3,ikp_)
        enk4 = evalsv(ie4,jkp4)

          !calculate M* M* M M
          mvm = zdotc(mbsiz,minmmat4(:,ie3,ie4),1,minmmat1(:,ie,ie1),1) * &
          &     zdotc(mbsiz,minmmat2(:,ie,ie4),1,minmmat3(:,ie3,ie1),1) 
  
          mvm = mvm / (om - (enk1-enk3+enk4)  - zi*eta  )
          sx=sx + wkq*wkq*wkq*mvm

         ! add singular terms
         !singc1 -> 1/q
         !singc2 -> 1/q^2

         !q->0, q'->0 (no core treatment)
         if (Gamma.and.Gamma_.and.(ie==ie4).and.(ie==ie1)) then            
         ! MvM(q)
         de = enk3-enk1
         call gc_q0_treatment_setup(prefac/de,pmatvv(ie3,ie1,:)) !gets  pq_q = sqrt(4pi/Omega) p*q / delta E
         mvm = prefac * pq_q / (om - (enk1-enk3+enk4)  - zi*eta  )
         ! MvM(q')
         de = enk3-enk4
         call gc_q0_treatment_setup(prefac/de,pmatvv(ie3,ie4,:))
         sigma2 = sigma2 + singc2 * prefac * dconjg(pq_q) * mvm * wkq  ! wkq from k' sum
!         (q->0)
!         else if (Gamma_.and.(ie==ie4)) then
!           de = enk3-enk1
!           call gc_q0_treatment_setup(prefac/de,pmatvv(ie3,ie1,:))
!           mvm = zdotc(mbsiz,minmmat4(:,ie3,ie4),1,minmmat1(:,ie,ie1),1) / (om - (enk1-enk3+enk4)  - zi*eta  )
!           sigma1 = sigma1 - singc1 * prefac * pq_q * mvm * wkq * wkq !from q, k'   
!         !q'->0
!         else if (Gamma.and.(ie==ie1)) then
!           de = enk3-enk4
!           call gc_q0_treatment_setup(prefac/de,pmatvv(ie3,ie4,:))
!           mvm =  zdotc(mbsiz,minmmat2(:,ie,ie4),1,minmmat3(:,ie3,ie1),1) / (om - (enk1-enk3+enk4)  - zi*eta  )
!           sigma1 = sigma1 - singc1 * prefac * dconjg(pq_q) * mvm * wkq * wkq
        end if

      
      end do 
      end do 
      end do 

  
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

!
!sum over v1,c3,v4 !!  __________________________________________________________________________________________
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ie1,ie3,ie4,enk1,enk3,enk4,mvm)
!$OMP DO
#endif    
        do ie1 = 1, nomax
        do ie3 = numin, nstse
        do ie4 = 1, nomax

        !KS energies   
        enk1 = evalsv(ie1,jkp1)  ! evalks(j,ik)
        enk3 = evalsv(ie3,ikp_)
        enk4 = evalsv(ie4,jkp4)
          
        !calculate M* M* M M
         mvm = zdotc(mbsiz,minmmat4(:,ie3,ie4),1,minmmat1(:,ie,ie1),1) * &
         &     zdotc(mbsiz,minmmat2(:,ie,ie4),1,minmmat3(:,ie3,ie1),1) 
  
         mvm = mvm / (om - (enk1-enk3+enk4)  + zi*eta  )
         sx=sx + wkq*wkq*wkq*mvm

         !singc1 -> 1/q
         !singc2 -> 1/q^2

         ! add singular terms
     if (Gamma.and.Gamma_.and.(ie==ie4).and.(ie==ie1)) then                 
         ! MvM(q)
         de = enk3-enk1
         call gc_q0_treatment_setup(prefac/de,pmatvv(ie3,ie1,:)) !gets  pq_q = sqrt(4pi/Omega) p*q / delta E
         mvm = prefac * pq_q / (om - (enk1-enk3+enk4)  + zi*eta  )
         ! MvM(q')
         de = enk3-enk4
         call gc_q0_treatment_setup(prefac/de,pmatvv(ie3,ie4,:))
         sigma2 = sigma2 + singc2 * prefac * dconjg(pq_q) * mvm * wkq  ! wkq from k' sum
        end if

      end do 
      end do 
      end do 
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

    selfex_GvvGG_nn(ie,iom) = selfex_GvvGG_nn(ie,iom) + sx  + sigma2 !+ sigma1

  end do ! ie!
 end do ! iom	


    if ((Gamma).and.(Gamma_)) then 
      close(fid_pmatvv)
      deallocate(pmatvv)
    end if

    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
    
	deallocate(minmmat1,minmmat4)
	
    end do  !q'
	
	deallocate(minmmat2,minmmat3)	
	
	end do !ispn
    
	
	
	
    return

 end subroutine 

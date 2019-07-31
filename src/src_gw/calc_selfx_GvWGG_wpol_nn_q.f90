subroutine calc_selfx_GvWGG_wpol_nn_q(ikp,ik_,iq)

!!USES:
    use modmain,        only : nstfv, apwordmax, lmmaxapw, natmtot, nspnfv, &
    &                          pi, idxas, zzero, nmatmax, zone, occsv, zi
    use mod_gamma_correction
    
!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq, ikp,ik_
    
!!LOCAL VARIABLES:            
    !integer(4) :: ik, jkp, jk, ik_, iq_,jkp4,jkp1
    integer(4) :: ik, jk, iq_ , jkp1 ,jkp4 ,ikp_
    real(8)    :: tstart, tend, t0, t1
    integer(4) :: ie,iom, ie1, ie3, ie4, im, s1
    integer(4) :: ia, is, ias, ic, icg
    real(8)    :: q0eps(3), modq0, enk1,enk3,enk4 , prefac ,  om, theta0=1.d-6 !numerical zero
    
    complex(8) :: sx,  mvm ,vc , zt1 !, sigma1, sigma2
    complex(8), allocatable :: evecsv(:,:,:)
    complex(8), allocatable :: minmmat1(:,:,:), minmmat2(:,:,:) , minmmat3(:,:,:) , minmmat4(:,:,:)
    complex(8), allocatable :: mw(:,:),mw_(:,:), mwt(:)
    
    integer :: recl
   
    logical :: Gamma_ 
   
    integer :: ispn !k, l, ist, jst

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


    write(*,*) "ik ", ik,"ik_ ", ik_, "iq ", iq,"iq_ ", iq_

     !n,n1,k,q'
     allocate(minmmat1(mbsiz,ibgw:nbgw,1:nstse))
     call calc_minmkq_gc(ik,iq_,ispn,minmmat1,'nv')

     !n3,n4,k',q'  with n3=v, n4=c
     allocate(minmmat4(mbsiz,1:nstse,1:nstse))
     call calc_minmkq_gc(ik_,iq_,ispn,minmmat4,'vv')

    allocate(mw(nstse,nvck))
    allocate(mw_(nstse,nvck))
    allocate(mwt(nvck))

    
   do ie = ibgw, nbgw
      
     !M_n n4 (kq)^* w_s1  
     call zgemm( 'c', 'n', nomax, nvck, mbsiz, &
      &           zone, minmmat2(1:mbsiz,ie,1:nomax), mbsiz, &
      &           wvck(1:mbsiz,1:nvck), mbsiz, &
      &           zzero, mw, nomax)
        
   !freq sum
   do iom = 1, freq%nomeg   
    om = freq%freqs(iom)   
   
    sx=zzero

 !sum over v_1,s,c_3,n_4 __________________________________________________________________________________________

      do ie3 = numin, nstse !unocc

        !M_n3 n1 (k'q)^* w_s1  
        call zgemm( 'c', 'n', nomax, nvck, mbsiz, &
        &           zone, minmmat3(1:mbsiz,ie3,1:nomax), mbsiz, &
        &           wvck(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw_, nomax)

        enk3 = evalsv(ie3,ikp_)
        
      do ie1 = 1, nomax ! occ 
         enk1 = evalsv(ie1,jkp1)  ! evalks(j,ik)
         
      do ie4 = numin, nstse
       enk4 = evalsv(ie4,jkp4)
       mvm = zdotc(mbsiz,minmmat4(:,ie3,ie4),1,minmmat1(:,ie,ie1),1)  
       mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
       do s1 = 1, nvck    
       
		zt1 = 1 / (enk1-enk3-tvck(s1)) * 1 / ( om - (enk1-enk3+enk4 ) + zi*eta) 
		
		if (abs(enk1-enk3-tvck(s1)) < theta0 ) zt1 = zzero
		
        mwt(s1) =  zt1*mvm*mw(ie4,s1)
       end do 
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif  
       !sum vck
       sx = sx - wkq*wkq*wkq*zdotc(nvck,mw_(ie1,:),1,mwt,1) 
      end do  !ie4 

       do ie4 = 1, nomax! occ
        enk4 = evalsv(ie4,jkp4)
        mvm = zdotc(mbsiz,minmmat4(:,ie3,ie4),1,minmmat1(:,ie,ie1),1)  
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
       do s1 = 1, nvck    
       
		zt1 = 1 / (enk1-enk3-tvck(s1)) * 1 / ( om - (enk4+tvck(s1) ) - zi*eta)  &
		&   + 1 / (enk1-enk3+tvck(s1)) * ( 1 / ( om - (enk4-tvck(s1) ) - zi*eta) - 1 / ( om - (enk1-enk3+enk4) - zi*eta) ) 
		
        if (abs(enk1-enk3-tvck(s1)) < theta0 ) zt1 = 1 / (enk1-enk3+tvck(s1)) * ( 1 / ( om - (enk4-tvck(s1) ) - zi*eta) - 1 / ( om - (enk1-enk3+enk4) - zi*eta) ) 
        if (abs(enk1-enk3+tvck(s1)) < theta0 ) zt1 = 1 / (enk1-enk3-tvck(s1)) * 1 / ( om - (enk4+tvck(s1) ) - zi*eta)		
		
        mwt(s1) =  zt1*mvm*mw(ie4,s1)
       end do 
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif  
       !sum vck
       sx = sx - wkq*wkq*wkq*zdotc(nvck,mw_(ie1,:),1,mwt,1) 
      end do  !ie4 

     end do !ie1
     end do !ie3
      
     
  
!sum over c_1,s,v_3,n_4 __________________________________________________________________________________________
 
     do ie3 = 1, nomax 

        !M_n3 n1 (k'q)^* w_s1  
        call zgemm( 'c', 'n', nomax, nvck, mbsiz, &
        &           zone, minmmat3(1:mbsiz,ie3,1:nomax), mbsiz, &
        &           wvck(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw_, nomax)

        enk3 = evalsv(ie3,ikp_)
        
     do ie1 = numin, nstse  
        enk1 = evalsv(ie1,jkp1)  ! evalks(j,ik)
         
     do ie4 = 1, nomax
       enk4 = evalsv(ie4,jkp4)
       mvm = zdotc(mbsiz,minmmat4(:,ie3,ie4),1,minmmat1(:,ie,ie1),1)  
       mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
       do s1 = 1, nvck    
       
		zt1 = 1 / (enk1-enk3+tvck(s1)) * 1 / ( om - (enk1-enk3+enk4) - zi*eta)   
		
		if (abs(enk1-enk3+tvck(s1)) < theta0 ) zt1=zzero
		
        mwt(s1) =  zt1*mvm*mw(ie4,s1)
       end do 
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif  
       !sum vck
       sx = sx + wkq*wkq*wkq*zdotc(nvck,mw_(ie1,:),1,mwt,1) 
    end do  !ie4

     do ie4 = numin,nstse
       enk4 = evalsv(ie4,jkp4)
       mvm = zdotc(mbsiz,minmmat4(:,ie3,ie4),1,minmmat1(:,ie,ie1),1)  
       mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
       do s1 = 1, nvck    
              
		zt1 = 1 / (enk1-enk3+tvck(s1)) * 1 / ( om - (enk4-tvck(s1) ) + zi*eta)  &
		&   + 1 / (enk1-enk3-tvck(s1)) * ( 1 / ( om - (enk4+tvck(s1) ) + zi*eta) - 1 / ( om - (enk1-enk3+enk4) + zi*eta) ) 
		
        if (abs(enk1-enk3-tvck(s1)) < theta0 ) zt1 = 1 / (enk1-enk3+tvck(s1)) * 1 / ( om - (enk4-tvck(s1) ) + zi*eta)
        if (abs(enk1-enk3+tvck(s1)) < theta0 ) zt1 = 1 / (enk1-enk3-tvck(s1)) * ( 1 / ( om - (enk4+tvck(s1) ) + zi*eta) - 1 / ( om - (enk1-enk3+enk4) + zi*eta) ) 
		
        mwt(s1) =  zt1*mvm*mw(ie4,s1)
       end do 
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif  
       !sum vck
       sx = sx + wkq*wkq*wkq*zdotc(nvck,mw_(ie1,:),1,mwt,1)       
      
       end do  !ie4
      
       end do !ie1
      end do !ie3
        
      selfex_GvWGG_nn(ie,iom) = selfex_GvWGG_nn(ie,iom) + sx  

    end do !iom
   end do !ie



deallocate(mw,mw_,mwt)
    
	deallocate(minmmat1,minmmat4)
	
    end do  !q'
	
	deallocate(minmmat2,minmmat3)	
	
	end do !ispn


    return

end subroutine 

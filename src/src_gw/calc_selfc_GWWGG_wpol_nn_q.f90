subroutine calc_selfc_GWWGG_wpol_nn_q(ikp,ik_,iq)

!!USES:
    use modmain,        only : nstfv, apwordmax, lmmaxapw, natmtot, nspnfv, &
    &                          pi, idxas, zzero, nmatmax, zone, occsv, zi
    use mod_gamma_correction
    use ieee_arithmetic
    
!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq, ikp,ik_
    
!!LOCAL VARIABLES:            
    !integer(4) :: ik, jkp, jk, ik_, iq_,jkp4,jkp1
    integer(4) :: ik, jk, iq_ , jkp1 ,jkp4 ,ikp_
    real(8)    :: tstart, tend, t0, t1
    integer(4) :: ie,iom, ie1, ie3, ie4, im, s1, s2
    integer(4) :: ia, is, ias, ic, icg
    real(8)    :: q0eps(3), modq0, enk1,enk3,enk4 , prefac ,de, om
    
    real(8)    ::  wkq_, theta0=1.d-5 !numerical zero
    
    complex(8) :: sx,  mvm ,vc , zt1, mwm   !, sigma1, sigma2
    complex(8), allocatable :: evecsv(:,:,:), wvck_q(:,:), tvck_q(:)
    complex(8), allocatable :: minmmat1(:,:,:), minmmat2(:,:,:) , minmmat3(:,:,:) , minmmat4(:,:,:)
    complex(8), allocatable :: mw1(:,:),mw2(:,:),mw3(:,:),mw4(:,:), mwt(:)
    
    integer :: recl
   
    logical :: Gamma_ 
   
    integer :: ispn !k, l, ist, jst

    ! external routines 
    complex(8), external :: zdotc 

    prefac = sqrt(4.d0*pi/omega) 
    
    wkq_=wkq*wkq*wkq 
    
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
     
     !save W(q)
     allocate(wvck_q(mbdim,nvck))
     allocate(tvck_q(nvck))
     wvck_q(:,:) = wvck(:,:)
     tvck_q(:) = tvck(:)

    !loop over q'
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

    !get tvck,wvck corresponding to q'
    call get_wpol(iq_)
	
    write(*,*) "ik ", ik,"ik_ ", ik_, "iq ", iq,"iq_ ", iq_

     !n,n1,k,q'
     allocate(minmmat1(mbsiz,ibgw:nbgw,1:nstse))
     call calc_minmkq_gc(ik,iq_,ispn,minmmat1,'nv')

     !n3,n4,k',q'  with n3=v, n4=c
     allocate(minmmat4(mbsiz,1:nstse,1:nstse))
     call calc_minmkq_gc(ik_,iq_,ispn,minmmat4,'vv')


    allocate(mw1(nstse,nvck))
    allocate(mw2(nstse,nvck))
    allocate(mw3(nstse,nvck))
    allocate(mw4(nstse,nvck))
    allocate(mwt(nvck))


do ie = ibgw, nbgw
   
   call zgemm( 'c', 'n', nomax, nvck, mbsiz, &
        &           zone, minmmat1(1:mbsiz,ie,1:nomax), mbsiz, &
        &           wvck(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw1, nomax)
   call zgemm( 'c', 'n', nomax, nvck, mbsiz, &
        &           zone, minmmat2(1:mbsiz,ie,1:nomax), mbsiz, &
        &           wvck_q(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw2, nomax)
		
 do iom = 1, freq%nomeg   
  om = freq%freqs(iom)
 
  sx = zzero
  
!________________________________________________________________________________________________________________________________________________________________________   
!ie3 occ ________________________________________________________________________________________________________________________________________________________________________         
!________________________________________________________________________________________________________________________________________________________________________          
 do ie3 = 1,nomax 
    call zgemm( 'c', 'n', nomax, nvck, mbsiz, &
        &           zone, minmmat3(1:mbsiz,ie3,1:nomax), mbsiz, &
        &           wvck_q(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw3, nomax)
    call zgemm( 'c', 'n', nomax, nvck, mbsiz, &
        &           zone, minmmat4(1:mbsiz,ie3,1:nomax), mbsiz, &
        &           wvck(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw4, nomax)
    enk3 = evalsv(ie3,ikp_)
      
      
!sum over v_1,s_2,s_1,v_3,v_4 __________________________________________________________________________________________  
   do ie1 = 1,nomax 
    enk1 = evalsv(ie1,jkp1) 
   do ie4 = 1,nomax 
    enk4 = evalsv(ie4,jkp4)
   do s2 = 1, nvck 
     mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
     mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
     do s1 = 1, nvck
 
        zt1 = 1 / (tvck_q(s1)-enk3+enk1) * 1 / (enk3-enk4-tvck(s2))          * ( 1 / ( om - (enk3-tvck_q(s1)-tvck(s2)) - zi*eta) -  1 / ( om - (enk4-tvck_q(s1)) - zi*eta) )   &        
          & + 1 / (tvck_q(s1)-enk3+enk1) * 1 / (tvck_q(s1)-tvck(s2)+enk1-enk4) * ( 1 / ( om - (enk4-tvck_q(s1)) - zi*eta)          -  1 / ( om - (enk1-tvck(s2)) - zi*eta) )  
        
	    !exclude possible NaN and INF situations 
        if (abs(tvck_q(s1)-enk3+enk1) < theta0 ) cycle
        if (abs(enk3-enk4-tvck(s2))  < theta0 ) cycle
        
        !exclude zero term, where the prefactor leads to NaN otherwise
        if (abs(tvck_q(s1)-tvck(s2)+enk1-enk4) < theta0 ) zt1 = 1 / (tvck_q(s1)-enk3+enk1) * 1 / (enk3-enk4-tvck(s2)) * ( 1 / ( om - (enk3-tvck_q(s1)-tvck(s2)) - zi*eta) -  1 / ( om - (enk4-tvck_q(s1)) - zi*eta) ) 
        
        mwt(s1) = zt1*mw2(ie4,s1)*mwm  
    end do !s1
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 

     end do !s2
     end do !ie4
     end do !ie1
     
!sum over v_1,s_2,s_1/s_1',v_3,c_4 _______________________________________________________________________________________________________________________________________________________________________________    
     do ie1 = 1,nomax
      enk1 = evalsv(ie1,jkp1) 
     do ie4 = numin,nstse
      enk4 = evalsv(ie4,jkp4)
     do s2 = 1, nvck
      mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
      mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
    do s1 = 1, nvck    
    
		if (abs(tvck(s2)-enk3+enk4) < theta0 ) cycle
        if (abs(tvck(s2)-tvck_q(s1)+enk4-enk1)  < theta0 ) cycle
        if (abs(tvck(s2)+tvck_q(s1)+enk4-enk1) < theta0 ) cycle
        
        zt1 = 1 / (tvck(s2)-enk3+enk4) * 1 / (tvck(s2)-tvck_q(s1)+enk4-enk1)  &
        &   * ( 1 / ( om - (enk4-tvck_q(s1)) - zi*eta) -  1 / ( om - (enk1-tvck(s2)) - zi*eta) )  &
        &   + 1 / (tvck(s2)-enk3+enk4) * 1 / (tvck(s2)+tvck_q(s1)+enk4-enk1)  &
        &   * ( 1 / ( om - (enk4+tvck_q(s1)) - zi*eta) -  1 / ( om - (enk1-tvck(s2)) + zi*eta) ) 
        mwt(s1) = zt1*mw2(ie4,s1)*mwm  
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
      end do !s2 
      end do !ie4
      end do !ie1
     
!sum over c_1,s_2/s_2',s_1,v_3,v_4 _______________________________________________________________________________________________________________________________________________________________________________    
      do ie1 = numin,nstse
        enk1 = evalsv(ie1,jkp1) 
      do ie4 = 1,nomax
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck
       mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2))  
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
   do s1 = 1, nvck 
   
   		if (abs(tvck_q(s1)-enk3+enk1) < theta0 ) cycle
        if (abs(tvck_q(s1)-tvck(s2)+enk1-enk4) < theta0 ) cycle
        if (abs(tvck_q(s1)+tvck(s2)+enk1-enk4) < theta0 ) cycle
        
        zt1 = 1 / (tvck_q(s1)-enk3+enk1) * 1 / (tvck_q(s1)-tvck(s2)+enk1-enk4) &
        &   * ( 1 / ( om - (enk1-tvck(s2)) - zi*eta) -  1 / ( om - (enk4-tvck_q(s1)) - zi*eta) )  &
        &   + 1 / (tvck_q(s1)-enk3+enk1) * 1 / (tvck_q(s1)+tvck(s2)+enk1-enk4) &
        &   * ( 1 / ( om - (enk1+tvck(s2)) - zi*eta) -  1 / ( om - (enk4-tvck_q(s1)) + zi*eta) ) 
        
        mwt(s1) = zt1*mw2(ie4,s1)*mwm  
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
      end do !s2
      end do !ie4
      end do !ie1
     
!sum over  c_1,s_2,s_1/s_1',v_3,c_4 _______________________________________________________________________________________________________________________________________________________________________________    
      do ie1 = numin,nstse
       enk1 = evalsv(ie1,jkp1) 
      do ie4 = numin,nstse
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck
       mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
       mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
   do s1 = 1, nvck 
   
		if (abs(tvck(s2)-enk3+enk4) < theta0 ) cycle
        if (abs(tvck(s2)+tvck_q(s1)+enk4-enk1) < theta0 ) cycle
     
        zt1 = 1 / (tvck(s2)-enk3+enk4) * 1 / (tvck(s2)-tvck_q(s1)+enk4-enk1) &
        &   * ( 1 / ( om - (enk1-tvck(s2)) - zi*eta) -  1 / ( om - (enk4-tvck_q(s1)) + zi*eta) )  &
        &   + 1 / (tvck(s2)-enk3+enk4) * 1 / (tvck(s2)+tvck_q(s1)+enk4-enk1) &
        &   * ( 1 / ( om - (enk1-tvck(s2)) + zi*eta) -  1 / ( om - (enk4+tvck_q(s1)) + zi*eta) )  
        
        if (abs(tvck(s2)-tvck_q(s1)+enk4-enk1) < theta0 ) zt1 = 1 / (tvck(s2)-enk3+enk4) * 1 / (tvck(s2)+tvck_q(s1)+enk4-enk1) * ( 1 / ( om - (enk1-tvck(s2)) + zi*eta) -  1 / ( om - (enk4+tvck_q(s1)) + zi*eta) ) 
        
        mwt(s1) = zt1*mw2(ie4,s1)*mwm  
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
      end do 
      end do !ie4    
      end do !ie1
     

!sum over  c_1,s_2',s_1/s_1',v_3,c_4 _______________________________________________________________________________________________________________________________________________________________________________    
      do ie1 = numin,nstse
       enk1 = evalsv(ie1,jkp1) 
      do ie4 = numin,nstse
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck
        mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
     do s1 = 1, nvck 
     
        if (abs(enk3-enk1-tvck_q(s1)) < theta0 ) cycle
        if (abs(tvck(s2)+tvck_q(s1)+enk1-enk4) < theta0 ) cycle
        if (abs(enk3-enk1+tvck_q(s1)) < theta0 ) cycle

        zt1 = 1 / (enk3-enk1-tvck_q(s1)) * 1 / (tvck(s2)+tvck_q(s1)+enk1-enk4) &
        &   * ( 1 / ( om - (enk1+tvck(s2)) + zi*eta) -  1 / ( om - (enk4-tvck_q(s1)) + zi*eta) )  &
        &   + 1 / (enk3-enk1+tvck_q(s1)) * 1 / (tvck(s2)-tvck_q(s1)+enk1-enk4) &
        &   * ( 1 / ( om - (enk1+tvck(s2)) + zi*eta) -  1 / ( om - (enk4+tvck_q(s1)) + zi*eta) )  
        
        if (abs(tvck(s2)-tvck_q(s1)+enk1-enk4) < theta0 ) zt1 = 1 / (enk3-enk1-tvck_q(s1)) * 1 / (tvck(s2)+tvck_q(s1)+enk1-enk4) * ( 1 / ( om - (enk1+tvck(s2)) + zi*eta) -  1 / ( om - (enk4-tvck_q(s1)) + zi*eta) ) 
        
        mwt(s1) = zt1*mw2(ie4,s1)*mwm
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
      end do 
      end do !ie4
      end do !ie1

!sum over  c_1,s_2',s_1',v_3,c_4 _______________________________________________________________________________________________________________________________________________________________________________    
      do ie1 = numin,nstse
       enk1 = evalsv(ie1,jkp1) 
      do ie4 = numin,nstse
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck    
        mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
     do s1 = 1, nvck 
     
        if (abs(enk3-enk1+tvck_q(s1)) < theta0 ) cycle
        if (abs(tvck(s2)+enk3-enk4)  < theta0 ) cycle
        
      zt1 = 1 / (enk3-enk1+tvck_q(s1)) * 1 / (tvck(s2)+enk3-enk4) &
      &   * ( 1 / ( om - (enk1-enk3+enk4) + zi*eta) -  1 / ( om - (enk1+tvck(s2)) + zi*eta) )  
      mwt(s1) = zt1*mw2(ie4,s1)*mwm  
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
     end do 
     end do !ie4
     end do !ie1
   
 end do !ie3-occ
     
!________________________________________________________________________________________________________________________________________________________________________   
!ie3 unocc ________________________________________________________________________________________________________________________________________________________________________         
!________________________________________________________________________________________________________________________________________________________________________   
 do ie3 = numin,nstse  
    call zgemm( 'c', 'n', nomax, nvck, mbsiz, &
        &           zone, minmmat3(1:mbsiz,ie3,1:nomax), mbsiz, &
        &           wvck(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw3, nomax)
    call zgemm( 'c', 'n', nomax, nvck, mbsiz, &
        &           zone, minmmat4(1:mbsiz,ie3,1:nomax), mbsiz, &
        &           wvck(1:mbsiz,1:nvck), mbsiz, &
        &           zzero, mw4, nomax)
    enk3 = evalsv(ie3,ikp_)

!sum over c_1,s_2',s_1',c_3,c_4 _______________________________________________________________________________________________________________________________________________________________________________    

      do ie1 = numin,nstse 
        enk1 = evalsv(ie1,jkp1) 
      do ie4 = numin,nstse
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck
        mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
   do s1 = 1, nvck 
   
        if (abs(tvck(s2)-enk4+enk3) < theta0 ) cycle
        if (abs(tvck_q(s1)-enk1+enk3) < theta0 ) cycle

        zt1 = 1 / (tvck(s2)-enk4+enk3) * 1 / (tvck_q(s1)-enk1+enk3) & 
        &   * ( 1 / ( om - (enk4+tvck_q(s1)) + zi*eta) - 1 / ( om - (enk3+tvck_q(s1)+tvck(s2)) + zi*eta) )  & 
        &   + 1 / (tvck_q(s1)-enk1+enk3) * 1 / (tvck_q(s1)-tvck(s2)+enk4-enk1) &
        &   * ( 1 / ( om - (enk4+tvck_q(s1)) + zi*eta) - 1 / ( om - (enk1+tvck(s2)) + zi*eta) )  
        
        if (abs(tvck_q(s1)-tvck(s2)+enk4-enk1) < theta0 ) zt1 = 1 / (tvck(s2)-enk4+enk3) * 1 / (tvck_q(s1)-enk1+enk3) * ( 1 / ( om - (enk4+tvck_q(s1)) + zi*eta) - 1 / ( om - (enk3+tvck_q(s1)+tvck(s2)) + zi*eta) ) 
        
        mwt(s1) = zt1*mw2(ie4,s1)*mwm 
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
      end do 
    end do 
   end do 
   
!sum over c_1,s_2',s_1/s_1',c_3,v_4 _______________________________________________________________________________________________________________________________________________________________________________    
      do ie1 = numin,nstse
        enk1 = evalsv(ie1,jkp1) 
      do ie4 = 1,nomax
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck
        mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
   do s1 = 1, nvck        
   
		if (abs(tvck(s2)-enk4+enk3) < theta0 ) cycle
        if (abs(tvck(s2)+tvck_q(s1)+enk1-enk4)  < theta0 ) cycle
        if (abs(tvck(s2)-tvck_q(s1)+enk1-enk4)  < theta0 ) cycle
        
        zt1 = 1 / (tvck(s2)-enk4+enk3) * 1 / (tvck(s2)+tvck_q(s1)+enk1-enk4)  &
        &   * ( 1 / ( om - (enk4-tvck_q(s1)) + zi*eta) -  1 / ( om - (enk1+tvck(s2)) - zi*eta) )  &
        &   + 1 / (tvck(s2)-enk4+enk3) * 1 / (tvck(s2)-tvck_q(s1)+enk1-enk4)  &
        &   * ( 1 / ( om - (enk4+tvck_q(s1)) + zi*eta) -  1 / ( om - (enk1+tvck(s2)) + zi*eta) ) 
        mwt(s1) = zt1*mw2(ie4,s1)*mwm 
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
     end do 
     end do 
     end do 

!sum over v_1,s_2/s_2',s_1',c_3,c_4 _______________________________________________________________________________________________________________________________________________________________________________    
      do ie1 = 1,nomax
        enk1 = evalsv(ie1,jkp1) 
      do ie4 = numin,nstse
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck
        mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
   do s1 = 1, nvck 
   
   		if (abs(enk1-enk3-tvck_q(s1)) < theta0 ) cycle
        if (abs(tvck_q(s1)+tvck(s2)+enk4-enk1)  < theta0 ) cycle
        if (abs(tvck_q(s1)-tvck(s2)+enk4-enk1)  < theta0 ) cycle
   
        zt1 = 1 / (enk1-enk3-tvck_q(s1)) * 1 / (tvck_q(s1)+tvck(s2)+enk4-enk1) &
        &   * ( 1 / ( om - (enk4+tvck_q(s1)) - zi*eta) -  1 / ( om - (enk1-tvck(s2)) + zi*eta) )  &
        &   + 1 / (enk1-enk3-tvck_q(s1)) * 1 / (tvck_q(s1)-tvck(s2)+enk4-enk1) &
        &   * ( 1 / ( om - (enk4+tvck_q(s1)) + zi*eta) -  1 / ( om - (enk1+tvck(s2)) + zi*eta) ) 
        mwt(s1) = zt1*mw2(ie4,s1)*mwm
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
     end do 
     end do 
     end do 

!sum over  v_1,s_2/s_2',s_1',c_3,v_4 _______________________________________________________________________________________________________________________________________________________________________________    
      do ie1 = 1,nomax 
        enk1 = evalsv(ie1,jkp1) 
      do ie4 = 1,nomax 
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck
        mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
   do s1 = 1, nvck 
   
    	if (abs(enk1-enk3-tvck_q(s1)) < theta0 ) cycle
        if (abs(tvck_q(s1)+tvck(s2)+enk4-enk1)  < theta0 ) cycle
        
        zt1 = 1 / (enk1-enk3-tvck_q(s1)) * 1 / (tvck_q(s1)+tvck(s2)+enk4-enk1) &
        &   * ( 1 / ( om - (enk1-tvck(s2)) - zi*eta) -  1 / ( om - (enk4+tvck_q(s1)) - zi*eta) )  &
        &   +  1 / (enk1-enk3-tvck_q(s1)) * 1 / (tvck_q(s1)-tvck(s2)+enk4-enk1) &
        &   * ( 1 / ( om - (enk1+tvck(s2)) - zi*eta) -  1 / ( om - (enk4+tvck_q(s1)) + zi*eta) )  
        
        if (abs(tvck_q(s1)-tvck(s2)+enk4-enk1) < theta0 ) zt1 = 1 / (enk1-enk3-tvck_q(s1)) * 1 / (tvck_q(s1)+tvck(s2)+enk4-enk1) * ( 1 / ( om - (enk1-tvck(s2)) - zi*eta) -  1 / ( om - (enk4+tvck_q(s1)) - zi*eta) )
        
        mwt(s1) = zt1*mw2(ie4,s1)*mwm
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
     end do 
     end do 
     end do 
     
!sum over  v_1,s_2/s_2',s_1,c_3,v_4 _______________________________________________________________________________________________________________________________________________________________________________  
      do ie1 = 1,nomax 
        enk1 = evalsv(ie1,jkp1) 
      do ie4 = 1,nomax 
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck
        mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
   do s1 = 1, nvck 
   
		if (abs(enk4-enk3+tvck(s2)) < theta0 ) cycle
        if (abs(enk4-enk3-tvck(s2)) < theta0 ) cycle
        if (abs(tvck_q(s1)+tvck(s2)+enk1-enk4)  < theta0 ) cycle
        
        zt1 = 1 / (enk4-enk3+tvck(s2)) * 1 / (tvck_q(s1)-tvck(s2)+enk1-enk4) &
        &   * ( 1 / ( om - (enk4-tvck_q(s1)) - zi*eta) -  1 / ( om - (enk1-tvck(s2)) - zi*eta) )  &
        &   +  1 / (enk4-enk3-tvck(s2)) * 1 / (tvck_q(s1)+tvck(s2)+enk1-enk4) &
        &   * ( 1 / ( om - (enk4-tvck_q(s1)) - zi*eta) -  1 / ( om - (enk1+tvck(s2)) - zi*eta) )  
        
        if (abs(tvck_q(s1)-tvck(s2)+enk1-enk4) < theta0 ) zt1 = 1 / (enk4-enk3-tvck(s2)) * 1 / (tvck_q(s1)+tvck(s2)+enk1-enk4) * ( 1 / ( om - (enk4-tvck_q(s1)) - zi*eta) -  1 / ( om - (enk1+tvck(s2)) - zi*eta) )
        
        mwt(s1) = zt1*mw2(ie4,s1)*mwm
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
     end do 
     end do 
     end do 

!sum over  v_1,s_2,s_1,c_3,v_4 _______________________________________________________________________________________________________________________________________________________________________________    
      do ie1 = 1,nomax 
        enk1 = evalsv(ie1,jkp1) 
      do ie4 = 1,nomax 
       enk4 = evalsv(ie4,jkp4)
      do s2 = 1, nvck
        mwm = mw4(ie4,s2)*dconjg(mw1(ie1,s2)) 
        mwt=zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1,zt1)
!$OMP DO
#endif    
     do s1 = 1, nvck 
     
		if (abs(enk1-enk3+tvck_q(s1))  < theta0 ) cycle
        if (abs(enk3-enk4-tvck(s2))  < theta0 ) cycle
     
        zt1 = 1 / (enk1-enk3+tvck_q(s1)) * 1 / (enk3-enk4-tvck(s2)) &
        &   * ( 1 / ( om - (enk1-tvck(s2)) - zi*eta) -  1 / ( om - (enk1-enk3+enk4) - zi*eta) )  
        mwt(s1) = zt1*mw2(ie4,s1)*mwm 
     end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
       !sum vck=s1     
       sx = sx - wkq_*zdotc(nvck,mw3(ie1,:),1,mwt,1) 
     end do 
     end do 
     end do 
   
 end do !ie3-unocc
     
 selfec_GWWGG_nn(ie,iom) = selfec_GWWGG_nn(ie,iom) + sx

 end do !iom
end do   !ie


     deallocate(mw1,mw2,mw3,mw4,mwt)
	 deallocate(minmmat1,minmmat4)
	 
	 !W(q')
     deallocate(tvck)
     deallocate(wvck)

    end do  !q'
	
	 deallocate(minmmat2,minmmat3)	
	 !W(q)
 	 deallocate(tvck_q)
     deallocate(wvck_q)
	
	end do !ispn


    return

end subroutine 

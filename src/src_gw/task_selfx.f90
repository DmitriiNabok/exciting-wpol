
subroutine task_selfx()

  use modinput
  use modmain,               only : zzero, evalsv, efermi
  use modgw
  use mod_mpi_gw
  use m_getunit
            
  implicit none
  integer(4) :: ikp, iq, fid, ik
  integer(4) :: recl
  integer :: ie

  !===========================================================================
  ! Initialization
  !===========================================================================
  call init_gw
  call clean_gndstate
    
  ! occupancy dependent BZ integration weights
  call kintw

  if (allocated(selfex)) deallocate(selfex)
  allocate(selfex(ibgw:nbgw,nkpt))
  selfex(:,:) = 0.d0

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
    write(*,*) '(task_selfx): q-point cycle, iq = ', iq
    
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
    call calcselfx(iq)
      
    ! clean unused data
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
      
  end do ! iq
    
  if (allocated(kiw)) deallocate(kiw)
  if (allocated(ciw)) deallocate(ciw)
    
#ifdef MPI
    call mpi_sum_array(0,selfex,nbandsgw,kset%nkpt,mycomm_row)
#endif
    
  if (myrank==0) call write_selfexnn(kset,ibgw,nbgw,selfex)
  
  ! clear memory  
  deallocate(selfex)
  call delete_freqgrid(freq)
  call delete_k_vectors(kset)
  call delete_G_vectors(Gset)
  call delete_Gk_vectors(Gkset)
  call delete_kq_vectors(kqset)
  call delete_Gk_vectors(Gqset)
  call delete_Gk_vectors(Gqbarc)
  
  return
end subroutine

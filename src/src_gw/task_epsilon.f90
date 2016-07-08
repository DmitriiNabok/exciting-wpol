      
subroutine task_epsilon(iflag)

    use modinput
    use modmain,               only : zzero, evalsv, efermi
    use modgw
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
            
    implicit none
    ! input
    integer, intent(in) :: iflag
    ! local
    integer(4) :: iq, fid, iom
    integer    :: im

    character(80) :: fname
    
    !===========================================================================
    ! Initialization
    !===========================================================================

    ! prepare GW global data
    call init_gw
    
    ! clean not used anymore global exciting variables
    call clean_gndstate
    
    if (.not.input%gw%rpmat) then
      !========================================================
      ! calculate momentum matrix elements and store to a file
      !========================================================
      call calcpmatgw
    end if
    
    ! occupancy dependent BZ integration weights
    call kintw
    
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
    
    ! each process does a subset
    do iq = iqstart, iqend
    
      write(*,*)
      write(*,*) '(task_gw): q-point cycle, iq = ', iq
    
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
      
      !========================================
      ! Set v-diagonal MB and reduce its size
      !========================================
      call setbarcev(input%gw%barecoul%barcevtol)
      
      !===================================
      ! Calculate the dielectric function
      !===================================
      call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
      
      call calcepsilon(iq,iomstart,iomend)

      ! ! Calculate the screened Coulomb potential W_{ij}
      if (iflag < 0) then
        call calcinveps(iomstart,iomend)
        ! W = v^{1/2} (e^{-1}-1) v^{1/2}
        do iom = iomstart, iomend
        do im = 1, mbsiz
          epsilon(im,im,iom) = barcev(im)*epsilon(im,im,iom)
        end do
        end do
      end if

      ! store q-dependent Wij
      call getunit(fid)
      write(fname,'("WMAT-mat-q",I4.4,".OUT")') iq
      open(fid, File=trim(fname), Action='WRITE')
      do iom = 1, freq%nomeg, 100
      do im = 1, mbsiz
        write(fid,'(i8, 2f16.6)') im, epsilon(im,im,iom)
      end do
      write(fid,*)
      end do
      close(fid)

      write(fname,'("WMAT-iom-q",I4.4,".OUT")') iq
      open(fid, File=trim(fname), Action='WRITE')
      do im = mbsiz-100, mbsiz, 10
      do iom = iomstart, iomend
        write(fid,'(i8, 2f16.6)') iom, epsilon(im,im,iom)
      end do
      write(fid,*)
      end do
      close(fid)

      call delete_coulomb_potential      
      call delete_dielectric_function(Gamma)
      if (allocated(kcw)) deallocate(kcw)
      if (allocated(unw)) deallocate(unw)
      
      ! clean unused data
      if (allocated(mpwipw)) deallocate(mpwipw)
      if (allocated(barc)) deallocate(barc)
      
    end do ! iq
    
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
end subroutine


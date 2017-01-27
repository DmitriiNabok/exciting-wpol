!BOP
!
!!ROUTINE: calcbarcmb
!
!!INTERFACE:
!
subroutine calcbarcmb(iq)
!
!!DESCRIPTION:
!
!This subroutine calculates the matrix of the bare coulomb potential
!
!!USES:
    use modinput
    use modgw
    use mod_mpi_gw, only: myrank
    use mod_wpol_diagonalization
    
!!INPUT PARAMETERS: 
    implicit none
    integer, intent(in) :: iq ! index of the q-point

!!LOCAL VARIABLES:
    integer :: imix, jmix, igq, jgq
    real(8) :: tstart, tend, t0, t1
      
!!REVISION HISTORY:
! 
! Created Jan 2014 by DIN
!
!EOP
!BOC      
    call timesec(tstart)
        
!=============================================================================== 
! Setup the bare Coulomb potential matrix in MB representation 
!===============================================================================

    if (allocated(barc)) deallocate(barc)
    allocate(barc(matsiz,matsiz))
    barc(:,:) = 0.d0
    
    select case (trim(input%gw%barecoul%basis))
    
    case('pw')
    
      call calcmpwmix(iq)
      call calcbarcmb_pw(iq)
      
    case('mb')
    
      if (Gamma) then
        !------------------------------------------------
        ! Matrix elements for the singular q=0, L=0 case
        !------------------------------------------------
        !call timesec(t0)
        call barcq0
        !call timesec(t1)
        !write(*,*) 'barcq0', t1-t0 
      end if
      
      !if (vccut) call calcmpwmix(iq)
        
      !-----------------------------------------------------------
      ! Matrix elements between MT and MT mixed product functions
      !-----------------------------------------------------------
      !call timesec(t0)
      call calcbarcmb_mt_mt(iq)
      !call timesec(t1)
      !write(*,*) 'calcbarcmb_mt_mt', t1-t0 

      !-----------------------------------------------------------
      ! Matrix elements between an atomic mixed function and an IPW
      !-----------------------------------------------------------
      !call timesec(t0)
      call calcbarcmb_ipw_mt(iq)
      !call timesec(t1)
      !write(*,*) 'calcbarcmb_ipw_mt', t1-t0
    
      !-----------------------------------------------------------
      ! Matrix elements between two IPW's
      !-----------------------------------------------------------
      !call timesec(t0)
      call calcbarcmb_ipw_ipw(iq)
      !call timesec(t1)
      !write(*,*) 'calcbarcmb_ipw_ipw', t1-t0

    case default
    
      write(*,*) 'ERROR(calcbarcmb): Unknown basis type!'
      stop
      
    end select

    if (input%gw%debug) then
      write(fdebug,*) "### barc ###"
      do imix = 1, matsiz, matsiz/10
        do jmix = 1, matsiz, matsiz/10
          write(fdebug,'(2i5,4e16.6)') imix, jmix, &
          &    barc(imix,jmix), barc(imix,jmix)-conjg(barc(jmix,imix)) 
        end do
      end do
    endif !debug    
    
!===============================================================================
! Diagonalize the bare coulomb matrix
!===============================================================================

    !call timesec(t0)

    if (allocated(vmat)) deallocate(vmat)
    allocate(vmat(matsiz,matsiz))
      
    if (allocated(barcev)) deallocate(barcev)
    allocate(barcev(matsiz))
    
if (.false.) then 
    vmat(1:matsiz,1:matsiz) = barc(1:matsiz,1:matsiz)
    deallocate(barc)  
    call mkl_zheev ( matsiz, vmat, barcev )
else
    vmat(1:matsiz,1:matsiz) = barc(1:matsiz,1:matsiz)
    deallocate(barc)
    call mkl_zheevr ( matsiz, vmat, barcev )
end if

    !call timesec(t1)
    !write(*,*) 'barc diagonalization', t1-t0

!----------------------    
! debug info
!----------------------
    if (input%gw%debug) then
      !----------------------    
      ! Memory usage info
      !----------------------
      msize = sizeof(barcev)*b2mb+sizeof(vmat)*b2mb
      write(*,'("calcbarcmb: rank, size(Coulomb potential) (Mb):",i4,f12.2)') myrank, msize
      write(fdebug,*) "### barcev ###"
      do imix = 1, matsiz
        write(fdebug,'(i5,e16.6)') imix, barcev(imix) 
      end do
      write(fdebug,*) "### vmat ###"
      do imix = 1, matsiz, matsiz/10
        do jmix = 1, matsiz, matsiz/10
          write(fdebug,'(2i5,2e16.6)') imix, jmix, vmat(imix,jmix)
        end do
      end do
    endif !debug

    call timesec(tend)
    time_barcmb = time_barcmb+tend-tstart
    
end subroutine ! calcbarcmb
!EOC

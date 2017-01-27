!BOP
!
!!ROUTINE: diagsgi
!!INTERFACE:
!
subroutine diagsgi(iq)
!
! !DESCRIPTION:
!
! This subroutine generates the overlap matrix of the product functions in the
! interstitial region and diagonalizes it. 
! The output is the matrix $S_{\vec{G}i}$.
!
!!USES:
    use modinput
    use modmain, only : cfunig
    use modgw,   only : Gset, Gqset, sgi, &
    &                   fdebug, time_diagsgi
    use mod_wpol_diagonalization
      
!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq      

!!LOCAL VARIABLES:
    integer(4) :: ig, ngq
    integer(4) :: igq   ! Counter: Runs over igq's
    integer(4) :: jgq   ! Counter: Runs over igq's
    integer(4), dimension(3) :: iig    ! integer coordinates of G-G'
    real(8)   :: tstart, tend
    
    complex(8) :: cfact
    real(8), allocatable :: epsipw(:)
!
! !REVISION HISTORY:
!
! Created Dec. 2003. by RGA
! Last Modification May. 2006 by RGA
! Revisited 10.05.2011 by DIN
!
!EOP
!BOC
    call timesec(tstart)
    
    ngq = Gqset%ngk(1,iq)

    if(allocated(sgi)) deallocate(sgi)
    allocate(sgi(ngq,ngq))
    
    !-----------------------------------------------------------
    ! Calculate the overlap matrix between product plane waves
    !-----------------------------------------------------------
    sgi = 0.0d0
    do igq = 1, ngq
      ig = Gset%ivgig(0,0,0)
      sgi(igq,igq) = conjg(cfunig(ig))
      ! Non diagonal elements:
      do jgq = igq+1, ngq
        iig(:) = Gset%ivg(:,Gqset%igkig(igq,1,iq))- &
        &        Gset%ivg(:,Gqset%igkig(jgq,1,iq))
        ig = Gset%ivgig(iig(1),iig(2),iig(3))
        ! if (ig>??) stop 'G basis is too small for the selected G_{max}^{MB}'
        sgi(igq,jgq) = conjg(cfunig(ig))
        sgi(jgq,igq) = conjg(sgi(igq,jgq))
      end do ! jgq
    end do ! igq
      
    if (input%gw%debug) then
      write(fdebug,*) "diagsgi: iq, ngq=", iq, ngq
      write(fdebug,*) "### sgi-0 ###"
      do jgq = 1, ngq, ngq/10
        do igq = 1, ngq, ngq/10
          write(fdebug,'(2i5,4f12.6)') igq, jgq, sgi(igq,jgq)
        end do
      end do
    end if 

    !-----------------------------------------------------------
    ! Diagonalize sgi
    !-----------------------------------------------------------
    allocate(epsipw(ngq))
    
if (.false.) then
    call mkl_zheev ( ngq, sgi, epsipw )
else
    call mkl_zheevr ( ngq, sgi, epsipw )
end if

    if (input%gw%debug) then
      write(fdebug,*) "### sgi-1 ###"
      do jgq = 1, ngq, ngq/10
        do igq = 1, ngq, ngq/10
          write(fdebug,'(2i5,4f12.6)') igq, jgq, sgi(igq,jgq)
        end do
      end do
      write(fdebug,*)
      write(fdebug,*) "### epsipw ###"
      do igq = 1, ngq
        write(fdebug,'(i5,2f12.6)') igq, epsipw(igq)
      end do 
    end if
    
    !-----------------------------------------------------------
    ! Normalize sgi
    !-----------------------------------------------------------
    do igq = 1, ngq
      cfact = cmplx(1.0d0/sqrt(dabs(epsipw(igq))),0.0d0,8)
      sgi(:,igq) = cfact*sgi(:,igq)
    enddo
    deallocate(epsipw)

    call timesec(tend)
    time_diagsgi = time_diagsgi+tend-tstart
      
    return
end subroutine
!EOC


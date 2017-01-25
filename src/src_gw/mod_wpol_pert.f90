

module mod_wpol_pert
 
  use mod_wpol_diagonalization
   
  implicit none
 
contains

  ! uses perturbation theory to calculate the eigenvalues and eigenvectors of the D^2+DMMD matrix 
  subroutine wpol_pert(nvck,d,dmmd,tvck,cutoff)

    implicit none
 
    integer,    intent(in)    :: nvck
    real(8),    intent(in)    :: d(nvck),cutoff
    complex(8), intent(inout) :: dmmd(nvck,nvck)
    real(8),    intent(out)   :: tvck(nvck)   
    !local
    integer :: degsub(nvck),ndeg
    integer, allocatable :: nsub(:),tmp(:)
	real(8) :: maxdmmd
    
    
    !the cutoff to determine the subsets
    write(*,*) 'Cutoff: ',cutoff
    
    !maximal amount of degeneracies is none
    allocate(nsub(nvck))
    
    ! degsub contains all degenerate subsets
    ! nsub contains how many elements each subset has
    call find_degeneracy(nvck,d,degsub,nsub,ndeg,cutoff)
    write(*,*) 'subsets: ',ndeg
    
    allocate(tmp(ndeg))
    tmp(:)=nsub(1:ndeg)
    deallocate(nsub)
    allocate(nsub(ndeg))
    call move_alloc(tmp,nsub) 

    ! apply perturbation theory and diagonalize the subsets if necessary
    call apply_pert_theory_1st(nvck,degsub,nsub,ndeg,d,dmmd,tvck,maxdmmd)
    
	!write convergence to file
	open(73,file='pert_test.dat',form='FORMATTED',status='UNKNOWN',action='WRITE',position='append')
	write(73,*) 'Cutoff: ', cutoff
	write(73,*) 'Maximum dmmd/DeltaE: ',maxdmmd
	write(73,*) 'Subsets: ',ndeg
	close(73)
	
    deallocate(nsub)
 
  end subroutine

  
!--------------------------------------------------------------------------------------------------- 

  ! find sets of degenerate eigenvalues
  subroutine find_degeneracy(nvck,d,degsub,nsub,ndeg,cutoff)
    implicit none
    integer, intent(in)  :: nvck
    real(8), intent(in)  :: d(nvck),cutoff
    integer, intent(out) :: degsub(nvck),nsub(nvck),ndeg
    ! local
    integer :: i,n
    real(8) :: delta       
    
    !sort d  
    call sortidx(nvck,d,degsub) 
    
    nsub = 0 
    nsub(1)=1
    ndeg = 1
        
    !create degsub, indices are grouped for each subset
    !nsub contains the number of elements for each subset (nsub(1)=number of elements in first subset)
    !ndeg contains the number of actual degenerate subsets
    do i = 2, nvck
     
     delta=abs(d(degsub(i-1))**2-d(degsub(i))**2)

     if (delta<cutoff) then
      !add a new element to the subset
      nsub(ndeg)=nsub(ndeg)+1
     else
      !create a new subset
      ndeg=ndeg+1
      nsub(ndeg)=1
     end if
    
    end do
    
  end subroutine   
    


 
!---------------------------------------------------------------------------------------------------   
  
  ! applies 1st order degenerate perturbation theory
  subroutine apply_pert_theory_1st(nvck,degsub,nsub,ndeg,d,dmmd,tvck,maxdmmd)
    implicit none
    integer,    intent(in)    :: nvck,ndeg
    integer,    intent(in)    :: degsub(nvck),nsub(nvck)
    real(8),    intent(in)    :: d(nvck)
    complex(8), intent(inout) :: dmmd(nvck,nvck)
    real(8),    intent(out)   :: tvck(nvck)
    real(8),    intent(out)   :: maxdmmd
	
    !local
    real(8), allocatable :: alpha(:)
    real(8) :: delta
    complex(8) :: tmp
    complex(8), allocatable :: beta(:,:)
    !complex(8) :: U(nvck,nvck) 
    integer :: x,nmax
    integer :: i,j,m,n,l,p

    x=0
	maxdmmd=0.d0
    
    do i=1,ndeg
    !degsub(sum_i nsub(i-1):nsub(i)) is the current subset to consider
    nmax=nsub(i)
    if (nmax/=1) then !apply deg. pert. theory for degenerate subset 
     
     allocate(beta(nmax,nmax))
     allocate(alpha(nmax)) 
     
     !assamble degenerate subset matrix
     do m = 1, nmax
      do n = 1, nmax
       beta(m,n) = dmmd( degsub(x+m), degsub(x+n) )
      end do
     end do 
     
     !diagonalize      
     call mkl_zheev(nmax,beta,alpha)
     
     ! energies
     do m=1,nmax
      tvck(degsub(x+m)) = d(degsub(x+m))*d(degsub(x+m))+alpha(m)
     end do

     ! vectors 1st order
     do m=1,nmax
      l=degsub(x+m)
      !do 1..x,x+nsub(ndeg)..nvck
      do j=1,x
       p=degsub(j)
	   tmp=dmmd(p,l)*beta(m,m)
	   dmmd(p,l)=0.d0
       do n = 1, nmax
        !U(p,l) = U(p,l) + dmmd(p,degsub(x+n))*beta(n,m)
		dmmd(p,l) = dmmd(p,l) + dmmd(p,degsub(x+n))*beta(n,m)
       end do
	   dmmd(p,l)=dmmd(p,l)+tmp
       delta = (d(l)-d(p))*d(l)-d(p)*(d(p)-d(l)) 
       dmmd(p,l) = dmmd(p,l) /delta
	   if (abs(dmmd(p,l))>maxdmmd) then
	    maxdmmd=dmmd(p,l)
	    !write(*,*) 'Warning: ',p,l,abs(U(p,l))
       end if
      end do
      do j=x+nmax+1,nvck
       p=degsub(j)
	   tmp=dmmd(p,l)*beta(m,m)
	   dmmd(p,l)=0.d0	   
       do n = 1, nmax
        !U(p,l) = U(p,l) + dmmd(p,degsub(x+n))*beta(n,m)
		dmmd(p,l) = dmmd(p,l) + dmmd(p,degsub(x+n))*beta(n,m)
       end do
	   dmmd(p,l)=dmmd(p,l)+tmp
       delta = (d(l)-d(p))*d(l)-d(p)*(d(p)-d(l)) 
       dmmd(p,l) = dmmd(p,l) /delta
	   if (abs(dmmd(p,l))>maxdmmd) then
	    maxdmmd=dmmd(p,l)
	    !write(*,*) 'Warning: ',p,l,abs(U(p,l))
       end if
      end do
      end do
     
     ! vectors 0th order
     do m = 1, nmax
      do n=1, nmax
       !U(degsub(x+m),degsub(x+n)) = U(degsub(x+m),degsub(x+n)) + beta(m,n)
	   dmmd(degsub(x+m),degsub(x+n)) =  beta(m,n)
      end do
     end do
	 
     deallocate(beta,alpha)
    
    else !apply normal pert. theory for non-degenerate subset 
     
     !get the index for the non-degenerate subset
     l=degsub(x+1)
     !energies
     tvck(l) = d(l)*d(l)+real(dmmd(l,l))     
     ! vectors
     do j = 1, nvck
      if (l /= j) then
       delta =  (d(l)-d(j))*d(l)-d(j)*(d(j)-d(l)) 
       dmmd(j,l) =  dmmd(j,l)/delta
	   if (abs(dmmd(j,l))>maxdmmd) then
	    maxdmmd=dmmd(j,l)
       end if
      else
       dmmd(j,l) = 1.d0
      end if
     end do
    
    end if !
    !count x
    x=x+nmax
    end do !ndeg
    

  end subroutine   
  
  
!---------------------------------------------------------------------------------------------------   

end module







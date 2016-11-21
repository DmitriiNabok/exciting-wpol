

module mod_wpol_pert
 
  use mod_wpol_diagonalization
   
  implicit none
 
contains

  ! uses perturbation theory to calculate the eigenvalues and eigenvectors of the D^2+DMMD matrix 
  subroutine wpol_pert(nvck,d,dmmd,tvck)

    implicit none
 
    integer,    intent(in)    :: nvck
    real(8),    intent(in)    :: d(nvck)
    complex(8), intent(inout) :: dmmd(nvck,nvck)
    real(8),    intent(out)   :: tvck(nvck)   
    !local
    integer, allocatable :: degsub(:,:), tmp(:,:)
    integer :: nsub,ndeg,m,n
    real(8) :: maxdmmd  
    
    maxdmmd=maxval(abs(dmmd))
    
    ! degsub(:,:) contains all degenerate subsets
    allocate(degsub(nvck,nvck))
    call find_degeneracy(nvck,d,degsub,nsub,ndeg,maxdmmd)
    
    !reallocate degsub to max. amount of degeneracy found
    allocate(tmp(nsub,ndeg))
    do m= 1, nsub
    do n = 1, ndeg
      tmp(m,n) = degsub(m,n)
    end do
    end do
    deallocate(degsub)

    allocate(degsub(nsub,ndeg))
    call move_alloc(tmp,degsub)
    
    write(*,*) 'Subsets: ',nsub
    write(*,*) 'Max amount of degeneracy: ',ndeg

    ! apply perturbation theory and diagonalize the subsets if necessary
    call apply_pert_theory_1st(nvck,degsub,nsub,ndeg,d,dmmd,tvck)
    
    deallocate(degsub)
 
  end subroutine

  
!--------------------------------------------------------------------------------------------------- 

  ! find sets of degenerate eigenvalues
  subroutine find_degeneracy(nvck,d,degsub,nsub,ndeg,maxdmmd)
    implicit none
    integer, intent(in)  :: nvck
    real(8), intent(in)  :: d(nvck),maxdmmd
    integer, intent(out) :: degsub(nvck,nvck)
    integer, intent(out) :: nsub,ndeg
    ! local
    integer :: i,m,n
    real(8) :: delta,cutoff,first,last
    integer :: idx(nvck)
        
    
    !sort d  
    call sortidx(nvck,d,idx) 
    
    degsub = 0
    ndeg = 1
    nsub = 1 
    m=1
    n=1
    
    !the cutoff to determine the subsets
    cutoff=maxdmmd
    !cutoff=0.06d0
    write(*,*) 'Cutoff: ',cutoff
    
    degsub(m,n)=idx(1)
    
    !create degsub, with indices of the degernerate subsets
    do i = 1, nvck-1
     
     delta=abs(d(idx(i))**2-d(idx(i+1))**2)

     if (delta<cutoff) then
      
      n=n+1
      if (n>ndeg) ndeg=n
      !cutoff=cutoff-delta
      !last = d(idx(i+1))**2
     
     else
     
      n=1
      m=m+1
      if (m>nsub) nsub=m
      !first = d(idx(i+1))**2
      !if (abs(first-last)>cutoff) write(*,*) 'the case', abs(first-last)
      !cutoff=maxdmmd
      
     end if
    
     degsub(m,n)=idx(i+1)
    
    end do
    
    
  end subroutine   
    


 
!---------------------------------------------------------------------------------------------------   
  
  ! applies 1st order degenerate perturbation theory
  subroutine apply_pert_theory_1st(nvck,degsub,nsub,ndeg,d,dmmd,tvck)
    implicit none
    integer,    intent(in)    :: nvck
    integer,    intent(in)    :: degsub(nsub,ndeg)
    real(8),    intent(in)    :: d(nvck)
    integer,    intent(in)    :: nsub,ndeg
    complex(8), intent(inout) :: dmmd(nvck,nvck)
    real(8),    intent(out)   :: tvck(nvck)
    
    !local
    real(8),    allocatable :: alpha(:),delta
    complex(8), allocatable :: beta(:,:) !, ap(:)
    complex(8) :: U(nvck,nvck) 
    integer :: i,j,n,m,l,kappa,nmax

    
    U = 0.d0
    kappa=0
    
    do m = 1, nsub

     !count number of degeneracies
     nmax = 0
     do n = 1, ndeg
      if (degsub(m,n)/=0) nmax = nmax+1
     end do
     write(*,*) 'subset dimension: ',nmax
     
     kappa=kappa+nmax
     
     !degeneracy -------------------------------------------------------------------------- 
     if (nmax>1) then
     
      !allocate(ap(nmax*(nmax+1)/2))
      allocate(beta(nmax,nmax))
      allocate(alpha(nmax)) 
      
      ! set up matrix for each degenerate subset in packed storage
      do j = 1, nmax
       do i = 1, j
        beta(i,j) = dmmd( degsub(m,i), degsub(m,j) )
        !ap(i+j*(j-1)/2) = dmmd( degsub(m,i), degsub(m,j) )
       end do
      end do 
      
      ! diagonalize
      call mkl_zheev(nmax,beta,alpha)
      !call mkl_zhpev(nmax,ap,beta,alpha)
      
      ! apply perturbation theory for degeneracy
      do i = 1, nmax
       
       l = degsub(m,i)
        
       ! energies
       tvck(l) = d(l)*d(l)+alpha(i)
       
       ! vectors 0th order
       do j = 1, nmax
        U(degsub(m,j),l) = U(degsub(m,j),l) + beta(j,i)
       end do
       
       ! vectors 1st order
       do j = 1, nvck
        if ( ALL(degsub(m,:) /= j) ) then 
        
         do n = 1, nmax
          U(j,l) = U(j,l) + dmmd(j,degsub(m,n))*beta(n,i)
         end do
         
         delta = (d(l)-d(j))*d(l)-d(j)*(d(j)-d(l)) 
         
         U(j,l) = U(j,l) /delta
         
         if (abs(U(j,l))>1.d0) then
          write(*,*) 'warning (deg): ',j,l,d(l),d(j),U(j,l)
         end if
        
        end if
       end do
               
      end do

      deallocate(beta,alpha)
      !deallocate(ap)

     else
     !no degeneracies, apply normal perturbation theory------------------------------------------
     
     
      l = degsub(m,1)
     
      ! energies
      if (aimag(dmmd(l,l)) > 1.d-8) then
       write(*,*) 'Warning: imaginary part of DMMD(i,i) not vanishing'
      end if
      
      tvck(l) = d(l)*d(l)+real(dmmd(l,l))     
     
      ! vectors
      do j = 1, nvck
       if (l /= j) then
          
        delta =  (d(l)-d(j))*d(l)-d(j)*(d(j)-d(l)) 
          
        U(j,l) = dmmd(j,l)/delta
        
       else
        U(j,l) = 1.d0
       end if
      end do 
      
     end if
     
    end do
     
    write(*,*) kappa
    dmmd(:,:)=U(:,:)

  end subroutine   
  
  
!---------------------------------------------------------------------------------------------------   

end module







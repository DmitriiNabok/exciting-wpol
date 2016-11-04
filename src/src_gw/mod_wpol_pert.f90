
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
    integer(8) :: bla
    complex(8) :: U(nvck,nvck)
    
    ! degsub(:,:) contains all degenerate subsets
    allocate(degsub(nvck,nvck))
    call find_degeneracy(nvck,d,dmmd,degsub,ndeg,nsub)
    
    !reallocate degsub to max. amount of degeneracy found
    allocate(tmp(nsub,ndeg))
    do n = 1, ndeg
     do m = 1, nsub
      tmp(m,n) = degsub(m,n)
     end do
    end do
    deallocate(degsub)

    allocate(degsub(nsub,ndeg))
    call move_alloc(tmp,degsub)
    
    write(*,*) 'Subsets: ',nsub
    write(*,*) 'Max amount of degeneracy: ',ndeg
!      do m=1,nsub   
!      if (ANY(degsub(m,:)/=0)) then
!      !write(*,*) degsub(m,:)
!      bla=bla+1
!      end if
!      end do
!      write(*,*) bla
!     bla=0
!     do m=1,nsub
!     do n=1,ndeg
!      if (degsub(m,n)/=0) then
!      bla=bla+1
!      end if
!     end do
!     end do
!     write(*,*) bla
    
    ! apply perturbation theory and diagonalize the subsets if necessary
    call apply_pert_theory_1st(nvck,degsub,nsub,ndeg,d,dmmd,U,tvck)
    !call apply_pert_theory_2nd(nvck,degsub,nsub,ndeg,d,dmmd,U,tvck) !this is not implemented correctly
    
    deallocate(degsub)
    
    dmmd = U
 
  end subroutine

  
!--------------------------------------------------------------------------------------------------- 

  ! find sets of degenerate eigenvalues
  subroutine find_degeneracy(nvck,d,dmmd,degsub,ndeg,nsub)
    implicit none
    integer, intent(in)  :: nvck
    real(8), intent(in)  :: d(nvck)
    integer, intent(out) :: degsub(nvck,nvck)
    integer, intent(out) :: nsub,ndeg
    complex(8), intent(in) :: dmmd(nvck,nvck)
    ! local
    integer :: i,j,n,m,a,nmax,iwork(nvck,nvck)
    real(8) :: delta

    !     
    iwork = 0
    nmax=1
    
    !first iteration
    do i = 1, nvck
    
     n = 1
        
       iwork(i,1) = i 
        
       do j = 1, nvck
        if (i/=j) then
         
         delta = (d(i)-d(j))*d(i) - d(j)*(d(j)-d(i))
         
         if ( (abs(delta) < 10*abs(dmmd(i,j))) ) then 
         
          !degeneracies in subset counter
          n = n+1
          if (n > nmax) then
           nmax = n
          end if
          
          !degeneracy in i,j
          !write into degeneracy array
          !iwork(i,1) = i
          iwork(i,n) = j
          
         end if
        
        end if       
       end do

    end do
    
    
    !
    degsub=0
    nsub = 1
    ndeg = 1
    m=1
    
    !second iteration
    do i = 1, nvck
    
     if (ALL(degsub/=i)) then
     
     if (m > nsub) then
      nsub = m
     end if 
     
     !length of non-zero elements in iwork(i,:)
     a=0
     do n=1,nmax
      if (iwork(i,n)/=0) then
       degsub(m,n)=iwork(i,n)
       a=a+1
      end if
     end do     
     !for each element in degsub(m,:) check all entries in corresponding iwork
     !if there are elements not in degsub(m,:) append them to degsub
     do n = 1, nvck
      do j = 1, nvck
         
       if ( (degsub(m,n)/=0) .and. &
          & (ALL(degsub(m,:)/=iwork(degsub(m,n),j))) .and. &
          & (iwork(degsub(m,n),j)/=0) ) then
         
        degsub(m,a+1)=iwork(degsub(m,n),j)
          
        a=a+1
        if (a > ndeg) then
         ndeg = a
        end if
         
       end if
         
      end do
     end do
     
     !number of subsets counter
     m = m+1         
     
     end if

    end do

  end subroutine

 
!---------------------------------------------------------------------------------------------------   
  
  ! applies 1st order degenerate perturbation theory
  subroutine apply_pert_theory_1st(nvck,degsub,nsub,ndeg,d,dmmd,U,tvck)
    implicit none
    integer,    intent(in)    :: nvck
    integer,    intent(in)    :: degsub(nsub,ndeg)
    real(8),    intent(in)    :: d(nvck)
    integer,    intent(in)    :: nsub,ndeg
    complex(8), intent(in)    :: dmmd(nvck,nvck)
    real(8),    intent(out)   :: tvck(nvck)
    complex(8), intent(out)   :: U(nvck,nvck)
    
    !local
    real(8),    allocatable :: alpha(:),delta
    complex(8), allocatable :: beta(:,:)
    integer :: i,j,n,m,l,kappa,nmax
    
    U = 0.d0
    kappa=0
    
    do i = 1, nsub
     
     
     !count number of degeneracies
     nmax = 0
     do m = 1, ndeg
      if (degsub(i,m)/=0) nmax = nmax+1
     end do
     write(*,*) 'subset dimension: ',nmax
     
     kappa=kappa+nmax
     
     ! no degeneracy, apply normal perturbation theory
     if (nmax==1) then
     
      l = degsub(i,1)
     
      ! energies
      if (aimag(dmmd(i,i)) > 1.d-8) then
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
        
        
     !degeneracy ------------------------------------------------------------
     else
           
      allocate(beta(nmax,nmax))
      allocate(alpha(nmax)) 
      
      ! set up matrix for each degenerate subset
      do j = 1, nmax
       do m = 1, nmax
        beta(j,m) = dmmd( degsub(i,j), degsub(i,m) )
       end do
      end do 
      
      ! diagonalize
      call mkl_zheev(nmax,beta,alpha)
      
      ! apply perturbation theory for degeneracy
      do m = 1, nmax
       
       l = degsub(i,m)
        
       ! energies
       tvck(l) = d(l)*d(l)+alpha(m)
       
       ! vectors
       do j = 1, nmax
        U(degsub(i,j),l) = U(degsub(i,j),l) + beta(j,m)
       end do
       
       do j = 1, nvck
        if ( ALL(degsub(i,:) /= j) ) then 
        
         do n = 1, nmax
          U(j,l) = U(j,l) + dmmd(j,degsub(i,n))*beta(n,m)
         end do

         delta =  (d(l)-d(j))*d(l)-d(j)*(d(j)-d(l)) 
         
         U(j,l) = U(j,l)/delta
         
        if (abs(U(j,l))>1.d0) then
         write(*,*) 'warning (deg): ',j,l,d(l),d(j),abs(U(j,l))
        end if
        
        end if
       end do
               
      end do

      deallocate(beta,alpha)

     end if
     
    end do
    
    write(*,*) kappa

  end subroutine   
  
  
!---------------------------------------------------------------------------------------------------   
!---------------------------------------------------------------------------------------------------   
  
  ! applies 2nd order degenerate perturbation theory
  subroutine apply_pert_theory_2nd(nvck,degsub,nsub,ndeg,d,dmmd,U,tvck)
    implicit none
    integer,    intent(in)    :: nvck
    integer,    intent(in)    :: degsub(nsub,ndeg)
    real(8),    intent(in)    :: d(nvck)
    integer,    intent(in)    :: nsub,ndeg
    complex(8), intent(in)    :: dmmd(nvck,nvck)
    real(8),    intent(out)   :: tvck(nvck)
    complex(8), intent(out)   :: U(nvck,nvck)
    
    !local
    real(8),    allocatable :: alpha(:),delta,delta2
    complex(8), allocatable :: beta(:,:)
    complex(8) :: pert
    integer :: i,j,n,m,l,kappa,nmax,p,q
    
    U = 0.d0
    kappa=0
    
    do i = 1, nsub
     
     
     !count number of degeneracies
     nmax = 0
     do m = 1, ndeg
      if (degsub(i,m)/=0) nmax = nmax+1
     end do
     write(*,*) 'subset dimension: ',nmax
     
     kappa=kappa+nmax
     
     ! no degeneracy, apply normal perturbation theory
     if (nmax==1) then
     
      l = degsub(i,1)
     
      ! energies
      if (aimag(dmmd(i,i)) > 1.d-8) then
       write(*,*) 'Warning: imaginary part of DMMD(i,i) not vanishing'
      end if
      
      tvck(l) = d(l)*d(l)+real(dmmd(l,l))     
      
      
      do j = 1, nvck
       if (l /= j) then
        
        delta =  (d(l)-d(j))*d(l)-d(j)*(d(j)-d(l)) 
        
        !energies 2nd order
        tvck(l) = tvck(l)+abs(dmmd(j,l))**2/delta
        
        !first order vectors  and 2nd order negative part 
        pert = dmmd(j,l)*dmmd(l,l)/delta**2
        
        if (abs(pert) < 1.d0) then
        U(j,l) = dmmd(j,l)/delta - pert
        else
        U(j,l) = dmmd(j,l)/delta 
        end if
        
        !2nd order vectors (positive part)
        do m=1,nvck
         if (l /= m) then
          
          delta2 = (d(l)-d(m))*d(l)-d(m)*(d(m)-d(l)) 
          
          pert = dmmd(j,m)*dmmd(m,l)/(delta*delta2)
          
          if (abs(pert)<1.d0 ) then
          U(j,l) = U(j,l) + dmmd(j,m)*dmmd(m,l)/(delta*delta2)
          end if
          
         end if     
        end do

       else
       
        U(j,l) = 1.d0
        
       end if
      end do 
        
        
     !degeneracy ------------------------------------------------------------
     else
           
      allocate(beta(nmax,nmax))
      allocate(alpha(nmax)) 
      
      ! set up matrix for each degenerate subset
      do j = 1, nmax
       do m = 1, nmax
        beta(j,m) = dmmd( degsub(i,j), degsub(i,m) )
       end do
      end do 
      
      ! diagonalize
      call mkl_zheev(nmax,beta,alpha)
      
      ! apply perturbation theory for degeneracy
      do m = 1, nmax
       
       l = degsub(i,m)
        
       ! energies 1st order
       tvck(l) = d(l)*d(l)+alpha(m)
       
       ! vectors 0th order
       do j = 1, nmax
        U(degsub(i,j),l) = U(degsub(i,j),l) + beta(j,m)
       end do
       
       do j = 1, nvck
        if ( ALL(degsub(i,:) /= j) ) then 
         
         do n = 1, nmax
          U(j,l) = U(j,l) + dmmd(j,degsub(i,n))*beta(n,m)
         end do
         
         delta =  (d(l)-d(j))*d(l)-d(j)*(d(j)-d(l)) 
         
         !energies 2nd order
         tvck(l) = tvck(l) + abs(U(j,l))**2/delta
         
         ! vector 1st order
         U(j,l) = U(j,l)/delta 
         
         pert=0.d0
         ! vector 2nd order (negative)
         do p=1,nmax
          do q=1,nmax
           !< n|dmmd|n>
           pert =pert + U(j,l)*CONJG(beta(p,m))*dmmd(degsub(i,p),degsub(i,q))*beta(q,m)/delta
          end do
         end do
         
         if (abs(pert)<1.d0 ) then
         
         U(j,l) = U(j,l) - pert
         
         end if
         
         ! vector 2nd order (pos)
         do p=1,nvck
          if ( ALL(degsub(i,:) /= p) ) then 
          
          delta2 = (d(l)-d(p))*d(l)-d(p)*(d(p)-d(l)) 
          
          pert = 0.d0
          do q=1,nmax
           pert = pert + dmmd(j,p)*dmmd(p,degsub(i,q))*beta(q,m)  
          end do
          
          pert = pert/(delta*delta2)
          
          if (abs(pert)<1.d0 ) then
          
          U(j,l) = U(j,l) + pert   
              
          end if
          end if
         end do
         
         
        if (abs(U(j,l))>1.d0) then
         write(*,*) 'no wonder (deg)',j,l,d(l),d(j),abs(U(j,l))
        end if
        
        end if
       end do
               
      end do

      deallocate(beta,alpha)

     end if
     
    end do
    
    write(*,*) kappa

  end subroutine   
  
  
!---------------------------------------------------------------------------------------------------   

end module







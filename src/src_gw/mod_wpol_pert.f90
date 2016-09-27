
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
    integer, allocatable :: iwork(:,:), tmp(:,:)
    integer :: nmax,m
    real(8) :: deg_cutoff
    complex(8) :: U(nvck,nvck)
    logical :: pert_cond

    pert_cond = .true.
    deg_cutoff = 0.05d0
 
    do while (pert_cond)
    
      ! iwork(i,:) contains all degenerate subsets for d(i)*d(i) 
      allocate(iwork(nvck,nvck))
      call find_degeneracy(nvck,d,iwork,nmax,deg_cutoff)
    
      ! reallocate iwork to max. amount of degeneracy found (nmax)
      allocate(tmp(nvck,nmax))
      do m = 1, nmax
        tmp(:,m) = iwork(:,m)
      end do
      deallocate(iwork)

      allocate(iwork(nvck,nmax))
      call move_alloc(tmp,iwork)
    
      ! apply perturbation theory and diagonalize the subsets if necessary
      call apply_pert_theory_1st(nvck,iwork,nmax,d,dmmd,U,tvck,pert_cond)
    
      ! increase the cutoff, in case the necessary perturbation condition |dmmd(j,i)/(E_i-E_j)| < 1 is not met
      ! pert_cond will be set to true then again in apply_pert_theory
      deg_cutoff = deg_cutoff+0.1d0 
      deallocate(iwork)    
    
    end do
    
    dmmd = U
 
  end subroutine

  
!--------------------------------------------------------------------------------------------------- 

  ! find sets of degenerate eigenvalues
  subroutine find_degeneracy(nvck,d,iwork,nmax,deg_cutoff)
    implicit none
    integer, intent(in)  :: nvck
    real(8), intent(in)  :: d(nvck)
    integer, intent(out) :: iwork(nvck,nvck)
    integer, intent(out) :: nmax
    real(8), intent(in)  :: deg_cutoff
    ! local
    integer :: i,j,n
    real(8) :: t1
         
    iwork = 0
    nmax = 1
    
    do i = 1, nvck
      if ( ALL(iwork/=i) ) then
        n = 2
        do j = 1, nvck
          t1 = (d(i)-d(j))*d(i) - d(j)*(d(j)-d(i))
          if ( (abs(t1) < deg_cutoff) .and. &
          &    (i /= j)               .and. &
          &    (ALL(iwork/=j)) )      then 
            iwork(i,1) = i
            iwork(i,n) = j
            n = n+1
            if (n > nmax) then
              nmax = n
            end if
          end if
        end do
      end if 
    end do
  
  end subroutine

 
!---------------------------------------------------------------------------------------------------   
  
  ! applies 1st order degenerate perturbation theory
  subroutine apply_pert_theory_1st(nvck,iwork,nmax,d,dmmd,U,tvck,pert_cond)
    implicit none
    integer,    intent(in)    :: nvck
    integer,    intent(in)    :: iwork(nvck,nmax)
    real(8),    intent(in)    :: d(nvck)
    integer,    intent(in)    :: nmax
    complex(8), intent(in)    :: dmmd(nvck,nvck)
    real(8),    intent(out)   :: tvck(nvck)
    complex(8), intent(out)   :: U(nvck,nvck)
    logical,    intent(inout) :: pert_cond
    
    !local
    real(8),    allocatable :: alpha(:)
    complex(8), allocatable :: degsub(:,:)
    integer :: i,j,n,m,l,kappa
    
    U = 0.d0
    
    kappa = 0
    
    pert_cond = .false. 
    
    do i = 1, nvck
     
      ! no degeneracy, apply normal perturbation theory
      if ( ALL(iwork /= i) ) then
      
        write(*,*) 'no-degeneracy: ', i
      
        ! energies
        if (aimag(dmmd(i,i)) > 1.d-8) then
          write(*,*) 'Warning: imaginary part of DMMD(i,i) not vanishing'
        end if
        tvck(i) = d(i)*d(i)+real(dmmd(i,i))
      
        ! vectors
        do j = 1, nvck
          if (i /= j) then
            U(j,i) = dmmd(j,i)/( (d(i)-d(j))*d(i)-d(j)*(d(j)-d(i)) )
            ! check for necessary perturbation condition
            if (.not. (abs(U(j,i))<1.d0) ) then
              write(*,*) 'Warning: Degeneracy parameter too small! (non-degenerate)'
              pert_cond = .true. 
            end if
          else
            U(j,i) = 1.d0
          end if
        end do 
    
      !degeneracy ------------------------------------------------------------
      else if ( ANY(iwork(i,:) /= 0) ) then
     
        ! find number of degenerate values
        n = 0
        do m = 1, nmax
          if (iwork(i,m)/=0) n = n+1
        end do
        write(*,*) 'subset dimension: ',n
      
        kappa = kappa+n
      
        allocate(degsub(n,n))
        allocate(alpha(n)) 
      
        ! set up matrix for each degenerate subset
        do j = 1, n
          do m = 1, n
            degsub(j,m) = dmmd( iwork(i,j), iwork(i,m) )
          end do
        end do 
      
        ! diagonalize
        call mkl_zheev(n,degsub,alpha)
        !write(*,*) 'degsub: ',degsub
        !write(*,*) 'alpha: ', alpha
      
        ! apply perturbation theory for degeneracy
        do m = 1, n

          ! energies
          tvck(iwork(i,m)) = d(iwork(i,m))*d(iwork(i,m))+alpha(m)
       
          ! vectors
          do l = 1, n
            U(iwork(i,l),iwork(i,m)) = U(iwork(i,l),iwork(i,m)) + degsub(l,m)
          end do
       
          do j = 1, nvck
            if ( ALL(iwork(i,:) /= j) ) then 
              do l = 1, n
                U(j,iwork(i,m)) = U(j,iwork(i,m)) + dmmd(j,iwork(i,l))*degsub(l,m)
              end do
              U(j,iwork(i,m)) = U(j,iwork(i,m))/( (d(iwork(i,m))-d(j))*d(iwork(i,m))-d(j)*(d(j)-d(iwork(i,m))) )
         
              ! check for necessary perturbation condition
              if (.not.(abs(U(j,iwork(i,m)))<1.d0 )) then
                write(*,*) 'Warning: Degeneracy parameter too small! (in subset)'
                pert_cond = .true. 
              end if 
            end if
          end do
               
        end do

        deallocate(degsub,alpha)
     
      end if
    
    end do
    
    write(*,*) 'kappa: ', kappa

  end subroutine   
  
  
!---------------------------------------------------------------------------------------------------   
  
  !applies 2nd order degenerate perturbation theory - not implemented yet
  subroutine apply_pert_theory_2nd(nvck,iwork,nmax,d,dmmd,U,tvck,pert_cond)
  
    implicit none
    
    integer, intent(in) :: nvck
    integer, intent(in) :: iwork(nvck,nmax)
    real(8), intent(in) :: d(nvck)
    integer, intent(in) :: nmax
    complex(8), intent(in) :: dmmd(nvck,nvck)
    real(8), intent(out) :: tvck(nvck)
    complex(8),intent(out) :: U(nvck,nvck)
    logical,intent(inout) :: pert_cond 
    
    !local
    real(8),allocatable :: alpha(:)
    complex(8),allocatable :: degsub(:,:)
    integer :: i,j,n,m,l,kappa
    
    U=0.d0
    
    kappa=0
    
    pert_cond=.false. 
    
    do i=1,nvck
     
     !no degeneracy, apply normal perturbation theory -----------------------
     if (ALL(iwork/=i) ) then
      
      write(*,*) 'no-degeneracy: ',i
      
      !energies
      if (aimag(dmmd(i,i))>1.d-8) then
       write(*,*) 'Warning: imaginary part of DMMD(i,i) not vanishing'
      end if
      tvck(i)=d(i)*d(i)+real(dmmd(i,i))
      
      !vectors
      do j=1,nvck
       if (i/=j) then
        U(j,i) = dmmd(j,i)/( (d(i)-d(j))*d(i)-d(j)*(d(j)-d(i)) )
        
        !check for necessary perturbation condition
        if (.not.(abs(U(j,i))<1.d0 )) then
         write(*,*) 'Warning: Degeneracy parameter too small! (non-degenerate)'
         pert_cond=.true. 
        end if         
        
!         if (.not.(abs(dmmd(j,i)/( (d(i)-d(j))*d(i)-d(j)*(d(j)-d(i)) ))<1.d0 )) then
!          write(*,*) 'non-deg: dmmd, d_i^2-d_j^2', dmmd(j,i), d(i)*d(i)-d(j)*d(j)
!          write(*,*) 'Warning: Degeneracy parameter too small! (not in subset)'
!          pert_cond=.true. 
!         end if 
        
       else
        U(j,i)=1.d0
       end if
      end do 
    
    
     !degeneracy ------------------------------------------------------------
     else if (ANY(iwork(i,:)/=0)) then
     
      !find number of degenerate values
      n=0
      do m=1,nmax
       if (iwork(i,m)/=0) then
        n=n+1
       end if
      end do
      
      write(*,*) 'subset dimension: ',n
      
      kappa=kappa+n
      
      allocate(degsub(n,n))
      allocate(alpha(n)) 
      
      !set up matrix for each degenerate subset
      do j=1,n
       do m=1,n
        degsub(j,m) = dmmd( iwork(i,j), iwork(i,m) )
       end do
      end do 
      
      !diagonalize
      call mkl_zheev(n,degsub,alpha)
      
       !write(*,*) 'degsub: ',degsub
       !write(*,*) 'alpha: ', alpha
      
      !apply perturbation theory for degeneracy
      do m=1,n

       !energies
       tvck(iwork(i,m)) = d(iwork(i,m))*d(iwork(i,m))+alpha(m)
       
       !vectors
       do l=1,n
        U(iwork(i,l),iwork(i,m)) = U(iwork(i,l),iwork(i,m)) + degsub(l,m)      
       end do
       
       do j=1,nvck
        if (ALL(iwork(i,:)/=j)) then 
                
         do l=1,n
          U(j,iwork(i,m)) = U(j,iwork(i,m)) + dmmd(j,iwork(i,l))*degsub(l,m)
         end do
         
         
         U(j,iwork(i,m)) = U(j,iwork(i,m))/( (d(iwork(i,m))-d(j))*d(iwork(i,m))-d(j)*(d(j)-d(iwork(i,m))) )
         
         !check for necessary perturbation condition
         if (.not.(abs(U(j,iwork(i,m)))<1.d0 )) then
          write(*,*) 'Warning: Degeneracy parameter too small! (in subset)'
          pert_cond=.true. 
         end if 
         
!          if (.not.(abs(U(j,iwork(i,m))/( (d(iwork(i,m))-d(j))*d(iwork(i,m))-d(j)*(d(j)-d(iwork(i,m))) ))<1.d0 )) then
!           write(*,*) 'deg: dmmd, d_i^2-d_j^2', U(j,iwork(i,m)), ( (d(iwork(i,m))-d(j))*d(iwork(i,m))-d(j)*(d(j)-d(iwork(i,m))) )
!           write(*,*) 'Warning: Degeneracy parameter too small! (subset)'
!           pert_cond=.true. 
!          end if 
    
        end if
        
       end do
               
      end do

      deallocate(degsub,alpha)
     
     end if
    end do 
    
    
    write(*,*) 'kappa: ', kappa

  end subroutine   

end module







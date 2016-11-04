
module mod_rank1_update

  use modmain, only : zzero, zone
  implicit none
  real(8), parameter :: eps = 1.d-8

contains

  subroutine test_rank1_update()
    use mod_wpol_diagonalization, only : mkl_zheev
    implicit none

    integer, parameter :: N = 4
    real(8)    :: d(N), rho 
    complex(8) :: z(N), H(N,N), Q(N,N)
    integer :: i, j

    write(*,*)
    write(*,*) '--------------------------'
    write(*,*) 'Rank-1 update'
    write(*,*) '--------------------------'

    rho = 1.d0
    d(:) = (/1.d0, 2.d0, 3.d0, 4.d0/)                 ! diagonal matrix
    z(:) = (/(1.d0,1.d0), (1.d0,1.d0), zzero, zzero/) ! column vector
    z(:) = z(:) / sqrt(sum(z*conjg(z)))

    write(*,*)
    write(*,*) 'rho = ', rho
    write(*,*) 'd = ', d
    write(*,*) 'z =', z

    ! H = D + rho*z*z'
    H(:,:) = zzero
    do i = 1, N
      H(i,i) = d(i)
      do j = i, N
        H(i,j) = H(i,j) + rho*z(i)*conjg(z(j))
        H(j,i) = conjg(H(i,j))
      end do
    end do
    write(*,*)    
    write(*,*) 'H ='
    do i = 1, N
      write(*,'(4f12.6)') dble(H(i,:))
    end do
    write(*,*)

    ! Rank-1 update
    Q(:,:) = zzero
    do i = 1, N
      Q(i,i) = zone
    end do
    call sevr(N, d, z, rho, Q)

    write(*,*) 'eigenvalue:'
    write(*,*) d
    write(*,*) 'eigenvectors (real):'
    do i = 1, N
      write(*,'(8f12.6)') Q(i,:)
    end do 

    !--------------------------
    ! Direct diagonalization
    !--------------------------
    write(*,*)
    write(*,*) '--------------------------'
    write(*,*) 'Direct diagonalization'
    write(*,*) '--------------------------'
    call mkl_zheev(N, H, d)
    write(*,*) 'eigenvalue:'
    write(*,*) d
    write(*,*) 'eigenvectors (real):'
    do i = 1, N
      write(*,'(8f12.6)') H(i,:)
    end do 

    return
  end subroutine

!==============================================================================
! Symmetric eigenproblem rank-one modification
!      Q*diag(d)*Q' = Qin*(diag(din) + rho*z*z')*Qin'
!
! inputs
!   din -- input eigenvalues
!     z -- updating vector, normalized in 2-norm
!   rho -- scalar factor in updating
!   Qin -- input unitary matrix, optional, if not provided
!          the outputs are the same as if Qin = eye
! outputs
!     d -- updated eigenvalues
!     Q -- updated unitary (eigenvector) matrix
! ifail -- 0 for success, 1 for fail
!
! dependency
!  deflate -- introduce zeros in z for equal diagonal elements in din
!      slv -- secular equation solver for eigenvalues and eigenvectors
!
! Reference:
!  W. Xu and S. Qiao.
!  A divide-and-conquer method for the Takagi factorization.
!  Technical Report No. CAS 05-01-SQ. Department of Computing and
!  Software, McMaster University, Hamilton, Ontario L8S 4K1, Canada.
!  February 2005.
!
! W. Xu and S. Qiao   McMaster Univ.   April 2005

  subroutine sevr(n, d, z, rho, Q)
    implicit none
    integer,    intent(in)     :: n
    real(8),    intent(inOut)  :: d(n)
    complex(8), intent(inOut)  :: z(n)
    real(8),    intent(in)     :: rho
    complex(8), intent(inOut)  :: Q(n,n)
    ! local
    integer :: k, i
    integer,    allocatable :: inzeros(:), indx(:)
    real(8),    allocatable :: dd(:), lam(:)
    complex(8), allocatable :: zz(:), G(:,:)

    real(8) :: t1, t2

    ! Matrix deflation
    allocate(inzeros(n), indx(n))
    call timesec(t1)
    call deflate(n, d, z, rho, Q, k, inzeros, indx)
    call timesec(t2)
    write(*,*) "Time(deflate): ", t2-t1
    allocate(dd(k),zz(k))
    do i = 1, k
      dd(i) = d(inzeros(i)) ! entries are distinct and sorted in decreasing order
      zz(i) = z(inzeros(i)) ! nonzero entries
    end do

    ! compute eigenvalues and eigenvectors
    allocate(lam(k))
    lam(:) = 0.d0
    allocate(G(k,k))
    G(:,:) = zzero

! #ifdef USEOMP
! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
! !$OMP DO
! #endif
    call timesec(t1)
    do i = 1, k
      call slv(i, k, dd, zz, abs(rho), lam(i), G(:,i))
    end do
    call timesec(t2)
    write(*,*) "Time(slv): ", t2-t1
! #ifdef USEOMP
! !$OMP END DO
! !$OMP END PARALLEL
! #endif

    deallocate(dd, zz)

    ! update eigenvalues
    do i = 1, k
      d(inzeros(i)) = lam(i)
    end do
    deallocate(lam)
    d(:) = sign(1.d0,rho)*d(:)

    ! update eigenvectors
    Q(:,indx(inzeros(1:k))) = matmul( Q(:,indx(inzeros(1:k))), G)
    Q(:,:) = Q(:,indx)

    deallocate(inzeros, indx)
    return
  end subroutine

!==============================================================================
! Deflation stage in the symmetric eigenvalue rank-one
! modification problem: diag(d) + rho*z*z'.
! Sort sign(rho)*d and introduce zeros into z for equal 
! entries in d
!
! Inputs
!       d -- eigenvalue vector to be updated
!       z -- symmetric rank-one update vector, normalized
!     rho -- scalar factor in update
!       Q -- unitary matrix
! Outputs
!       d -- eigenvalue vector after deflation
!            d sorted in decreasing order
!       z -- symmetric rank-one update vector after deflation
! inzeros -- indices to the nonzero entries of z
!            d(inzeros) entries are distinct
!            z(inzeros) entries are nonzeros
!    indx -- permutation vector for Q
!       Q -- updated unitary matrix
!
! On return
!
! Q(:indx)*sign(rho)*(sign(rho)*diag(d) + abs(rho)*z*z')*Q(:,indx)'
! =
! Qin*(diag(din) + rho*zin*zin')*Qin'
!
! Reference:
!  W. Xu and S. Qiao.
!  A divide-and-conquer method for the Takagi factorization.
!  Technical Report No. CAS 05-01-SQ. Department of Computing and
!  Software, McMaster University, Hamilton, Ontario L8S 4K1, Canada.
!  February 2005.
!
! W. Xu and S. Qiao  McMaster Univ.   April 2005
!
  subroutine deflate(n, d, z, rho, Q, k, inzeros, indx)
    use modmain, only : zzero, zone
    implicit none
    integer,    intent(in)     :: n
    real(8),    intent(inOut)  :: d(n)
    complex(8), intent(inOut)  :: z(n)
    real(8),    intent(in)     :: rho
    complex(8), intent(inOut)  :: Q(n,n)
    integer,    intent(out)    :: k
    integer,    intent(out)    :: inzeros(n)
    integer,    intent(out)    :: indx(n)
    ! local
    integer :: j
    real(8) :: signrho, told, normz, dmax, tolz
    real(8) :: cc
    complex(8) :: r, c, s, M(2,2)
    real(8),    allocatable :: d1(:)
    complex(8), allocatable :: z1(:)

    inzeros(:) = 0  ! indices of nonzeros in z
    
    ! sort signrho*din in decreasing order
    signrho = sign(1.d0,rho)
    call sortidx(n, -signrho*d, indx)

    allocate(d1(n))
    d1(:) = d(indx)
    ! write(*,*) 'd1=', d1

    ! permute z correspondingly
    allocate(z1(n))
    z1(:) = z(indx)
    ! write(*,*) 'z1=', z1

    dmax = max( abs(d1(1)), abs(d1(n)) )
    
    ! tolz = eps
    tolz = (n*eps*dmax) / abs(rho)  ! tolerance for numerical zeros in z

    k = 0
    j = 1
    do while ( j < n )
      if ( abs(z1(j)) > tolz ) then  ! z(j) not very small
        normz = sqrt( z1(j)*conjg(z1(j)) + z1(j+1)*conjg(z1(j+1)) )
        ! told = eps
        told = (dmax*eps)/(abs(z1(j))/normz);
        ! tolerance for d(j) = d(j+1) numerically
        if ( abs(d1(j) - d1(j+1)) < told ) then  ! d(j) ~ d(j+1)
          ! find a rotation to eliminate z(j) using z(j+1)
          call rotate(z1(j+1), z1(j), c, s)
          z1(j) = zzero
          z1(j+1) = normz
          ! update Q
          M(1,:) = (/  conjg(c), s /)
          M(2,:) = (/ -conjg(s), c /)
          Q(:,indx(j:j+1)) = matmul(Q(:,indx(j:j+1)), M)
        else  ! d1(j) does not equal d1(j+1)
          k = k+1
          inzeros(k) = j
        end if
      end if
      j = j+1
    end do

    if (abs(z1(n)) > tolz) then  ! check the last z(n)
      k = k+1
      inzeros(k) = n
    end if

    d(:) = d1(:)
    z(:) = z1(:)
    deallocate(d1, z1)

    return
  end subroutine

!==============================================================================
!
! Generate cosine and sine of Givens rotation
! such that
!  ( c' s )(a) = (*)
!  (-s' c )(b)   (0)
! INPUT:   a, b  complex numbers
! OUTPUT:  c, s  cosine and sine of the
!          rotation which eliminates b
!
! Reference:
! G.Golub and C.Van Loan, Matrix Computations,
! 2ed Edition, The John Hopkins University Press,
! 1989, p.202.
!
! Implemented by S. Qiao, McMaster University
!
  subroutine rotate(a, b, c, s)
    implicit none
    complex(8), intent(in)  :: a, b
    complex(8), intent(out) :: c, s
    ! local
    complex(8) :: t1, t2

    if (abs(b) < eps) then
      if (abs(a) < eps) then
        c = zone
      else
        c = a / abs(a)
      end if
      s = zzero
      return
    end if

    if (abs(a) < eps) then
      c = zzero
      s = b / abs(b)
      return
    end if

    if ( abs(a) > abs(b) ) then
      t1 = b / a
      t2 = sqrt(zone + t1*conjg(t1))
      c  = a / abs(a) / t2
      s  = t1 * c
    else
      t1 = a / b
      t2 = sqrt(zone + t1*conjg(t1))
      s  = b / abs(b) / t2
      c  = t1 * s
    end if

    return
  end subroutine

!==============================================================================
!
! This function finds the root lam of the rational function
!          1 + rho*sum(|z(1:n)|^2/(d(1:n) - lam))
! in the interval
!      (d(i),d(i-1)) for i>1  or  (d(1),d(1) + rho) for i=1.
! It is used in the rank-one modification of symmetric eigenproblem:
!         diag(d) + rho*z*z'
! So, this function computes the i-th eigenvalue lam of the rank-one
! modified matrix and its corresponding eigenvector.
!
! inputs:
!   i -- interval (d(i),d(i-1)) where the eigenvalue to be computed
!   d -- the original eigenvalues, the entries are sorted
!        in decreasing order
!   z -- the updating vector, normalized (2-norm)
! rho -- a positive scalar factor in the updating
!
! outputs:
!   lam -- ith eigenvalue after updating
!     g -- the corresponding eigenvector
!
! References
! J.R. Bunch, C.P. Nielsen and D.C. Sorensen,
! Rank-one modification of the symmetric eigenproblem,
! Numer. Math. 31, 31-48 (1978).
!
! W. Xu and S. Qiao, McMaster Univ. Nov. 2004
! revised  April 2005

  subroutine slv(i, n, d, z, rho, lam, g)
    implicit none
    integer,    intent(in)  :: i
    integer,    intent(in)  :: n
    real(8),    intent(in)  :: d(n)
    complex(8), intent(in)  :: z(n)
    real(8),    intent(in)  :: rho
    real(8),    intent(out) :: lam
    complex(8), intent(out) :: g(n)
    ! local
    integer :: j, maxiter, niter, im1, ifail
    real(8) :: del, lambda, eps1, t
    real(8), allocatable :: delta(:), z2(:)
    real(8) :: a, b, c
    complex(8) :: w, temp
    complex(8) :: psi, phi, dpsi, dphi
    real(8), external :: dznrm2

    maxiter = 1000 ! maximum number of iterations
    ifail   = 1    ! exit status
    im1     = i-1    

    ! |z|^2
    allocate(z2(n))
    do j = 1, n
      z2(j) = z(j)*conjg(z(j))
    end do

    allocate(delta(n))
    del = d(i)
    delta(:) = (d(:) - del) / rho

    ! calculate an initial guess
    if ( i > 1 ) then
      del = d(im1)
    else
      del =  d(1) + rho;
    end if

    a = 0.d0
    do j = i+1, n
      a = a + z2(j) / (d(j) - del)
    end do
    a = rho*a
    ! write(*,*) 'a=', a

    b = 0.d0
    do j = 1, i-2
      b = b + z2(j) / (d(j) - del)
    end do
    b = rho*b
    ! write(*,*) 'b=', b

    a = a + b + 1.d0
    
    if ( i == 1 ) then
      t = z2(i) / a
    else
      t = a*delta(im1)
      b = t + z2(i) + z2(im1)
      if ( abs(b) > eps ) then
        t = 2.d0*z2(i)*delta(im1)/(b + sqrt(abs(b*b - 4.d0*t*z2(i))))
      else
        t = (b - sqrt(b*b - 4.d0*t*z2(i)))/(2.d0*a)
      end if
    end if

    ! Iterative solution
    lambda = 0.d0
    niter  = 0
    do while (niter < maxiter)
      niter = niter + 1
      if ( i > 1 ) then
        ! t is too close to endpoint
        if (t > 0.9*delta(im1)) t = 0.9*delta(im1) ! back off a little
      end if
      delta(:) = delta(:) - t
      lambda = lambda + t
      lam = d(i) + rho*lambda  ! adjust the root
    
      ! evaluate psi and its derivative dpsi
      psi = 0.d0
      dpsi = 0.d0
      do j = i, n        
        t = z2(j) / delta(j)
        psi = psi + t
        dpsi = dpsi + t/delta(j)
      end do

      ! evaluate phi and its derivative dphi
      phi = 0
      dphi = 0
      if ( i > 1 ) then
        do j = 1, im1
          t = z2(j) / delta(j)
          phi = phi + t
          dphi = dphi + t / delta(j)
        end do
      end if
    
      ! test for convergence
      w = 1.d0 + phi + psi
      ! eps1 = eps
      eps1 = eps*n*(abs(phi)+abs(psi)+1.d0)
      if (abs(w) < eps1) then    ! converged
        g(:) = z(:) / delta      ! eigenvector
        ! t = sqrt(sum(g*conjg(g)))
        t = dznrm2(n, g, 1)
        g(:) = g(:) / t
        ifail = 0                ! success
        deallocate(z2, delta)
        return
      end if
    
      ! calculate the new estimate
      if (i == 1) then
        t = (w*psi) / dpsi
      else
        del = delta(im1)
        temp = psi / dpsi  
        a = 1.d0 + phi - del*dphi
        if (abs(a) < eps) return ! fail            
        b = (del*(1.d0+phi) + psi*temp) / a + temp
        c = (2.d0*temp*del*w) / a
        t = c / (b + sqrt(abs(b*b - 2.d0*c)))
      end if

    end do ! while

    if (ifail /= 0) then
      write(*,*) 'WARNING(mod_rank1_update::slv) Failed to reach convergence!'
      ! stop
    end if
 
     deallocate(z2,delta)
    return
  end subroutine

!------------------------------------------------------------------------------
! SVD of a rectangular complex matrix A
!
! input:
!       -- rectangular matrix (complex)    : A(m,n)
!
! output:
!       -- right singular vectors (complex): A(m,n) = VT(m,n)
!       -- singular eigenvalues (real)     : S(min(m,n))
!
  subroutine generate_update_vectors(m,n,A,S)
    implicit none
    integer,    intent(in)    :: m
    integer,    intent(in)    :: n
    complex(8), intent(inOut) :: A(m,n)
    real(8),    intent(Out)   :: S(min(m,n))
    ! local
    integer :: lwork, lrwork, info
    real(8),    allocatable :: rwork(:)
    complex(8), allocatable :: U(:,:), VT(:,:), work(:)
    external zgesvd

    lrwork = 5*min(m,n)
    allocate(rwork(lrwork))
    !
    ! Query the optimal workspace
    !
    lwork = -1
    allocate(work(1),U(1,1),VT(1,1))
    call zgesvd( 'N', 'O', m, n, A, m, S, U, m, VT, n, &
    &             work, lwork, rwork, info )
    lwork  = int(work(1))
    deallocate(work)
    !
    ! Compute SVD
    !
    allocate(work(lwork))
    call zgesvd( 'N', 'O', m, n, A, m, S, U, m, VT, n, &
    &             work, lwork, rwork, info )

    deallocate(work,rwork)
    deallocate(U,VT)
    !
    ! check for convergence.
    !
    if( info > 0 ) then
      write(*,*)'The algorithm computing svd failed to converge.'
      stop
    end if

    return
  end subroutine

end module
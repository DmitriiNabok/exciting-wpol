
module mod_lanczos

  implicit none


contains

!------------------------------------------------------------------------------
  subroutine lanczos_simple(n,H,niter,eval,evec)
    use modmain, only : zzero, zone
    implicit none
    ! input
    integer,    intent(in) :: n
    complex(8), intent(in) :: H(n,n)
    integer,    intent(in) :: niter
    ! output
    real(8),    intent(out) :: eval(niter)
    complex(8), intent(out) :: evec(n,niter)
    ! parameters
    real(8), parameter :: eps = 1.d-8
    ! local
    integer :: i
    real(8),    allocatable :: a(:), b(:)
    complex(8), allocatable :: zr(:)
    complex(8), allocatable :: Q(:,:), evec1(:,:)

    ! external subroutines
    real(8),    external :: dznrm2
    complex(8), external :: zdotc
    external             :: zgemv

    allocate(a(niter),b(niter))
    a(:) = 0.d0
    b(:) = 0.d0

    allocate(Q(n,niter))
    ! starting vector = unity
    Q(:,:) = zzero
    Q(1,1) = zone

    allocate(zr(n))
    do i = 1, niter

      ! r = H*Q(i)
      call zgemv('n', n, n, zone, H, n, Q(:,i), 1, zzero, zr, 1)

      ! r = r - T(i,i-1)*Q(i-1)
      if (i > 1) zr(:) = zr(:) - b(i-1)*Q(:,i-1) 
      
      ! a = r'*Q(i)
      a(i) = dble( zdotc(n, zr, 1, Q(:,i), 1) )

      ! r = r - T(i,i)*Q(i)
      zr(:) = zr(:) - a(i)*Q(:,i) 

      if (i+1 > niter) exit ! stop iteration loop

      b(i) = dznrm2(n, zr, 1)
      if (b(i) > eps) then 
        Q(:,i+1) = zr(:) / b(i)
      else
        write(*,*) 'Warning(mod_lanczos::lanczos_simple): Zero basis vector! Quit execution...'
        exit ! stop iteration loop
      end if

    end do
    deallocate(zr)

    allocate(evec1(niter,niter))
    call diagonalize_tridiagonal(niter, a, b, eval, evec1)
    call zgemm( 'n', 'n', n, niter, niter, &
    &           zone, Q, n, evec1, niter,  &
    &           zzero, evec, n)
    deallocate(evec1)
    deallocate(a,b,Q)

    return
  end subroutine

!------------------------------------------------------------------------------
  subroutine lanczos_band(n,H,niter,blks,Q0,eval)
    use modmain, only : zone, zzero
    implicit none
    ! input / output
    integer,    intent(in)    :: n
    complex(8), intent(inout) :: H(n,n)
    integer,    intent(in)    :: niter
    integer,    intent(in)    :: blks
    complex(8), intent(in)    :: Q0(n,blks)
    real(8),    intent(out)   :: eval(blks*niter)
    ! parameters
    real(8), parameter :: eps = 1.d-8
    ! local
    integer :: i, j, m, iab
    real(8) :: norm
    complex(8), allocatable :: zr(:), evec1(:,:)
    complex(8), allocatable :: T(:,:), Q(:,:), AB(:,:)
    ! external subroutines
    real(8),    external :: dznrm2
    complex(8), external :: zdotc
    external             :: zgemv

    ! Lanczos basis size
    m = blks*niter

    if (m > N) then
      write(*,*) 'ERROR(mod_lanczos::lanczos_band): Inconsistent parameters!'
      write(*,*) '    (block size)*(number of iteration) > (matrix size)'
      stop
    end if

    allocate(zr(n))
    allocate(T(m,m))
    T(:,:) = 0.d0

    ! starting basis vectors
    allocate(Q(n,m))
    Q(:,1:blks) = Q0(:,1:blks)

    ! Block Lanczos iterations
    do i = 1, m

      ! r = H*q_i
      call zgemv('n', n, n, zone, H, n, Q(:,i), 1, zzero, zr, 1)

      ! r = r - \sum T_ij*q_j
      do j = max(1,i-blks), i-1
        zr(:) = zr(:) - T(j,i)*Q(:,j)
      end do

      do j = i, min(m,i+blks-1)
        ! T_ij = r^{+}*q_j
        T(i,j) = zdotc(n, zr, 1, Q(:,j), 1)
        ! r = r - T_{ij}*q_j
        zr(:) = zr(:) - conjg(T(i,j))*Q(:,j)
      end do

      if (blks+i > m) cycle

      !-----------------------------
      ! Prepare for next iteration
      !-----------------------------

      ! T_{i,i+n} = norm(r)
      norm = dznrm2(n, zr, 1)
      T(i,blks+i) = norm

      ! q_{i+n} = r / T(i,i+n)
      if (norm > eps) then
        Q(:,blks+i) = zr / norm
      else
        write(*,*) 'WARNING(mod_lanczos::lanczos_band): Zero basis vector! Quit execution...'
        exit
        ! one should do something here ...
        ! reduce block size?
      end if

    end do ! iteration loop
    deallocate(zr)

    ! Lapack 'U' banded matrix format
    allocate(AB(blks+1,m))
    AB(:,:) = 0.d0
    do j = 1, m
      do i = max(1,j-blks), j
        iab = blks+1+i-j
        AB(iab,j) = T(i,j)
      end do
    end do
    deallocate(T)

    ! Solve eigenvalue problem
    allocate(evec1(m,m))
    call diagonalize_banded(blks+1,m,AB,eval,evec1)
    deallocate(AB)

    ! Eigenvector of the original matrix
    ! evec H evec^{+} = (Q*T) eval*I (Q*T)^{+}
    !
    ! return it back to the input array
    call zgemm( 'n', 'n', n, m, m, &
    &           zone, Q, n, evec1, m,  &
    &           zzero, H(:,1:m), n)
    deallocate(Q)
    deallocate(evec1)

    return
  end subroutine


!------------------------------------------------------------------------------
  subroutine test_lanczos()
    use modmain, only : zzero, zone
    implicit none
    integer, parameter :: N = 4
    complex(8) :: A(N,N)
    real(8)    :: d(N), e(N-1)
    real(8)    :: d_ref(N), e_ref(N-1), eval_ref(N)
    complex(8) :: evec_ref(N,N)
    complex(8) :: Q(N,N), tau(N-1)
    integer    :: info, niter, blks, lwork
    integer    :: i, j, m, iab
    complex(8), allocatable :: work(:)
    real(8),    allocatable :: eval(:)
    complex(8), allocatable :: evec(:,:), evec1(:,:)
    complex(8), allocatable :: AB(:,:), T(:,:), V(:,:), Q0(:,:)

    real(8), external :: DNRM2, ZLANGE
    
    ! test matrix
    A(1,:) = (/  4.,  1., -2.,  2. /)
    A(2,:) = (/  1.,  2.,  0.,  1. /)
    A(3,:) = (/ -2.,  0.,  3., -2. /)
    A(4,:) = (/  2.,  1., -2., -1. /)

    ! tridiagonal elements
    d_ref(:) = (/ 4., 10./3., -33./25., 149./75. /)
    e_ref(:) = (/ 3., 5./3., 68./75. /)
    eval_ref(:) = (/ -2.19751697743942, 1.08436446377322, \
    &                 2.26853140643124, 6.84462110723497 /) 
    evec_ref(1,:) = (/ -0.1767051706244765,  -0.6422600067660329,  -0.2017110966209421,   0.7180459594506803 /)
    evec_ref(2,:) = (/ -0.1781004699063008,   0.5441878481946869,  -0.7894499124874400,   0.2211529881551809 /)
    evec_ref(3,:) = (/  0.2876680343197012,  -0.5202218509386151,  -0.5796341693167877,  -0.5573513771375095 /)
    evec_ref(4,:) = (/  0.9242849167461130,   0.1439822745761914,  -0.0102809986863304,   0.3533564748939931 /)

    write(*,*)
    write(*,*) '******************'
    write(*,*) '  Simple Lanczos'
    write(*,*) '******************'
    allocate(eval(N))
    allocate(evec(N,N))
    call lanczos_simple(N, A, N, eval, evec)
    do i = 1, N
      write(*,*) i, eval(i), eval_ref(i)
    end do
    write(*,*) '(lanczos_simple): diff(eval)=', &
    &           DNRM2(4, eval-eval_ref, 1)
    write(*,*) '(lanczos_simple): diff(evec)=', &
    &           ZLANGE( 'F', N, N, evec-evec_ref, N, work )
    deallocate(eval)
    deallocate(evec)

    !-----------------------
    ! 
    !-----------------------
    write(*,*)
    write(*,*) '******************'
    write(*,*) '   Band Lanczos'
    write(*,*) '******************'

    niter = 4
    blks = 1
    m = blks*niter

    allocate(eval(m))
    ! allocate(evec(N,m))
    allocate(Q0(N,blks))
    Q0(:,:) = zzero
    do i = 1, blks
      Q0(i,i) = zone
    end do
    ! call lanczos_band(N,A,niter,blks,Q0,eval,evec)
    call lanczos_band(N, A, niter, blks, Q0, eval)
    do i = 1, N
      write(*,*) i, eval(i), eval_ref(i)
    end do
    write(*,*) '(lanczos_band): diff(eval)=', &
    &           DNRM2(4, eval-eval_ref, 1)
    write(*,*) '(lanczos_band): diff(evec)=', &
    &           ZLANGE( 'F', N, N, A-evec_ref, N, work )
    do i = 1, N
      write(*,*) i
      write(*,'(4f12.6)') dble(evec(i,:))
      write(*,'(4f12.6)') dble(evec_ref(i,:))
    end do
    deallocate(eval)
    ! deallocate(evec)

    return
  end subroutine

!------------------------------------------------------------------------------
  subroutine diagonalize_banded(m,n,AB,eval,evec)
    use modmain, only : zzero, zone
    implicit none
    ! input
    integer,    intent(in) :: m, n
    complex(8), intent(in) :: AB(m,n)
    ! output
    real(8),    intent(out) :: eval(n)
    complex(8), intent(out) :: evec(n,n)
    ! local
    integer :: info
    real(8),    allocatable :: rwork(:)
    complex(8), allocatable :: work(:)

    external ZHBEV

    allocate(work(n))
    allocate(rwork(3*n-2))

    !-----------------------------------------
    ! Diagonalize the banded matrix
    !-----------------------------------------

    call ZHBEV( 'V',      & ! compute eigenvalues and eigenvectors
                'U',      & ! lower triangular part of matrix is stored in banded_matrix
                 n,       & ! dimension of matrix
                 m-1,     & ! number of superdiagonals in banded matrix
                 AB,      & ! matrix in banded storage
                 m,       & ! leading dimension of banded_matrix
                 eval,    & ! eigenvalues of matrix
                 evec,    & ! eigenvectors of matrix
                 n,       & ! dimension of eigenvector matrix
                 work, rwork, info )  ! work arrays and info

    if ( info /= 0) then
      write(*,*) 'ERROR(mod_lanczos::diagonalize_banded): Failed to diagonalize banded matrix!'
      stop
    end if

    return
  end subroutine


  subroutine diagonalize_tridiagonal(N, D, E, W, Z)
    implicit none
    integer,    intent(in)  :: N
    real(8),    intent(in)  :: D(N)
    real(8),    intent(in)  :: E(N)
    real(8),    intent(out) :: W(N)
    complex(8), intent(out) :: Z(N,N)
    
    ! local
    real(8) :: VL, VU, ABSTOL
    integer :: IL, IU, M, LWORK, LIWORK, INFO
    integer, allocatable :: ISUPPZ(:), IWORK(:)
    real(8), allocatable :: WORK(:)

    ABSTOL = 1.d-8

    allocate(ISUPPZ(2*N))

    LWORK = 18*N
    allocate(WORK(LWORK))

    LIWORK = 10*N
    allocate(IWORK(LIWORK))

    call zstegr  (  'V', &
                    'A', &
                     N,  &
                     D,  &
                     E,  &
                     VL, &
                     VU, &
                     IL, &
                     IU, &
                     ABSTOL, &
                     M, &
                     W, &
                     Z, &
                     N, &
                     ISUPPZ, &
                     WORK, &
                     LWORK, &
                     IWORK, &
                     LIWORK, &
                     INFO &
                  )

    if (info /= 0) then
      write(*,*) 'ERROR(mod_lanczos::diagonalize_tridiagonal): ZSTEGR'
      stop
    end if

  end subroutine

  subroutine orthogonalize(m,n,A,blks,Q)
    implicit none
    integer,    intent(in)  :: m
    integer,    intent(in)  :: n
    complex(8), intent(in)  :: A(m,n)
    integer,    intent(in)  :: blks
    complex(8), intent(out) :: Q(n,blks)
    ! local
    integer :: i
    integer :: lmn, lwork, lrwork, info
    real(8),    allocatable :: S(:), rwork(:)
    complex(8), allocatable :: A_(:,:), U(:,:), VT(:,:), work(:)
    external zgesvd

    allocate(A_(m,n))
    A_(:,:) = A(:,:)

    lmn = min(m,n)
    allocate(S(lmn))
    allocate(U(m,m),VT(n,n))
    lrwork = 5*lmn
    allocate(rwork(lrwork))
    !
    ! Query the optimal workspace
    !
    lwork = -1
    allocate(work(1))
    call zgesvd( 'all', 'all', m, n, A_, m, S, U, m, VT, n, &
    &             work, lwork, rwork, info )
    lwork  = int(work(1))
    deallocate(work)
    !
    ! Compute SVD
    !
    allocate(work(lwork))
    call zgesvd( 'all', 'all', m, n, A_, m, S, U, m, VT, n, &
    &             work, lwork, rwork, info )
    !
    ! check for convergence.
    !
    if( info > 0 ) then
      write(*,*)'The algorithm computing svd failed to converge.'
      stop
    end if

    do i = 1, blks
      Q(:,i) = conjg(VT(i,:))
    end do

    deallocate(A_)
    deallocate(S)
    deallocate(U,VT)
    deallocate(rwork,work)

    return
  end subroutine  

end
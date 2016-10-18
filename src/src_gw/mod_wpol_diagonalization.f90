
module mod_wpol_diagonalization

  implicit none

contains

!--------------------------------------------------------------------------------
  ! subroutine eig_rank1_update(n,D,M,m,evec,eval)
  !   implicit none
  !   integer,    intent(in)  :: n
  !   real(8),    intent(in)  :: D(n)
  !   complex(8), intent(in)  :: M(m,n)
  !   integer,    intent(in)  :: m
  !   complex(8), intent(out) :: evec(n,n)
  !   real(8),    intent(out) :: eval(n)
  !   ! local
  !   integer :: i, j, info
  !   complex(8), allocatable :: V(:,:)
  !   real(8),    allocatable :: S(:), D_old(:), D_new(:), W(:)
  !   real(8),    allocatable :: Q(:,:), Z(:,:)

  !   stol = 1.d-2

  !   ! step 1: M^{+}M matrix diagonalization and reduction
  !   allocate(S(n),D_old(n),D_new(n))
  !   allocate(V(n,n))
  !   call mkl_svd(n,M,m,V,S)

  !   allocate(W(n))
  !   allocate(Q(n,n),Z(n,n))

  !   k = 0
  !   do i = 1, n
  !     if (S(i) > stol) k = k+1
  !   end do

  !   k = 1
  !   do i = 1, n
  !     call dlaed9( k, k, k, n, D_new, Q, n, S(i), D_old, W,
  !     &            Z, n, info )
  !   end do



  ! end subroutine

!--------------------------------------------------------------------------------
  subroutine get_unique(n,A,etol,k,res,pos,mlt)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n)
    real(8), intent(in) :: etol
    integer, intent(out) :: k
    real(8), intent(out) :: res(n) 
    integer, intent(out) :: pos(n)
    integer, intent(out) :: mlt(n)
    ! local
    integer :: i, j
    integer, allocatable :: indx(:)

    ! sort array
    allocate(indx(n))
    call sortidx(n,A,indx)

    k = 1
    res(k) = A(indx(n))
    pos(k) = indx(n)
    mlt(:) = 1
    outer: do i = n-1, 1, -1
      do j = 1, k
        if ( abs( res(j) - A(indx(i)) ) < etol ) then
          ! Found a match so start looking again
          mlt(j) = mlt(j) + 1
          cycle outer
        end if
      end do
      ! No match found so add it to the output
      k = k + 1
      res(k) = A(indx(i))
      pos(k) = indx(i)
    end do outer

    deallocate(indx)
    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine mkl_chol(n,A,eval)
    use modmain, only : zzero, zone
    implicit none
    integer,    intent(in)     :: n
    complex(8), intent(inout)  :: A(n,n)
    real(8),    intent(out)    :: eval(n)

    ! local
    ! RWORK dimension should be at least MAX( 1, 5*MIN(M,N) )
    integer :: i, iter
    integer :: info
    complex(8), allocatable :: R(:,:), J(:,:)
    external zpotrf

    write(*,*) 'Info(mod_zmatrix_diagonalization::mkl_chol):'
    write(*,*) '   Diagonalization using Cholesky decomposition'
    write(*,*)

    allocate(R(n,n))
    allocate(J(n,n))

    ! 0-step
    J(:,:) = A(:,:)

    A(:,:) = 0.d0
    do i = 1, n
      A(i,i) = 1.d0
    end do

    do iter = 1, 4

      call zpotrf( 'U', n, J, n, info )
      if (info /= 0) then
        write(*,*)'The algorithm computing Cholesky decomposition failed.'
        stop
      end if
      R(:,:) = J(:,:)

      ! update eigenvector
      call zgemm( 'n', 'n', n, n, n, &
      &           zone, R, n, &
      &           A, n, &
      &           zzero, J, n)
      A(:,:) = J(:,:)

      ! prepare for next iteration: R*R^{+}
      call zgemm( 'n', 'c', n, n, n, &
      &           zone, R, n, &
      &           R, n, &
      &           zzero, J, n)

    end do

    ! eigenvalues
    do i = 1, n
      eval(i) = dble(J(i,i))
      write(77,*) i, eval(i)
    end do

    deallocate(J)
    deallocate(R)

    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine mkl_svd(n,A,m,evec,eval)
    implicit none
    integer,    intent(in)  :: n
    complex(8), intent(in)  :: A(m,n)
    integer,    intent(in)  :: m
    complex(8), intent(out) :: evec(n,n)
    real(8),    intent(out) :: eval(n)

    ! local
    ! RWORK dimension should be at least MAX( 1, 5*MIN(M,N) )
    integer :: i
    integer :: lmn, lwork, lrwork, info
    real(8),    allocatable :: S(:), rwork(:)
    complex(8), allocatable :: U(:,:), VT(:,:), work(:)
    complex(8), allocatable :: sigma(:,:), A_(:,:)
    external zgesvd
    real(8) :: dif_norm, a_norm

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
    !
    ! Check
    !
    ! allocate(sigma(m,n))
    ! sigma(:,:) = 0.d0
    ! do i = 1, m
    !   sigma(i,i) = cmplx(S(i),0.d0)
    ! end do
    ! A_(1:m,1:n) = matmul ( U(1:m,1:m), matmul ( sigma(1:m,1:n), VT(1:n,1:n) ) )
    ! a_norm = abs ( sqrt ( sum ( A(1:m,1:n)**2 ) ) )
    ! dif_norm = abs ( sqrt ( sum ( ( A(1:m,1:n) - A_(1:m,1:n) )**2 ) ) )
    ! write ( *, '(a)' ) ' '
    ! write ( *, '(a,g14.6)' ) '  Frobenius Norm of A, A_NORM = ', a_norm
    ! write ( *, '(a)' ) ' '
    ! write ( *, '(a)' ) '  ABSOLUTE ERROR for A = U*S*V'':'
    ! write ( *, '(a,g14.6)' ) '  Frobenius norm of difference A-U*S*V'' = ', dif_norm
    ! write ( *, '(a)' ) ' '
    ! write ( *, '(a)' ) '  RELATIVE ERROR for A = U*S*V'':'
    ! write ( *, '(a,g14.6)' ) '  Ratio of DIF_NORM / A_NORM = ', dif_norm / a_norm

    !
    ! Eigenvalues: \lambda = S^{+}*S
    !
    eval(:) = 0.d0
    do i = 1, m
      eval(i) = S(i)*S(i)
    end do
    !
    ! Eigenvectors
    !
    do i = 1, n
      evec(:,i) = conjg(VT(i,:))
    end do

    deallocate(A_)
    deallocate(S)
    deallocate(U,VT)
    deallocate(rwork,work)

    return
  end subroutine


!--------------------------------------------------------------------------------
  subroutine mkl_zheev(ndim,zevec,deval)
    implicit none
    integer,    intent(in)    :: ndim
    complex(8), intent(inout) :: zevec(ndim,ndim)
    real(8),    intent(out)   :: deval(ndim)
    ! local
    integer :: info, lwork, lrwork
    integer,    allocatable :: iwork(:)
    real(8),    allocatable :: rwork(:)
    complex(8), allocatable :: work(:)
    !
    write(*,*) 'Info(mod_zmatrix_diagonalization::mkl_zheev):'
    write(*,*) '   Diagonalization with ZHEEV'
    write(*,*)
    !
    ! Query the optimal workspace
    !
    lwork = -1
    allocate(work(1), rwork(1))
    call zheev( 'vectors', 'lower', ndim, zevec, ndim, deval, &
    &            work, lwork, rwork, info )
    lwork  = int(work(1))
    lrwork = 3*ndim-2
    deallocate(work, rwork)
    !
    ! solve eigenproblem
    !
    allocate(work(lwork), rwork(lrwork))
    call zheev( 'vectors', 'lower', ndim, zevec, ndim, deval, &
    &            work, lwork, rwork, info )
    deallocate(work, rwork)
    if (info > 0) then
      write(*,*)
      write(*,*) 'Error(mod_zmatrix_diagonalization::mkl_zheev):'
      write(*,*) '    ZHEEV algorithm failed to compute eigenvalues'
      write(*,*) '    info = ', info
      stop
    end if
    return
  end subroutine

!--------------------------------------------------------------------------------
  subroutine mkl_zheevr(ndim,zevec,deval)
    implicit none
    integer,    intent(in)    :: ndim
    complex(8), intent(inout) :: zevec(ndim,ndim)
    real(8),    intent(out)   :: deval(ndim)
    ! local
    integer :: info, lwork, lrwork, liwork, lwmax
    integer :: il, iu, m
    real(8) :: abstol, vl, vu
    integer,    allocatable :: iwork(:), isuppz(:)
    real(8),    allocatable :: rwork(:)
    complex(8), allocatable :: work(:), z(:,:)
    !
    write(*,*) 'Info(mod_zmatrix_diagonalization::mkl_zheevr):'
    write(*,*) '   Diagonalization with ZHEEVR'
    write(*,*)
    ! Negative ABSTOL means using the default value
    abstol = -1.0 
    !
    ! Query the optimal workspace.
    !
    lwork  = -1
    lrwork = -1
    liwork = -1
    allocate( z(ndim,ndim), isuppz(2*ndim) )
    allocate( work(1), rwork(1), iwork(1) )
    call zheevr( 'vectors', 'all', 'lower', ndim, zevec, ndim, &
    &            vl, vu, il, iu, abstol, m, deval, z, ndim, isuppz, &
    &            work, lwork, rwork, lrwork, iwork, liwork, &
    &            info )
    lwork  = int(work(1))
    lrwork = int(rwork(1))
    liwork = iwork(1)
    deallocate( work, rwork, iwork )
    !
    ! solve eigenproblem.
    !
    allocate( work(lwork), rwork(lrwork), iwork(liwork) )
    call zheevr( 'vectors', 'all', 'lower', ndim, zevec, ndim, &
    &            vl, vu, il, iu, abstol, m, deval, z, ndim, isuppz, &
    &            work, lwork, rwork, lrwork, iwork, liwork, &
    &            info )
    zevec(:,:) = z(:,:)
    deallocate(work, rwork, iwork)
    deallocate(z, isuppz)
    if (info > 0) then
      write(*,*)
      write(*,*) 'Error(mod_zmatrix_diagonalization::mkl_zheevr):'
      write(*,*) '    ZHEEVR algorithm failed to compute eigenvalues'
      write(*,*) '    info = ', info
      stop
    end if
    return
  end subroutine

end module

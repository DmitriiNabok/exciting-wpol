
module mod_wpol_diagonalization

  implicit none

contains

!--------------------------------------------------------------------------------
  subroutine mkl_svd(m,n,A,zevec,deval)
    use mod_wpol_tests
    implicit none
    integer,    intent(in)  :: m, n
    complex(8), intent(in)  :: A(m,n)
    complex(8), intent(out) :: zevec(n,n)
    real(8),    intent(out) :: deval(n)

    ! local
    ! RWORK dimension should be at least MAX( 1, 5*MIN(M,N) )
    integer :: i
    integer :: lmn, lwork, lrwork, info
    real(8),    allocatable :: S(:), rwork(:)
    complex(8), allocatable :: U(:,:), VT(:,:), work(:)
    complex(8), allocatable :: sigma(:,:), A_(:,:)
    external zgesvd
    real(8) :: dif_norm, a_norm

    write(*,*) 'Info(mod_zmatrix_diagonalization::mkl_svd):'
    write(*,*) '   Diagonalization using SVD decomposition'
    write(*,*)
    
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
    allocate(sigma(m,n))
    sigma(:,:) = 0.d0
    do i = 1, m
      sigma(i,i) = cmplx(S(i),0.d0)
    end do

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
    deval(:) = 0.d0
    do i = 1, m
      deval(i) = S(i)*S(i)
    end do

    !
    ! Eigenvectors
    !
    do i = 1, n
      zevec(:,i) = conjg(VT(i,:))
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

!--------------------------------------------------------------------------------
  subroutine lanczos(ndim,h,nreig,nsteps,deval,zevec)
    use modmain, only: zzero, zone
    implicit none
    integer,    intent(in)  :: ndim
    complex(8), intent(in)  :: h(ndim,ndim)
    integer,    intent(in)  :: nreig
    integer,    intent(in)  :: nsteps
    real(8),    intent(out) :: deval(nsteps)
    complex(8), intent(out) :: zevec(ndim,nsteps)

    ! local
    integer :: iter, lflag, info(4), lprnt, lpset(5)
    real(8) :: hnrm
    complex(8), allocatable :: basis(:,:)
    complex(8), allocatable :: q(:), r(:)
    real(8),    allocatable :: rnrm(:), work(:)

    external hlzdrd

    lprnt  = 64+128         ! LPRNT  is the printing code

    lpset(1) = ndim
    lpset(2) = nreig      
    lpset(3) = nsteps
    lpset(4) = lprnt
    lpset(5) = ndim

    allocate(rnrm(nsteps))
    allocate(work((nsteps+10)*nsteps))
    allocate(q(ndim), r(ndim))
    allocate(basis(ndim,nsteps))

    lflag = 0
    hnrm  = 0.0d0
  
    do while (.true.)

      call hlzdrd (lflag, lpset, info, hnrm, deval, rnrm, work, q, r, basis, zevec)

      if (lflag < 0) then

        write(*,'(a,4i5)') 'Abnormal exit, execution finished', info(1), info(2), info(3), info(4)
        if ( info(4) /= 0 ) stop

      else if (lflag == 1) then

        call zgemv ('n', ndim, ndim, zone, h, ndim, q, 1, zzero, r, 1)

      else

        exit ! iter loop

      end if

    end do

    ! if (iter > 1000) then
    !   write(*,*) 'too many iterations'
    !   stop
    ! end if

    call hlzdrd (3, lpset, info, hnrm, deval, rnrm, work, q, r, basis, zevec)
    write (*,'(a,4i5)') 'Standard exit, execution finished', info(1), info(2), info(3), info(4)

    deallocate(basis)
    deallocate(q, r, rnrm, work)

    return
  end subroutine

end module
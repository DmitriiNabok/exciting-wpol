
module mod_wpol_diagonalization

  implicit none

contains

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
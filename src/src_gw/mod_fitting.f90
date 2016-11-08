
!********************************************************************
! Find an approximating polynomial of known degree for a given data
!********************************************************************

module mod_fitting

contains
 
  function polyfit(vx, vy, d)
    implicit none
    integer, intent(in)                   :: d
    integer, parameter                    :: dp = 8
    real(dp), dimension(d+1)              :: polyfit
    real(dp), dimension(:), intent(in)    :: vx, vy
 
    real(dp), dimension(:,:), allocatable :: X
    real(dp), dimension(:,:), allocatable :: XT
    real(dp), dimension(:,:), allocatable :: XTX
 
    integer :: i, j
 
    integer :: n, lda, lwork
    integer :: info
    integer,  dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work
 
    n = d+1
    lda = n
    lwork = n
 
    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n, size(vx)))
    allocate(X(size(vx), n))
    allocate(XTX(n, n))
 
    ! prepare the matrix
    do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
    end do
 
    XT  = transpose(X)
    XTX = matmul(XT, X)
 
    ! calls to LAPACK subs DGETRF and DGETRI
    call DGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
    call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
 
    polyfit = matmul( matmul(XTX, XT), vy)
 
    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)
 
  end function

  subroutine test_polyfit()
    implicit none
     
    integer, parameter      :: degree = 2
    integer, parameter      :: dp = 8
    integer                 :: i
    real(dp), dimension(11) :: x = (/ (i,i=0,10) /)
    real(dp), dimension(11) :: y = (/ 1,   6,  17,  34, &
                                     57,  86, 121, 162, &
                                    209, 262, 321 /)
    real(dp), dimension(degree+1) :: a
 
    a = polyfit(x, y, degree)
 
    write (*, '(F9.4)') a
 
  end subroutine
 
end module
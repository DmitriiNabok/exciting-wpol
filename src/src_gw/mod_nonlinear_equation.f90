
module mod_nonlinear_equation

  implicit none

contains

  subroutine newton(func, x0, nmax, ftol, x, n, debug)
  !
  ! Search for a root of equation f(x)=0 using Newton's method 
  !
  ! Input:
  !       f: the function to find a root of that returns f(x) and f'(x)
  !      x0: the initial guess
  !    nmax: maximum number of iterations
  !    xtol: required tolerance
  !   debug: if .true. prints debugging info
  !
  ! Output:
  !       x: estimated solution
  !       n: actual number of iterations
  !
    implicit none
    external :: func
    real(8),  intent(in)  :: x0
    integer,  intent(in)  :: nmax
    real(8),  intent(in)  :: ftol
    real(8),  intent(out) :: x
    integer,  intent(out) :: n
    logical,  intent(in)  :: debug
    ! local variables
    integer :: i
    real(8) :: deltax, fx, dfx

    ! initial guess
    x = x0
    if (debug) write(*,'(" Initial guess: x = ", e22.15)') x

    ! Newton iteration to find a zero of f(x)
    do i = 1, nmax

      ! Evaluate function and its derivative
      call func(x, fx, dfx)

      ! Checke for convergence
      if (abs(fx) < ftol) exit

      ! Compute Newton increment x:
      deltax = fx / dfx

      ! update x
      x = x - deltax

      if (debug) write(*,'(" After ", i3, " iterations, x = ", e22.15)') i, x

    end do

    if (i == nmax) then
      ! might not have converged
      call func(x, fx, dfx)
      if (abs(fx) > ftol) write(*,*) 'WARNING(mod_nonlinear_equation::newton): Not converged!'
    end if
    
    n = i

    return
  end subroutine

  subroutine test_newton()
    implicit none
    real(8) :: x, x0, fx, dfx
    real(8) :: x0vals(3)
    integer :: n, i
    logical :: debug

    write(*,*)
    write(*,*) "**************************************"
    write(*,*) " Test routine for computing zero of f"
    write(*,*) "**************************************"

    debug = .true.

    ! values to test as x0
    x0vals = (/1.d0, 2.d0, 100.d0 /)

    do i = 1, 3

        write(*,*)
        write(*,*) "--------------------------------------------"
        
        x0 = x0vals(i)

        call newton(f_sqrt, x0, 20, 1.d-14, x, n, debug)

        write(*,*) "Solver returns:"
        write(*,"(' x = ', e22.15, ' after', i3, ' iterations')") x, n

        call f_sqrt(x, fx, dfx)
        write(*,"(' f(x) = ', e22.15)") fx
        write(*,*)

        if (abs(x-2.d0) > 1d-14) then
          write(*,*) "ERROR(mod_nonlinear_equation::test_newton):"
          write(*,*) "    Test failed!"
          write(*,*) "    Unexpected result: x = ", x
        end if

    end do

    return
  end subroutine

  subroutine f_sqrt(x, fx, dfx)
    real(8), intent(in)  :: x
    real(8), intent(out) :: fx
    real(8), intent(out) :: dfx
    fx  = x**2 - 4.d0
    dfx = 2.d0 * x
    return
  end subroutine

end module

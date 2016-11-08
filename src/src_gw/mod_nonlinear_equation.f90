
module mod_nonlinear_equation

  implicit none

contains

!***********************************************
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
!***********************************************
  subroutine newton(func, x0, nmax, ftol, x, n, debug)
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

  !----------------------------------------------------------
  subroutine test_newton()
    implicit none
    real(8) :: x, x0, fx, dfx
    real(8) :: x0vals(3)
    integer :: n, i
    logical :: debug

    write(*,*)
    write(*,*) "**************************************"
    write(*,*) " Newton's solver"
    write(*,*) "**************************************"

    debug = .true.

    ! values to test as x0
    x0vals = (/1.d0, 2.d0, 100.d0 /)

    do i = 1, 3

        write(*,*)
        write(*,*) "--------------------------------------------"
        
        x0 = x0vals(i)

        call newton(f_test1, x0, 20, 1.d-14, x, n, debug)

        write(*,*) "Solver returns:"
        write(*,"(' x = ', e22.15, ' after', i3, ' iterations')") x, n

        call f_test1(x, fx, dfx)
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


!***********************************************
!*         Mueller's method subroutine         *
!* ------------------------------------------- *
!* This routine iteratively seeks the root of  *
!* a function Y(x) by fitting a parabola to    *
!* three points and calculating the nearest    *
!* root as described in Becket and Hurt, nume- *
!* rical calculations and algorithms.          *
!* writeS:                                     *
!* The routine requires an initial guess, x0,  *
!* a bound on the error on this guess, d and a *
!* convergence criteria, e. Also required is a *
!* limit on the number of iterations, n.       *
!* OUTPUTS:                                    *
!* The routine returns the value of the root   *
!* found, x and the number of iterations per-  *
!* formed, k.                                  *
!***********************************************
  subroutine mueller(func, x0, d, e, n, x, k) 
    implicit none
    external :: func
    real(8), intent(in)  :: x0
    real(8), intent(in)  :: d
    real(8), intent(in)  :: e
    integer, intent(in)  :: n
    real(8), intent(out) :: x
    integer, intent(out) :: k
    ! local
    real(8), parameter :: eps = 1.d-16
    real(8) :: a1, b, c1, d1, e1, e2, e3
    real(8) :: x1, x2, x3, xl, xl1

    ! Set up the three evaluation points
    k = 1
    x3 = x0
    x1 = x3-d
    x2 = x3+d
    
    ! Calculate Mueller parameters and guard against divide by zero
    if (abs(x2-x1) < eps) x2 = x2*1.0000001
100 if (abs(x2-x1) < eps) x2 = x2+0.0000001
    xl1 = (x3-x2)/(x2-x1)
    d1  = (x3-x1)/(x2-x1)
    if (k > 1) goto 200

    ! Get values of function
    call func(x1, e1)
    call func(x2, e2)
200 call func(x3, e3)
    a1 = xl1*xl1*e1-d1*d1*e2+(xl1+d1)*e3
    c1 = xl1*(xl1*e1-d1*e2+e3)
    b  = a1*a1-4.d0*d1*c1*e3

    ! Test for complex root, meaning the parabola is inverted
    if (b<0) b = 0.d0
    ! Choose closest root
    if (a1 < 0) a1 = a1-dsqrt(b)
    if (a1 > 0) a1 = a1+dsqrt(b)
    ! Guard against divide by zero
    if (abs(a1) < eps) a1 = eps
    xl = -2.d0*d1*e3/a1
    ! Calculate next estimate
    x = x3+xl*(x3-x2)

    ! Test for convergence
    if (dabs(x-x3)<e) return
    ! Test for number of iterations
    if (k >= n) return

    ! Otherwise, make another pass
    k = k+1
    ! Some saving
    x1 = x2 
    x2 = x3 
    x3 = x 
    e1 = e2 
    e2 = e3
    goto 100
  end subroutine

!*****************************************************
!*     Program to demonstrate Mueller's method       *
!* ------------------------------------------------- *
!* Reference: BASIC Scientific Subroutines, Vol. II  *
!* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
!*                                                   *
!* Example:  Find a real root of f(x)=(x+1)^5        *
!*                                                   *
!* Sample run:                                       *
!*                                                   *
!*  Input the initial guess:                         *
!*     X0 = 0                                        *
!*  What is the bond of this guess: 3                *
!*  Error criterion: 1e-6                            *
!*  Maximum number of iterations: 100                *
!*                                                   *
!*  The calculated zero is X = -0.999164E+00         *
!*  The associated value is Y =  0.000000E+00        *
!*  The number of steps was: 44                      *
!*                                                   *
!*****************************************************
  subroutine test_mueller()
    implicit none
    integer :: k, n
    real(8) :: d, e, x, x0, y

    write(*,*)
    write(*,*) "**************************************"
    write(*,*) " Mueller's solver"
    write(*,*) "**************************************"

    x0 = 0.d0
    d  = 3.d0
    e  = 1.d-6
    n  = 100

    call mueller(f_test2, x0, d, e, n, x, k)

    write(*,50) x
    call f_test2(x, y)
    write(*,60) y
    write(*,70) k

50  format(' The calculated zero is X = ',e13.6)
60  format(' The associated Y value is Y = ',e13.6)
70  format(' The number of steps was: ',i4//)
  
  end subroutine


!***********************************************
!*          Aitken Method Subroutine           *
!* ------------------------------------------- *
!* This subroutine calculates the zeroes of a  *
!* function Y(x) by iterations, and employs    *
!* Aitken acceleration to speed up convergence.*
!* An initial guess is required, x0, and two   *
!* convergence factors, c and e. e relates to  *
!* the accuracy of the estimate, and c is used *
!* to aid convergence. Also required is the    *
!* maximum number of iterations, m. c=-1 is a  *
!* normal value, if divergence occurs, smaller *
!* and/or positive values should be tried.     *
!* The root is returned in x, the number of    *
!* iterations in n.                            *
!***********************************************
  subroutine aitken(func, m, n, c, e, x, x0)
    implicit none
    external :: func
    integer  :: m, n
    real(8)  :: c, e, x, x0
    real(8)  :: x1, x2, xk, yy
    real(8), parameter :: eps = 1.d-16
    n = 0
    x = x0
    ! Get y
100 call func(x, yy)
    yy = x+c*yy
    ! Enough points for acceleration? 
    if (n > 0) goto 200
    x1 = yy; x = x1; n = n+1
    goto 100
200 x2 = yy; n = n+1
    ! Guard against a zero denominator
    if (abs(x2-2.0*x1+x0) < eps) x0 = x0+1.d-8
    ! Perform acceleration
    xk = (x2-x1)*(x2-x1)/(x2-2.0*x1+x0)
    x2 = x2-xk
    ! Test for convergence
    if (n >= m) return
    if (abs(xk) < e) return
    x0 = x1; x1 = x2; x = x1
    goto 100
  end subroutine

!*****************************************************
!*           Program to demonstrate Aitken           *
!*                method subroutine                  *
!* ------------------------------------------------- *
!* Reference: BASIC Scientific Subroutines, Vol. II  *
!* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
!*                                                   *
!* Example: Find a real root of f(x)=(x+1)^5         *
!*                                                   *
!* Sample run:                                       *
!*                                                   *
!*  Input the initial guess:                         *
!*     X0 = 0                                        *
!*  Convergence criterion: 0.000001                  *
!*  Convergence factor: -1                           *
!*  Maximum number of iterations: 100                *
!*                                                   *
!*  The calculated zero is X = -1                    *
!*  The associated Y value is Y =  0                 *
!*  The number of iterations was:  2                 *
!*                                                   *
!*****************************************************
  subroutine test_aitken()
    implicit none
    integer :: m, n
    real(8) :: c, e, x, x0, y, dy

    write(*,*)
    write(*,*) "**************************************"
    write(*,*) " Aitken's solver"
    write(*,*) "**************************************"

    x0 = 0.d0
    c  = -1.d0
    e  = 1.d-6
    m  = 100

    call aitken(f_test2, m, n, c, e, x, x0)

    write(*,50) x
    call f_test2(x, y)
    write(*,60) y
    write(*,70) n

50  format(' The calculated zero is X = ',e13.6)
60  format(' The associated Y value is Y = ',e13.6)
70  format(' The number of steps was: ',i4//)
  
  end subroutine


!*****************************************************
!*              Brent Method Function                *
!* ------------------------------------------------- *
!* The purpose is to find a real root of a real      *
!* function f(x) using Brent method.                 *
!*                                                   *
!* INPUTS:  x1,x2     : interval of root             *
!*          Tolerance : desired accuracy for root    *
!*          maxIter   : maximum number of iterations *
!*                                                   *
!* OUTPUTS: The function returns the root value      *
!*          ValueAtRoot : value of f(root)           *
!*          niter    : number of done iterations     *
!*          error    : =0, all OK                    *
!*                   : =1, no root found in interval *
!*                   : =2, no more iterations !      *
!*****************************************************
  subroutine brent( func,                     & ! user supplied function
  &                 x1, x2,                   & ! search interval
  &                 tolerance, maxIterations, & ! convergence criteria
                    root, valueAtRoot,        & ! solution
                    niter, error )
    implicit none
    external :: func
    real(8), intent(in)  :: x1, x2
    real(8), intent(in)  :: tolerance
    integer, intent(in)  :: maxIterations
    real(8), intent(out) :: root, valueAtRoot
    integer, intent(out) :: niter, error

    real(8), parameter :: FPP = 1.d-8
    real(8), parameter :: nearzero = 1.d-16
    
    real(8) :: AA, BB, CC, DD, EE, FA, FB, FC, Tol1, PP, QQ, RR, SS, xm
    integer :: i, done

    i = 0; done = 0; error = 0
    AA = x1; call func(AA, FA)
    BB = x2; call func(BB, FB) 
    if (RootNotBracketed(FA,FB)) then 
      error = 1
    else 
      FC = FB
      do while (done == 0 .and. i < maxIterations)
        if (RootNotBracketed(FC,FB)) then
          CC = AA; FC = FA; DD = BB - AA; EE = DD
        end if
        if (dabs(FC) < dabs(FB)) then
          AA = BB; BB = CC; CC = AA
          FA = FB; FB = FC; FC = FA
        endif
        Tol1 = 2.0 * FPP * dabs(BB) + 0.5 * Tolerance
        xm = 0.5 * (CC-BB)
        if ((dabs(xm) <= Tol1).or.(dabs(FA) < nearzero)) then
          ! A root has been found
          root = BB
          done = 1
          call func(root, valueAtRoot)
        else 
          if ((dabs(EE) >= Tol1).and.(dabs(FA) > dabs(FB))) then
            SS = FB / FA;
            if (dabs(AA - CC) < nearzero) then
              PP = 2.0 * xm * SS
              QQ = 1.0 - SS
            else 
              QQ = FA / FC
              RR = FB / FC
              PP = SS * (2.0 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1.0))
              QQ = (QQ - 1.0) * (RR - 1.0) * (SS - 1.0)
            endif
            if (PP > nearzero) QQ = -QQ
            PP = dabs(PP)
            if ((2.0 * PP) < min(3.0*xm *QQ-dabs(Tol1 * QQ), dabs(EE * QQ))) then
              EE = DD;  DD = PP/QQ
            else 
              DD = xm;   EE = DD;
            endif
          else 
            DD = xm;
            EE = DD;
          endif
          AA = BB;
          FA = FB;
          if (dabs(DD) > Tol1) then 
            BB = BB + DD;
          else 
            if (xm > 0) then 
              BB = BB + dabs(Tol1)
            else 
              BB = BB - dabs(Tol1)
            endif
          endif
          call func(BB, FB)
          i = i+1
        endif
    end do
      if (i >= maxIterations) error = 2
    endif
    niter = i
  
  contains

    ! TRUE if x1*x2 negative
    logical function RootNotBracketed(x1, x2)
      real(8) :: x1, x2
      if ((x1 > 0 .and. x2 > 0).or.(x1 < 0 .and. x2 < 0)) then 
        RootNotBracketed = .true.
      else
        RootNotBracketed = .false.
      end if
    end function RootNotBracketed

  end subroutine

!*****************************************************
!*      Program to demonstrate the real domain       *
!*               brent subroutine                    *
!* ------------------------------------------------- *
!* Reference:  BORLAND MATHEMATICAL LIBRARY          *
!*                                                   *
!*                F90 version by J-P Moreau, Paris.  *
!*                       (www.jpmoreau.fr)           *
!* ------------------------------------------------- *
!* Example:    Find a real root of f(x)=(x+1)^5      *
!*                                                   *
!* SAMPLE RUN:                                       *
!*                                                   *
!*  Input interval (X1,X2):                          *
!*                                                   *
!*        X1 = -2                                    *
!*        X2 =  0                                    *
!*                                                   *
!*  Convergence criterion: 1e-10                     *
!*  Maximum number of iterations: 10                 *
!*                                                   *
!*  The estimated root is:                           *
!*                                                   *
!*        X = -1.000000                              *
!*                                                   *
!*  The associated Y value is Y = 0.000000           *
!*                                                   *
!*  The number of iterations was: 2                  *
!*  The error code is: 0                             *
!*                                                   *
!*****************************************************
  subroutine test_brent()
    implicit none
    integer :: m, n, error
    real(8) :: e, x, x1, x2, y

    write(*,*)
    write(*,*) "**************************************"
    write(*,*) " Brent's solver"
    write(*,*) "**************************************"

    x1 = -2.d0
    x2 =  0.d0
    e  =  1.d-6
    m  = 100

    call brent(f_test2, x1, x2, e, m, x, y, n, error)

    write(*,50) x
    write(*,60) y
    write(*,70) n
    write(*,80) error

50  format(' The calculated zero is X = ',e13.6)
60  format(' The associated Y value is Y = ',e13.6)
70  format(' The number of steps was: ',i4//)
80  format(' Exit status: ',i4//)

  end subroutine

!==========================================================
! Test functions
!==========================================================

  subroutine f_test1(x, fx, dfx)
    real(8), intent(in)            :: x
    real(8), intent(out)           :: fx
    real(8), intent(out), optional :: dfx
    fx  = x**2 - 4.d0
    if (present(dfx)) dfx = 2.d0 * x
    return
  end subroutine

  subroutine f_test2(x, fx, dfx)
    real(8), intent(in)            :: x
    real(8), intent(out)           :: fx
    real(8), intent(out), optional :: dfx
    fx  = (x+1.d0)**5
    if (present(dfx)) dfx = 5.d0*(x+1.d0)**4
    return
  end subroutine  

end module

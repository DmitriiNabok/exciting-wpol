
subroutine calc_evalqp_wpol()

  use modinput
  use modgw
  use mod_mpi_gw
  use m_getunit
  use mod_vxc,        only: vxcnn, read_vxcnn
  use mod_selfenergy, only: selfex, selfec, read_selfexnn, read_selfecnn, &
  &                         evalks, evalqp, sigc, znorm, eferqp, &
  &                         init_selfenergy, delete_selfenergy
  use mod_nonlinear_equation
  implicit none
  ! local
  integer :: ik, ib, n, error
  real(8) :: enk, eqp, znk, egap, df
  real(8) :: s, ds, w, e
  complex(8) :: sigma, dsigma

  input%gw%taskname = 'g0w0'
  input%gw%skipgnd = .true.
  call init_gw

  if (rank == 0) then

    if (allocated(vxcnn)) deallocate(vxcnn)
    allocate(vxcnn(ibgw:nbgw,kset%nkpt))
    call read_vxcnn()

    call init_selfenergy(ibgw,nbgw,kset%nkpt,freq%nomeg)
    call read_selfexnn()
    call read_selfecnn()

    !---------------------
    ! Solve QP equation
    !---------------------

    ! ik = 1
    ! ib = 1
    ! eks = evalsv(ib, ik)
    ! call sigma(freq%nomeg, freq%freqs, dble(selfec(ib,:,ik)), eks, s, ds)

    do ik = 1, kset%nkpt
      do ib = ibgw, nbgw
    
        enk = evalsv(ib,ik)

        if (.true.) then
          !---------------------
          ! iterative solution
          !---------------------
          e = 0.5 ! trust region
          ! call mueller(f1, enk, e, 1.d-6, 1000, eqp, n)
          call brent(f1, enk-e, enk+e, 1.d-6, 1000, eqp, s, n, error)
          ! call newton(f2, enk, 1000, 1.d-6, eqp, n, .false.)
          
          evalks(ib,ik) = enk
          evalqp(ib,ik) = eqp
          
          call getSelfc(freq%nomeg, freq%freqs, selfec(ib,:,ik), eqp, sigma, dsigma)
          sigc(ib,ik) = sigma
          znorm(ib,ik)  = 1.0d0/(1.0d0-dble(dsigma))

        else

          !---------------------
          ! linearized version
          !---------------------
          call getSelfc(freq%nomeg, freq%freqs, selfec(ib,:,ik), enk, sigma, dsigma)
          sigc(ib,ik) = sigma

          ! Set the renormalization factor              
          znk = 1.0d0/(1.0d0-dble(dsigma))
          if ((znk > 1.d0) .or. (znk < 0.5d0)) then
            write(fgw,*) 'WARNING(calcevalqp):'
            write(fgw,100) ik, ib, enk, znk, sigma, dsigma
            100 format(' Suspicious Znk',' irk=',i4,' ie=',i4, &
            &          ' enk=',f8.3,' eV',' Znk=',f6.2," ReSc=",f8.3,  &
            &          ' ImSc=',f8.3," ReSc'=",f8.3," ImSc'=",f8.3)
            write(fgw,*)
            ! znk = 1.d0  ! set a default value
          endif
          znorm(ib,ik)  = znk

          evalks(ib,ik) = enk
          evalqp(ib,ik) = enk + &
          &               znk * dble( selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik) )

        end if

      end do
      
    end do    

    ! Calculate Fermi energy
    call fermi_exciting(input%groundstate%tevecsv, nvelgw, &
    &                   nbandsgw, kset%nkpt, evalqp(ibgw:nbgw,:), &
    &                   kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
    &                   eferqp, egap, df)
    
    call write_qp_energies('EVALQP.DAT')
    call bandstructure_analysis('G0W0', ibgw, nbgw, kset%nkpt, evalqp(ibgw:nbgw,:), eferqp)

    deallocate(vxcnn)
    call delete_selfenergy()

  end if ! rank == 0


contains

  subroutine f1(x, fx)
    real(8), intent(in)  :: x
    real(8), intent(out) :: fx
    ! local
    complex(8) :: s, ds
    call getSelfc(freq%nomeg, freq%freqs, selfec(ib,:,ik), x, s, ds)
    fx  = evalsv(ib,ik) + dble(selfex(ib,ik) - vxcnn(ib,ik) + s) - x
    return
  end subroutine

  subroutine f2(x, fx, dfx)
    real(8), intent(in)  :: x
    real(8), intent(out) :: fx
    real(8), intent(out) :: dfx
    ! local
    complex(8) :: s, ds
    call getSelfc(freq%nomeg, freq%freqs, selfec(ib,:,ik), x, s, ds)
    fx  = evalsv(ib,ik) + dble(selfex(ib,ik) - vxcnn(ib,ik) + s) - x
    dfx = dble(ds) - 1.d0
    return
  end subroutine

  subroutine getSelfc(n, x, y, x0, f, df)
    integer,    intent(in)  :: n
    real(8),    intent(in)  :: x(n)
    complex(8), intent(in)  :: y(n)
    real(8),    intent(in)  :: x0
    complex(8), intent(out) :: f
    complex(8), intent(out) :: df
    ! local
    integer    :: i, i0, ii, j
    integer    :: np, np2
    real(8)    :: dx, t1, t2, xx
    complex(8) :: dy
    real(8),    allocatable :: c(:)
    complex(8), allocatable :: ya(:)
    real(8),    external    :: polynom

    ! only 1 point
    if (n <= 1) then
      write(*,*) 'ERROR(calc_evalqp_wpol::sigma): Not enough data points!'
      stop
    end if

    xx = max(x0, x(1))
    xx = min(x0, x(n))

if (.true.) then
    !-----------------------
    ! linear interpolation
    !-----------------------
    do i = 1, n
      if (x(i) >= xx) then
        if (i == 1) then
          df = (y(2)-y(1)) / (x(2)-x(1))
          f  = y(1) + df * (xx-x(1))
        else if (i == n) then
          df = (y(n)-y(n-1)) / (x(n)-x(n-1))
          f  = y(n-1) + df * (xx-x(n-1))
        else
          df = 0.5*(y(i+1)-y(i-1)) / (x(i)-x(i-1))
          f  = y(i-1) + df * (xx-x(i-1))
        end if
        exit
      end if
    end do
else
    !-----------------------
    ! polynomial fitting
    !-----------------------
    np = 3
    np2 = np/2
    allocate(ya(np), c(np))
    do i = 1, n
      if (x(i) >= xx) then
        if (i <= np2) then
          i0 = 1
        else if (i > n-np2) then
          i0 = n-np+1
        else
          i0 = i-np2
        end if
        do j = 1, np
          ii = i0+j-1
          ya(j) = y(ii)
        end do
        t1 = polynom( 0, np, x(i0), dble(ya), c, xx)
        t2 = polynom( 0, np, x(i0), aimag(ya), c, xx)
        f  = cmplx(t1, t2, 8)
        t1 = polynom( 1, np, x(i0), dble(ya), c, xx)
        t2 = polynom( 1, np, x(i0), aimag(ya), c, xx)
        df = cmplx(t1, t2, 8)
        exit
      end if
    end do
    deallocate(ya, c)
end if

    return
  end subroutine


end
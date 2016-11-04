
subroutine calc_evalqp_wpol()

  use modinput
  use modgw
  use mod_mpi_gw
  use m_getunit
  use mod_vxc,                only: vxcnn, read_vxcnn
  use mod_selfenergy,         only: selfex, selfec, read_selfexnn, read_selfecnn
  use mod_nonlinear_equation, only: newton
  implicit none
  ! local
  integer :: ik, ib, n, i
  real(8) :: eks, eqp
  real(8) :: s, ds, w

  input%gw%skipgnd = .true.
  call init_gw

  if (rank == 0) then

    if (allocated(vxcnn)) deallocate(vxcnn)
    allocate(vxcnn(ibgw:nbgw,kset%nkpt))
    vxcnn(:,:) = 0.d0
    call read_vxcnn(kset,ibgw,nbgw,vxcnn)

    if (allocated(selfex)) deallocate(selfex)
    allocate(selfex(ibgw:nbgw,nkpt))
    selfex(:,:) = 0.d0
    call read_selfexnn(kset,ibgw,nbgw,selfex)
  
    if (allocated(selfec)) deallocate(selfec)
    allocate(selfec(ibgw:nbgw,freq%nomeg,kset%nkpt))
    selfec(:,:,:) = 0.d0
    call read_selfecnn(kset,freq,ibgw,nbgw,selfec)

    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(ibgw:nbgw,1:kset%nkpt))
    evalqp(ibgw:nbgw,1:kset%nkpt) = evalsv(ibgw:nbgw,1:kset%nkpt) + &
    &                               selfex(ibgw:nbgw,1:kset%nkpt) - &
    &                                vxcnn(ibgw:nbgw,1:kset%nkpt)

    !---------------------
    ! Solve QP equation
    !---------------------

    ! ik = 1
    ! ib = 1
    ! eks = evalsv(ib, ik)
    ! call sigma(freq%nomeg, freq%freqs, dble(selfec(ib,:,ik)), eks, s, ds)

    ! iterative solution
    do ik = 1, kset%nkpt
      do ib = ibgw, nbgw
    
        eks = evalsv(ib,ik)
    
        call newton(func, eks, 100, 1.d-6, eqp, n, .true.)
        write(*,*) i, eks, eqp

        stop

      end do
    end do
    


  end if ! rank == 0

contains

  subroutine func(x, fx, dfx)
    real(8), intent(in)  :: x
    real(8), intent(out) :: fx
    real(8), intent(out) :: dfx
    call sigma(freq%nomeg, freq%freqs, dble(selfec(ib,:,ik)), x, fx, dfx)
    fx  = fx - x - evalqp(ib,ik)
    dfx = dfx - 1.d0
    return
  end subroutine

  subroutine sigma(n, w, d, e, s, ds)
    integer, intent(in)  :: n
    real(8), intent(in)  :: w(n)
    real(8), intent(in)  :: d(n)
    real(8), intent(in)  :: e
    real(8), intent(out) :: s
    real(8), intent(out) :: ds
    ! local
    real(8), parameter   :: eps = 1.d-8
    integer :: iom, iom0, np, np2, i, j
    real(8), allocatable :: x(:), y(:), c(:)
    real(8), external :: polynom

    ! only 1 point
    if (n == 1) then
      write(*,*) 'ERROR(calc_evalqp_wpol::sigma): Not enough frequency points!'
      stop
    end if

    ! check the input w
    if ((e < w(1)).or.(e > w(n))) then
      write(*,*) 'ERROR(calc_evalqp_wpol::sigma): w is out of the frequency range!'
      stop
    end if    

    ! order of predictor-corrector polynomial
    np = 4
    np2 = np/2
    allocate(x(np), y(np), c(np))

    do iom = 1, n
      if (w(iom) >= e) then
        if (iom <= np2) then
          iom0 = 1
        else if (iom > n-np2) then
          iom0 = n - np + 1
        else
          iom0 = iom - np2
        end if
        do j = 1, np
          i = iom0+j-1
          x(j) = w(i)
          y(j) = d(i)
        end do
        s  = polynom(0, np, x, y, c, e)
        ds = polynom(1, np, x, y, c, e)

        ! write(*,*) iom
        ! write(*,*) x
        ! write(*,*) y
        ! write(*,*) s
        ! write(*,*) ds

        exit
      end if
    end do

    deallocate(x, y, c)
    return
  end subroutine

end

subroutine test()
  use mod_lanczos
  use mod_rank1_update
  use mod_nonlinear_equation
  implicit none

  ! call test_lanczos()
  ! call test_rank1_update()
  
  ! call test_newton()
  call test_mueller()
  call test_aitken()
  call test_brent()

  return
end
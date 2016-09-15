
subroutine gwtasklauncher()
   
  use modinput
  use modmain, only: task
  use inputdom

  if (associated(input%gw)) then
      call rereadinput
      task = 1
      if (associated(input%gw%gwplan)) then
        do i = 1, size(input%gw%gwplan%taskarray)
          input%gw%taskname = input%gw%gwplan%taskarray(i)%task%name
          call gw_main
        end do
      else
        call gw_main
      end if
  else
     write (*,*) "error gwtasklauncher"
     stop
  end if

  return
end subroutine gwtasklauncher

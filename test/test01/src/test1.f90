program test
use modreport
implicit none
testplan_name ="test1"
tdirectory="test01"

call inittestoutputfile(50) !file unit 50

! list test routines here and call testreport(testunit,input,output,passed)
! before leaving the routine
call test_comparefiles()
call test_array_assert_equal()
call test_readinput()
call test_gndstate_init()
call test_GEOMETRY()
call test_EFERMI()
call test_EQATOMS()
call test_LINENGY()

call finalizeoutput(50)

end program
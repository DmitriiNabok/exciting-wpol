subroutine test_comparefiles()
use modreport
implicit none
character(256)::file1,file2,diff
logical::passed
file1="../reference/GEOMETRY.OUT"
file2="../reference/GEOMETRY.OUT"
diff="../output/diff"
call comparefiles(file1,file2,diff,passed)
testunitname="comparefiles"
inputf="GEOMETRY.OUT"
outputf=trim(diff)
call testreport(passed)
end subroutine
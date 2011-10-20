
module modsymmetries
use inputdom
implicit none
type symmetries_type
 character(512)::HermannMauguinSymbol
 character(512)::HallSymbol
 character(512)::SchoenfliesSymbol
 character(512)::spaceGroupNumber
 character(512)::title
  type(lattice_type),pointer::lattice
  type(WyckoffPositions_type),pointer::WyckoffPositions
end type
type lattice_type
 real(8)::a
 real(8)::b
 real(8)::c
 real(8)::ab
 real(8)::ac
 real(8)::bc
 integer::ncell(3)
 real(8)::scale
 real(8)::stretch(3)
 real(8)::epslat
 logical::primcell
 character(512)::speciespath
end type
type WyckoffPositions_type
  type(wspecies_type_array),pointer::wspeciesarray(:)
end type
type wspecies_type
 character(512)::speciesfile
  type(wpos_type_array),pointer::wposarray(:)
end type

type wspecies_type_array
type(wspecies_type),pointer::wspecies
 end type
    type wpos_type
 real(8)::coord(3)
end type

type wpos_type_array
type(wpos_type),pointer::wpos
 end type
    
   type(symmetries_type)::symmetries
contains

function getstructsymmetries(thisnode)

implicit none
type(Node),pointer::thisnode
type(symmetries_type),pointer::getstructsymmetries
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructsymmetries)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at symmetries"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"HermannMauguinSymbol")
if(associated(np)) then
       call extractDataAttribute(thisnode,"HermannMauguinSymbol",getstructsymmetries%HermannMauguinSymbol)
       call removeAttribute(thisnode,"HermannMauguinSymbol")  
        else
        write(*,*)"Parser ERROR: The element 'symmetries' requires the attribute 'HermannMauguinSymbol' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"HallSymbol")
if(associated(np)) then
       call extractDataAttribute(thisnode,"HallSymbol",getstructsymmetries%HallSymbol)
       call removeAttribute(thisnode,"HallSymbol")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"SchoenfliesSymbol")
if(associated(np)) then
       call extractDataAttribute(thisnode,"SchoenfliesSymbol",getstructsymmetries%SchoenfliesSymbol)
       call removeAttribute(thisnode,"SchoenfliesSymbol")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"spaceGroupNumber")
if(associated(np)) then
       call extractDataAttribute(thisnode,"spaceGroupNumber",getstructsymmetries%spaceGroupNumber)
       call removeAttribute(thisnode,"spaceGroupNumber")  
endif

            len= countChildEmentsWithName(thisnode,"lattice")

        if(len.eq.0) then
        write(*,*)"Parser ERROR: The symmetries element must contain at least 1 lattice element"
        endif
        getstructsymmetries%lattice=>null()
Do i=0,len-1
getstructsymmetries%lattice=>getstructlattice(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"lattice"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"WyckoffPositions")
getstructsymmetries%WyckoffPositions=>null()
Do i=0,len-1
getstructsymmetries%WyckoffPositions=>getstructWyckoffPositions(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"WyckoffPositions"),0)) ) 
enddo

      len= countChildEmentsWithName (thisnode,"title")
if (len .lt. 1) then
  write(*,*) "Parser ERROR: "
  Write (*,*)"The Element: title must occur at least 1 times in the"
   Write (*,*) "symmetries element"
  stop
endif
Do i=1,len

getstructsymmetries%title=getvalueoftitle(&
      removechild(thisnode,item(getElementsByTagname(thisnode,&
      "title"),0)))
end do

      i=0
      len=0
      
      call  handleunknownnodes(thisnode)
end function

function getstructlattice(thisnode)

implicit none
type(Node),pointer::thisnode
type(lattice_type),pointer::getstructlattice
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructlattice)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at lattice"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"a")
if(associated(np)) then
       call extractDataAttribute(thisnode,"a",getstructlattice%a)
       call removeAttribute(thisnode,"a")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'a' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"b")
if(associated(np)) then
       call extractDataAttribute(thisnode,"b",getstructlattice%b)
       call removeAttribute(thisnode,"b")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'b' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"c")
if(associated(np)) then
       call extractDataAttribute(thisnode,"c",getstructlattice%c)
       call removeAttribute(thisnode,"c")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'c' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ab")
if(associated(np)) then
       call extractDataAttribute(thisnode,"ab",getstructlattice%ab)
       call removeAttribute(thisnode,"ab")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'ab' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ac")
if(associated(np)) then
       call extractDataAttribute(thisnode,"ac",getstructlattice%ac)
       call removeAttribute(thisnode,"ac")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'ac' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"bc")
if(associated(np)) then
       call extractDataAttribute(thisnode,"bc",getstructlattice%bc)
       call removeAttribute(thisnode,"bc")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'bc' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ncell")
getstructlattice%ncell=(/1,1,1/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"ncell",getstructlattice%ncell)
       call removeAttribute(thisnode,"ncell")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"scale")
getstructlattice%scale=1
if(associated(np)) then
       call extractDataAttribute(thisnode,"scale",getstructlattice%scale)
       call removeAttribute(thisnode,"scale")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"stretch")
getstructlattice%stretch=(/1.0d0,1.0d0,1.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"stretch",getstructlattice%stretch)
       call removeAttribute(thisnode,"stretch")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epslat")
getstructlattice%epslat=1.0d-6
if(associated(np)) then
       call extractDataAttribute(thisnode,"epslat",getstructlattice%epslat)
       call removeAttribute(thisnode,"epslat")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"primcell")
getstructlattice%primcell= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"primcell",getstructlattice%primcell)
       call removeAttribute(thisnode,"primcell")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"speciespath")
getstructlattice%speciespath= "http://xml.exciting-code.org/species/"
if(associated(np)) then
       call extractDataAttribute(thisnode,"speciespath",getstructlattice%speciespath)
       call removeAttribute(thisnode,"speciespath")  
endif

      i=0
      len=0
      
      call  handleunknownnodes(thisnode)
end function

function getstructWyckoffPositions(thisnode)

implicit none
type(Node),pointer::thisnode
type(WyckoffPositions_type),pointer::getstructWyckoffPositions

integer::len=1,i=0
allocate(getstructWyckoffPositions)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at WyckoffPositions"
#endif
      
            len= countChildEmentsWithName(thisnode,"wspecies")
     
allocate(getstructWyckoffPositions%wspeciesarray(len))
Do i=0,len-1
getstructWyckoffPositions%wspeciesarray(i+1)%wspecies=>getstructwspecies(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"wspecies"),0)) ) 
enddo

      i=0
      len=0
      
      call  handleunknownnodes(thisnode)
end function

function getstructwspecies(thisnode)

implicit none
type(Node),pointer::thisnode
type(wspecies_type),pointer::getstructwspecies
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructwspecies)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at wspecies"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"speciesfile")
if(associated(np)) then
       call extractDataAttribute(thisnode,"speciesfile",getstructwspecies%speciesfile)
       call removeAttribute(thisnode,"speciesfile")  
endif

            len= countChildEmentsWithName(thisnode,"wpos")
     
allocate(getstructwspecies%wposarray(len))
Do i=0,len-1
getstructwspecies%wposarray(i+1)%wpos=>getstructwpos(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"wpos"),0)) ) 
enddo

      i=0
      len=0
      
      call  handleunknownnodes(thisnode)
end function

function getstructwpos(thisnode)

implicit none
type(Node),pointer::thisnode
type(wpos_type),pointer::getstructwpos
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructwpos)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at wpos"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"coord")
if(associated(np)) then
       call extractDataAttribute(thisnode,"coord",getstructwpos%coord)
       call removeAttribute(thisnode,"coord")  
endif

      i=0
      len=0
      
      call  handleunknownnodes(thisnode)
end function
 
function getvalueoftitle(thisnode)
implicit none
type(Node),pointer::thisnode
 character(512)::getvalueoftitle

#ifdef INPUTDEBUG
  write(*,*)"we are at title"
#endif  
   call extractDataContent(thisnode,  getvalueoftitle)
end function

function countChildEmentsWithName(nodep,name)
  implicit none
  integer::countChildEmentsWithName
  type(Node),pointer ::nodep
  character(len=*),intent(in)::name
  type(NodeList),pointer::children
  type(Node),pointer::child
  
  integer::i
  children=>getChildNodes(nodep)
  countChildEmentsWithName=0
  do i=0,getlength(children)-1
    child=>item(children,i)
    if(name.eq.getNodeName(child)) countChildEmentsWithName=countChildEmentsWithName+1
  end do

end function  
    
end module


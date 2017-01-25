subroutine bandstructure_info(title,ib,nb,nkpt,eband,efermi)

    use modgw,     only : kset, fgw, hev
    use mod_bands, only : metallic
    use modmpi,    only : rank

    implicit none
    character(len=*), intent(in) :: title
    integer, intent(in)    :: ib, nb, nkpt
    real(8), intent(in)    :: eband(ib:nb,nkpt)
    real(8), intent(inout) :: efermi
    ! local variables
    real(8) :: ebmax, ebmin, egf, ego
    real(8) :: fermidos

    integer :: i, j, ik, jk
    integer :: bindx(2), kindx(2)
    real(8) :: egap, t1

    real(8), external :: dostet_exciting
   
    if (rank==0) then
      call boxmsg(fgw,'-',trim(title))
      write(fgw,'(a,f10.4)') " Fermi energy: ", efermi
    end if
    
    !-------------------------------------------------------------------
    ! check Fermi energy for correspondence to the specified band range 
    !-------------------------------------------------------------------
    ebmin = minval(eband)
    ebmax = maxval(eband)
    if ((ebmax < efermi) .or. (ebmin > efermi)) then 
        write(*,*) "ERROR(bandstructure_analysis): Fermi energy is outside the specified electronic bands energy range!"
        stop
    end if
    if (rank==0) then
      write(fgw,'(a,2f10.4)') " Energy range: ", ebmin, ebmax
    end if

    ! Determine the bandgap
    egap = 1000.d0
    do ik = 1, nkpt
    do i = ib, nb
    if ( eband(i,ik) <= efermi ) then

      do jk = 1, nkpt
      do j = nb, ib, -1
      if ( eband(j,jk) > efermi ) then

        t1 = eband(j,jk)-eband(i,ik)

        if (t1 < 0.d0) then
          metallic = .true.
          goto 100
        end if

        if ( ( t1 < egap ) .or. ( abs(t1-egap)<1.d-8 ) ) then
          egap = t1
          bindx(1) = i
          bindx(2) = j
          kindx(1) = ik
          kindx(2) = jk
        end if

      end if
      end do
      end do 

    end if
    end do
    end do

100 continue
    
    ! Calculate DOS at the fermi level
    fermidos = dostet_exciting(nb-ib+1,nkpt,eband, &
    &                          kset%ntet,kset%tnodes,kset%wtet,kset%tvol, &
    &                          efermi)
    
    ! check for VBM and CBM overlap (metal)
    if ( abs(fermidos)>1.d-4 ) then
        metallic = .true.
    else
        metallic = .false.
    end if
    
    if (rank==0) then
      
      ! write(fgw,*) egap
      ! write(fgw,*) bindx
      ! write(fgw,*) kindx

      if (metallic) then
        write(fgw,'(a,f8.4)') " DOS at Fermi level: ", fermidos
        write(fgw,*) "WARNING(bandstructure_analysis): Valence and Conduction bands overlap (metal)!"
      else
        if ( kindx(1) == kindx(2) ) then 
          ! direct gap
          write(fgw,10) ' Direct BandGap (eV):', egap*hev
          write(fgw,11) kset%vkl(:,kindx(1)), kindx(1)
        else
          ! indirect gap
          write(fgw,10) ' Indirect BandGap (eV):', egap*hev
          write(fgw,12) kset%vkl(:,kindx(1)), kindx(1), kset%vkl(:,kindx(2)), kindx(2)
          write(fgw,10) ' Direct Bandgap at VBM (eV):', (eband(bindx(1)+1,kindx(1))-eband(bindx(1),kindx(1)))*hev
          write(fgw,10) ' Direct Bandgap at CBM (eV):', (eband(bindx(2),kindx(2))-eband(bindx(2)-1,kindx(2)))*hev
        end if
      end if
      call linmsg(fgw,'-','')
      call flushifc(fgw)
    end if
    
10  format(a,T40,f10.4)
11  format(' at k      = ',3f8.3,' ik = ',i5)
12  format(' at k(VBM) = ',3f8.3,' ik = ',i5,/,&
&          '    k(CBM) = ',3f8.3,' ik = ',i5)
    
    return
end subroutine

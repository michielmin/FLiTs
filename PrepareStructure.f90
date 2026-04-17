subroutine PrepareStructure()
  use GlobalSetup
  use Constants
  use InOut
  implicit none
  integer :: i, j, imol, ispec(nmol), ipop

  call output("==================================================================")
  call output("PrepareStructure ...")

  do i = 0, nR
    C(i, 0)%dens = 1d-50
  end do
  do j = 0, nTheta
    C(0, j)%dens = 1d-50
    if (cylindrical) C(nR, j)%dens = 1d-50
  end do

  do imol = 1, nmol
    Mol(imol)%LTE = LTE
  end do

  ! now the data should be rearranged properly
  do imol = 1, nmol
    do i = nspec, 1, -1
      if (trim(Mol(imol)%name) == trim(mol_name0(i))) exit
    end do
    ispec(imol) = i
    if (ispec(imol) < 1 .or. ispec(imol) > nspec) then
      call output("Species "//trim(Mol(imol)%name)//" not found")
      !   if (.not. LTE) call output("Switching to LTE for this species")
      call output("removing this species")
      Mol(imol)%LTE = .true.
    end if
  end do

  call output("Rearranging data and determine line_width ...") 
  do i = 0, nR
    call tellertje(i + 1, nR + 1)
    do j = 0, nTheta
      !   allocate(C(i,j)%npop(nmol,maxlevels))
      do imol = 1, nmol
        ! this is the case where we do have a lamda file, but species is not in the fits file
        if (ispec(imol) < 1 .or. ispec(imol) > nspec) then
          C(i, j)%N(imol) = 1d-70
          C(i, j)%line_width(imol) = 1d5
          allocate (C(i, j)%npop(imol)%N(Mol(imol)%nlevels), source=0.d0)
        else
          if (C(i, j)%line_width(imol) < vres_profile*3d0) C(i, j)%line_width(imol) = 3d0*vres_profile
        end if
      end do
    end do
  end do

  call output("Rearranging data and computing maximum population levels ...")
  do i = 0, nR
    call tellertje(i + 1, nR + 1)
    do j = 0, nTheta
      C(i, j)%iT = C(i, j)%Tdust + 0.5d0
      if (C(i, j)%iT < 1) C(i, j)%iT = 1
      if (C(i, j)%iT > MAXT) C(i, j)%iT = MAXT
      allocate (C(i, j)%npopmax(nmol))
      do imol = 1, nmol
        C(i, j)%npopmax(imol) = 1
        if (ispec(imol) >= 1 .and. ispec(imol) <= nspec) then
          do ipop = size(C(i, j)%npop(imol)%N), 1, -1
            if (C(i, j)%npop(imol)%N(ipop) > 1d-150) exit
          end do
          if (ipop > C(i, j)%npopmax(imol)) C(i, j)%npopmax(imol) = ipop
        end if
      end do
    end do
  end do

  ! compute LTE where needed or requested
  call ComputeLTE()

  ! everything is read in now
  ! output the setup to the screen and the log file
  call output("==================================================================")

  rlines = 1d0/(sqrt((1d0 + vresolution/clight)/(1d0 - vresolution/clight)) - 1d0)

  call output("Mass of the star:   "//trim(dbl2string(Mstar, '(f13.4)'))//" Msun")
  call output("Minimum wavelength: "//trim(dbl2string(lmin, '(f13.4)'))//" micron")
  call output("Maximum wavelength: "//trim(dbl2string(lmax, '(f13.4)'))//" micron")
  call output("Resolution lines:   "//trim(dbl2string(rlines, '(f13.4)'))//" (dlam/lam)")
  call output("Velocity resolution:"//trim(dbl2string(vresolution/1d5, '(f13.4)'))//" km/s")
  call output("Inclination angle:  "//trim(dbl2string(inc, '(f13.4)'))//" degrees")

  call output("==================================================================")
  do i = 1, nmol
    call output("Line file: "//trim(linefile(i)))
  end do

  return
end subroutine PrepareStructure

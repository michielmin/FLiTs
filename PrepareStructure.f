	subroutine PrepareStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,k,imol,ispec,maxlevels
	
	do i=0,nR
		C(i,0)%dens=1d-50
	enddo
	do j=0,nTheta
		C(0,j)%dens=1d-50
		if(cylindrical) C(nR,j)%dens=1d-50
	enddo

	do imol=1,nmol
		Mol(imol)%LTE=LTE
	enddo

	maxlevels=0
c now the data should be rearranged properly
	do imol=1,nmol
		do ispec=1,nspec
			if(trim(Mol(imol)%name).eq.trim(mol_name0(ispec))) exit
		enddo
		if(ispec.gt.nspec) then
			call output("Species " // trim(Mol(imol)%name) // " not found")
c			if(.not.LTE) call output("Switching to LTE for this species")
			call output("removing this species")
			Mol(imol)%LTE=.true.
		else
			if(npop0(ispec).gt.maxlevels) maxlevels=npop0(ispec)
		endif
		if(Mol(imol)%nlevels.gt.maxlevels) maxlevels=Mol(imol)%nlevels
	enddo

	do i=0,nR
		do j=0,nTheta
			allocate(C(i,j)%N(nmol))
			allocate(C(i,j)%npop(nmol,maxlevels))
			allocate(C(i,j)%line_width(nmol))
		enddo
	enddo

	do imol=1,nmol
		do ispec=1,nspec
			if(trim(Mol(imol)%name).eq.trim(mol_name0(ispec))) exit
		enddo
		if(ispec.gt.nspec) then
			do i=0,nR
				do j=0,nTheta
					C(i,j)%N(imol)=1d-70
					C(i,j)%line_width(imol)=1d5
					do k=1,Mol(imol)%nlevels
						C(i,j)%npop(imol,k)=0d0
					enddo
				enddo
			enddo
		else
			do i=0,nR
				do j=0,nTheta
					C(i,j)%N(imol)=C(i,j)%N0(ispec)
					C(i,j)%line_width(imol)=C(i,j)%line_width0(ispec)
					do k=1,npop0(ispec)
						C(i,j)%npop(imol,k)=C(i,j)%npop0(ispec)%N(k)
					enddo
				enddo
			enddo
		endif
	enddo

	do i=0,nR
		do j=0,nTheta
			C(i,j)%iT=C(i,j)%Tdust+0.5d0
			if(C(i,j)%iT.lt.1) C(i,j)%iT=1
			if(C(i,j)%iT.gt.MAXT) C(i,j)%iT=MAXT
		enddo
	enddo

c compute LTE where needed or requested
	call ComputeLTE()
	
c everything is read in now
c output the setup to the screen and the log file
	call output("==================================================================")

	rlines=1d0/(sqrt((1d0+vresolution/clight)/(1d0-vresolution/clight))-1d0)

	call output("Mass of the star:   "//trim(dbl2string(Mstar,'(f13.4)'))//" Msun")
	call output("Minimum wavelength: "//trim(dbl2string(lmin,'(f13.4)'))//" micron")
	call output("Maximum wavelength: "//trim(dbl2string(lmax,'(f13.4)'))//" micron")
	call output("Resolution lines:   "//trim(dbl2string(rlines,'(f13.4)'))//" (dlam/lam)")
	call output("Velocity resolution:"//trim(dbl2string(vresolution/1d5,'(f13.4)'))//" km/s")
	call output("Inclination angle:  "//trim(dbl2string(inc,'(f13.4)'))//" degrees")

	call output("==================================================================")
	do i=1,nmol
		call output("Line file: "//trim(linefile(i)))
	enddo
	
	return
	end


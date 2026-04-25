  subroutine RaytraceLines()
    use GlobalSetup
    use Constants
    use InOut
    implicit none
    integer :: i, j,igrid, ilam, k, iblends, vmult, iv, nv, nl, imol, maxblend, ilines, nvmax, nvmin
    integer :: nb, ib0, nb0, ib, nltot
    integer(kind=4) :: counts, count_rate, count_max
    integer, external :: OMP_GET_THREAD_NUM
    integer, allocatable :: imol_blend(:), count(:)
    integer :: ilam_cube_gap, ilam_cube_start, ilam_cube_end
    real(kind=8), allocatable :: v_blend(:), flux4(:)
    real(kind=8) :: lam, flux0, starttime, stoptime, fact, lcmin
    real(kind=8) :: lmin_next, lam_velo
    real(kind=8), allocatable :: flux(:), flux_cont(:)
    type(Path), pointer :: PP
    type(Line) :: LL
    type(Blend), pointer :: Bl
    logical :: gas, doit
    logical, allocatable :: doit_ib(:), doit_ib0(:)
    real(kind=8) :: flux1, flux2, flux3, fc, f, dnu, lam_w, lam_w_min, lam_w_max
    real(kind=8) :: vlo_eff, vhi_eff
    real(kind=8) :: wl11, wl21, wl13, wl23, flux_l1, flux_l2, flux_c
    character(len=1000) :: comment
    character(len=200) :: callerstr
    real(kind=8) :: delpix,starcorr,lamgapmax    
    real(kind=8) :: time, utime, timegap, utimegap    

    call output("==================================================================")
    call output("Preparing the profiles")

    do i = 1, ngrids
      do j = 1, npoints(i)
        allocate (P(i, j)%cont_contr(P(i, j)%n))
        allocate (P(i, j)%exptau_dust(P(i, j)%n))
        allocate (P(i, j)%S_dust(P(i, j)%n))
      end do
    end do

    open (unit=20, file=outputFile, RECL=1500)
    open (unit=21, file=output_lineFluxFile, RECL=1500)
    write (21, '("# lam[mu]   Fline[W/m^2]    lmin[mu]    lmax[mu]    comment")')

    ! number of velocity points is determined by the width of the broadest possible line
    ! which is given by vmax (which should be the max projected Keplerian velocity)
    nv = int(vmax*1.1/vresolution) + 1
    nvprofile = int(vmax*vres_mult/vresolution)

    ! this should always give nvmin=-nv and nmvmax=nv
    call DetermineBlends(nv, maxblend, nvmin, nvmax)

    allocate (flux(nvmin:nvmax))
    allocate (flux_cont(nvmin:nvmax))
    allocate (imol_blend(maxblend))
    allocate (v_blend(maxblend))
    allocate (count(nmol))

    if (imagecube) then
      call clock(time, utime)
      call output("")
      call output("Preparing the image cube ...")
      call output("Allocating the image cube ...")
      dlam_cube = 0.5d0*(lmin + lmax)*vresolution/clight
      write(*,*) "dlam_cube: ", dlam_cube
      nlam_cube = int((lmax - lmin)/dlam_cube) ! number of center wl-points
      allocate (imcube(npix, npix, nlam_cube))      
      allocate (lam_cube(nlam_cube))
      allocate (lamb_cube(nlam_cube+1))
      ! for convenience store the center wl points
      ! lmin should be the first center
      do i = 1, nlam_cube
        lam_cube(i) = lmin + dlam_cube*(i - 1)
      end do

      ! store the edge wl points      
      do i = 1, nlam_cube+1
        lamb_cube(i) = (lmin - dlam_cube/2d0) + dlam_cube*(i-1)
      end do
      write(*,*) lam_cube
      write(*,*) lamb_cube


      imcube = 0d0      
      allocate (im_coord(npix))
      delpix = 2.d0*Rout*AU/real(npix)
      pixA=delpix**2
      starcorr=path2star%A/pixA
      ! create the center coordinates for the pixels
      call output("Calculate pixel coordinates ...")
      do i = 1, npix ! center coordinates of pixels in cm
        im_coord(i) = delpix/2d0 + delpix*(i - 1) - Rout*AU
      end do
            
      ! map the paths to the pixel, need to do it for each grid
      !do i = 1, ngrids
      call output("map grid "//int2string(1, '(i7)')//" to pixels")
      call map_pixels_to_path(1, npoints(1))
      !end do

      call output("Image cube has "//trim(int2string(nlam_cube, '(i7)'))//" channels (dlam="//trim(dbl2string(dlam_cube, '(F11.4)'))//") and "//trim(int2string(npix*npix, '(i7)'))//" pixels.")
      call output("Image cube size: "//trim(dbl2string(sizeof(imcube)/1.e9, '(F11.4)'))//" GB")
      call clock_write(time, utime, "Prepare image cube: ")
    end if

    lam = lmin
    ilam1Cont = 1

    do while (lam > lam_cont(ilam1Cont + 1) .and. ilam1Cont < nlamCont)
      ilam1Cont = ilam1Cont + 1
    end do

    allocate (profile(-nvprofile:nvprofile))
    allocate (profile_nz(-nvprofile:nvprofile))
    do i = 0, nR
      do j = 0, nTheta
        allocate (C(i, j)%line_abs(maxblend))
        allocate (C(i, j)%line_emis(maxblend))
      end do
    end do

    do k = -nvprofile, nvprofile
      profile(k) = ((real(k)*vresolution/vres_mult)/1d5)**2
      if (abs(profile(k)) < 10d0) then
        profile_nz(k) = .true.
      else
        profile_nz(k) = .false.
      end if
    end do
    profile(:) = exp(-profile(:))

    nl = 0

    call output("==================================================================")

    call clock(time, utime)

    call SYSTEM_CLOCK(counts, count_rate, count_max)
    starttime = DBLE(counts)/DBLE(count_rate)

    call output("Tracing "//trim(int2string(nlines, '(i6)'))//" lines")

    lcmin = lmin
    lmin_next = 0d0

    nltot = 0
    Bl => Blends
    do iblends = 1, nblends
      nltot = nltot + Bl%n
      if (iblends < nblends) Bl => Bl%next
    end do
    !call output("Total number of wl/velocity points: "//trim(int2string(nltot, '(i7)')))

    Bl => Blends
    do iblends = 1, nblends
      flux = 0d0

      nb = Bl%n

      do ib = 1, nb
        imol_blend(ib) = Bl%L(ib)%imol
        v_blend(ib) = Bl%v(ib)
      end do

      LL = Bl%L(1)
      lam = Bl%lam

      ! Only do the requested lines, lmin,lmax are set in the input file.
      if (lam > lmin .and. lam < lmax) then
        nl = nl + nb

        ! make sure that we do not end up with empty channels.
        if (imagecube) then
          ilam_cube_start = nlam_cube
          ilam_cube_end = 0

          call clock(timegap, utimegap)

          ! here lcmin is the wl that was done in the last blend, 
          ! Something is fishy here, bzt Bl%lmax is not the maximum wavelength that is actually done
          lamgapmax=lam*sqrt((1d0 + real(Bl%nvmax)*vresolution/clight)/(1d0 - real(Bl%nvmax)*vresolution/clight))
          do i = 1, nlam_cube ! check if lcmin is  not in the current cube bin            
            if ((lamb_cube(i+1) > lcmin .and. lamb_cube(i+1) < Bl%lmin) .or. &
                (iblends == nblends .and. lamb_cube(i) > lamgapmax)) then ! this is for the end
              ilam_cube_start = min(ilam_cube_start, i)
              ilam_cube_end = max(ilam_cube_end, i)
            end if
          end do
          ! I have some gap          
          if (ilam_cube_start <= ilam_cube_end) then            
            call output("  Fill imcube gap from lam "//dbl2string(lam_cube(ilam_cube_start))//" to"//dbl2string(lam_cube(ilam_cube_end))//" with continuum.")          
            !write(*,*) Bl%lmin, Bl%lmax,lamgapmax, lcmin, lamb_cube(ilam_cube_start),lam_cube(ilam_cube_start),lam_cube(ilam_cube_end)
            do i = ilam_cube_start, ilam_cube_end
              !if ((lam_cube(i) > lcmin .and. lam_cube(i) < Bl%lmin) .or. &
              !    (iblends == nblends .and. lam_cube(i) > Bl%lmax)) then ! this is for the end
                ilam_cube_gap = 1
                do while (lam_cube(i) > lam_cont(ilam_cube_gap + 1) .and. ilam_cube_gap < nlamCont)
                  ilam_cube_gap = ilam_cube_gap + 1
                end do
                call InterpolateLam(lam_cube(i), ilam_cube_gap)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(j,PP, flux0,vmult) &
!$OMP SHARED(P, npoints, lam_cube, i, igrid, ngrids)
!$OMP DO SCHEDULE(dynamic,1)
                do j = 1, npoints(1)  ! just do it for the first grid for imagecube
                  PP => P(1, j)
                  call ContContrPath(PP, flux0)
                  do vmult=-1,1,2
                    call AddImage(0, lam_cube(i), flux0*PP%A/2d0, PP, vmult, "continuum in gap")                  
                  end do
                end do
!$OMP END DO
!$OMP END PARALLEL
                ! do the star, as there are no lines, nb = 0
                PP => path2star                
                call Trace2StarLines(PP, flux0, 0, imol_blend, v_blend, 0)    
                ! for the star we correct the area as the pixel is usually much larger than the star            
                call AddImage(0, lam_cube(i), flux0*PP%A*starcorr, PP, 1, "trace star gap")
              !end if
            end do
            call clock_write(timegap, utimegap,"  Fill imcube gap: ")
          end if
        end if

        if (iblends == 1) then
          call output("  First line blend ("//trim(int2string(nb, '(i4)'))//" lines at lamcenter="//(dbl2string(Bl%lam, '(F10.3)'))//" micron)")
        end if
        

        ilam = ilam1Cont
        do while (lam > lam_cont(ilam + 1) .and. ilam < nlamCont)
          ilam = ilam + 1
          if (ilam == nlamCont) exit
        end do
        ! that writes out continuum only points for the spectrum, however, only the points available
        ! which could be only 1
        if (ilam > ilam1Cont) then
          do k = ilam1Cont + 1, ilam
            if (lam_cont(k) > lcmin .and. lam_cont(k) < Bl%lmin) then
              flux0 = 0d0
              do i = 1, ngrids
              do j = 1, npoints(i)
                PP => P(i, j)
                flux0 = flux0 + PP%flux_cont(k)*PP%A/real(ngrids)
              end do
              end do
              flux0 = flux0 + path2star%flux_cont(k)*path2star%A
              write (20, "(F15.8, 1pE16.8, 1pE16.8, 1pE16.8, 3X, A)") lam_cont(k), &
                flux0*1e23/(distance*parsec)**2, 0.0, flux0*1e23/(distance*parsec)**2, "continuum"
            end if
          end do
          ilam1Cont = ilam
        end if

        write(*,*) "Blend: ",iblends,lcmin, Bl%lmin, Bl%lmax, lam, nb,lmin_next

        lcmin = Bl%lmax
        write(*,*) "Blend: ",iblends,lcmin, Bl%lmin, Bl%lmax, lam, nb,lmin_next

        call InterpolateLam(lam, ilam)
        do i = 0, nR
          do j = 1, nTheta
            do ilines = 1, Bl%n
              LL = Bl%L(ilines)
              fact = clight*hplanck*C(i, j)%N(LL%imol)/(4d0*pi*C(i, j)%line_width(LL%imol)*sqrt(pi))
              C(i, j)%line_abs(ilines) = fact*(C(i, j)%npop(LL%imol)%N(LL%jlow)*LL%Blu - C(i, j)%npop(LL%imol)%N(LL%jup)*LL%Bul)
              C(i, j)%line_emis(ilines) = fact*C(i, j)%npop(LL%imol)%N(LL%jup)*LL%Aul
            end do
          end do
        end do

        flux2 = 0d0

        ! randomly select a grid.
        if (imagecube) then 
          i=1
        else
          i = int(ran1(idum)*real(ngrids) + 1)
        endif
        

        ! FIXME: check parallel implementation for imagecube
!$OMP PARALLEL IF(.true.) &
!$OMP DEFAULT(NONE) &
!$OMP PRIVATE(j,PP,iv,vmult,ib0,nb0,gas,ib,imol,flux0,flux4,doit,flux_c,doit_ib,doit_ib0,lam_velo,callerstr,vlo_eff,vhi_eff) &
!$OMP SHARED(i,P,ngrids,npoints,nb,LL,nv,Bl,vresolution,imol_blend,v_blend,flux,ilam) &
!$OMP SHARED(iblends,flux2,maxblend,lam,lmin_next,imagecube,nvmin,nvmax)
        allocate (doit_ib0(maxblend))
        allocate (doit_ib(maxblend))
        allocate (flux4(nvmin:nvmax))

#ifdef _OPENMP
        idum = -42 - OMP_GET_THREAD_NUM()     ! setting the seed
#endif

        flux4(:) = 0d0
!$OMP DO SCHEDULE(dynamic,1)
        do j = 1, npoints(i)
          !if (iblends == 1) call tellertje(j, npoints(i))
          PP => P(i, j)

          !write(88,*) i,j,PP%x/AU,PP%y/AU

          call ContContrPath(PP, flux_c)
!$OMP ATOMIC
          flux2 = flux2 + flux_c*PP%A
          doit = .false.
          doit_ib0 = .false.
          do ib = 1, nb
            if (PP%npopmax(Bl%L(ib)%imol) > Bl%L(ib)%jlow) then
              doit = .true.
              doit_ib0(ib) = .true.
            end if
          end do
          if (doit) then
            !write(*,*) "doit",i,j
            !write(*,*) Bl%nvmin,Bl%nvmax
            do iv = Bl%nvmin, Bl%nvmax
              lam_velo = lam*sqrt((1d0 + real(iv)*vresolution/clight)/(1d0 - real(iv)*vresolution/clight))
              if (lam_velo > lmin_next) then
                vmult = 1
                if (nb > 1) then
                  ib0 = Bl%ib0(iv)
                  nb0 = Bl%nb0(iv)
                  ! The disk is axisymmetric: each randomly-sampled path PP also represents its
                  ! mirror-image on the opposite side of the disk, which has all line-of-sight
                  ! velocities sign-flipped. We therefore trace the path twice -- once for each
                  ! side -- each weighted by half the solid-angle area (PP%A/2).
                  !   vmult = -1 : approaching (blueshifted) side, gas velocities in [-vmax, -vmin]
                  !   vmult = +1 : receding  (redshifted)  side, gas velocities in [ vmin,  vmax]
                  ! Note that to make this work the grid has to be prepared properly (i.e. only have positive y values)
                  ! Thanks to Claude
                  do vmult = -1, 1, 2
                    doit_ib(1:nb) = doit_ib0(1:nb)
                    gas = .false.
                    do ib = ib0, ib0 + nb0 - 1
                      imol = Bl%L(ib)%imol
                      ! Check whether channel iv (corrected for blend velocity offset Bl%v(ib))
                      ! overlaps the gas velocity range of this side. If not, skip line RT and
                      ! fall back to the cached continuum flux.
                      if (vmult == -1) then
                        vlo_eff = -PP%vmax(imol)   ! approaching side: velocity range is [-vmax, -vmin]
                        vhi_eff = -PP%vmin(imol)
                      else
                        vlo_eff = PP%vmin(imol)    ! receding side: velocity range is [vmin, vmax]
                        vhi_eff = PP%vmax(imol)
                      end if
                      if (((real(iv) - 0.5d0)*vresolution - Bl%v(ib)) < vhi_eff &
                          .and. ((real(iv) + 0.5d0)*vresolution - Bl%v(ib)) > vlo_eff) then
                        gas = .true.
                        exit
                      end if
                    end do
                    if (gas) then
                      call TraceFluxLines(PP, flux0, iv, vmult, imol_blend(ib0), v_blend(ib0), doit_ib(ib0), nb0, ib0)
                      if (imagecube) call AddImage(iv, lam, flux0*PP%A/2d0, PP, vmult, callerstr)
                    else
                      flux0 = flux_c
                      !write (callerstr, "(A,' ' ,i5,A,F10.3,A,F10.3)") "call nb ", iblends,"lmin:",lmin_next,"lv:",lam_velo
                      ! FIXME: It cann happen hat somehoe a iv of several blends produces a continuum in the same imcube channel
                      ! hence the continuum is contouned multipe times, the reason is the if (gas) thingy  
                      if (imagecube) call AddImage(iv, lam, flux0*PP%A/2d0, PP, vmult, "cont")
                    end if
                    flux4(iv) = flux4(iv) + flux0*PP%A/2d0
                    !write (callerstr, "(A,' ' ,i5,A,F10.3,A,F10.3)") "call nb ", iblends,"lmin:",lmin_next,"lv:",lam_velo
                    !if (imagecube) call AddImage(iv, lam, flux0*PP%A/2d0, PP, vmult, callerstr)
                  end do
                  ! channel fully outside the line, so just add the continuum contribution
                else if (((real(iv) + 0.5d0)*vresolution > PP%vmax(LL%imol) .and. &
                          (real(iv) - 0.5d0)*vresolution > PP%vmax(LL%imol)) &
                         .or. ((real(iv) + 0.5d0)*vresolution < PP%vmin(LL%imol) .and. &
                               (real(iv) - 0.5d0)*vresolution < PP%vmin(LL%imol))) then
                  flux0 = flux_c
                  do vmult = -1, 1, 2
                    flux4(iv*vmult) = flux4(iv*vmult) + flux0*PP%A/2d0
                    ! Seems to be required here and is done twice
                    ! FIXME: callerstring is only for debugging, remove it
                    write (callerstr, "(A,' ' ,i5,' ',i2)") "call nb else if", iv, vmult
                    if (imagecube) then
                      call AddImage(iv*vmult, lam, flux0*PP%A/2d0, PP, vmult, "cont")
                    end if
                  end do
                  ! channel covers single line (not blended)
                else
                  doit_ib = .true.
                  call TraceFluxLines(PP, flux0, iv, vmult, imol_blend, v_blend, doit_ib, nb, 1)
                  do vmult = -1, 1, 2
                    flux4(iv*vmult) = flux4(iv*vmult) + flux0*PP%A/2d0
                    write (callerstr, "(A,' ' ,i5,' ',i2)") "call nb else", iv, vmult
                    if (imagecube) then
                      call AddImage(iv*vmult, lam, flux0*PP%A/2d0, PP, vmult, callerstr)
                    end if
                  end do
                end if
              end if
            end do
          else
            !write(*,*) "call not doit", not populated
            ! not sure why here we do not use vmult treatment
            flux0 = flux_c
            flux4(Bl%nvmin:Bl%nvmax) = flux4(Bl%nvmin:Bl%nvmax) + flux0*PP%A
            if (imagecube.and..false.) then ! do the vmult treatment for the image cube, important where the flux is
              do vmult = -1, 1, 2
                do iv = Bl%nvmin, Bl%nvmax
                  lam_velo = lam*sqrt((1d0 + real(iv)*vresolution/clight)/(1d0 - real(iv)*vresolution/clight))
                  if (lam_velo > lmin_next) then
                    call AddImage(iv, lam, flux0*PP%A, PP, vmult, "call not doit")
                  endif
                end do
              end do 
            end if
          end if
        end do
!$OMP END DO
!$OMP CRITICAL
        flux(:) = flux(:) + flux4(:)
!$OMP END CRITICAL
        deallocate (doit_ib0)
        deallocate (doit_ib)
        deallocate (flux4)
!$OMP FLUSH
!$OMP END PARALLEL
        PP => path2star
        ! do the star
        if (.not. allocated(PP%im_ixy)) allocate (PP%im_ixy(2, 1))
        do iv = Bl%nvmin, Bl%nvmax
          if (nb > 1) then  ! FIXME: this if seems to be unnecessary
            call Trace2StarLines(PP, flux0, iv, imol_blend, v_blend, nb)
            flux(iv) = flux(iv) + flux0*PP%A
          else
            call Trace2StarLines(PP, flux0, iv, imol_blend, v_blend, nb)
            flux(iv) = flux(iv) + flux0*PP%A
          end if
          if (imagecube) call AddImage(iv, lam, flux0*PP%A*starcorr, PP, 1, "trace star")
        end do
        flux2 = flux2 + flux0*PP%A

        flux1 = 0d0
        flux3 = 0d0

        ! f = sqrt((1d0 + real(Bl%nvmin)*vresolution/clight)/(1d0 - real(Bl%nvmin)*vresolution/clight))
        ! wl11 = log10(lam_cont(ilam + 1)/(lam*f))/log10(lam_cont(ilam + 1)/lam_cont(ilam))
        ! ! wl11 can be > 1 if nv*vresolution is wider than the spacing between the continuum points.
        ! ! So the continuum grid is too fine, just assume the nearest continuum point then.
        ! ! > 1 would mean extrapolation, but in some cases that ends in NaN, limiting it to 1 means, no extrapolation
        ! wl11 = min(1d0, max(0d0, wl11))
        ! wl21 = 1d0 - wl11

        ! f = sqrt((1d0 + real(Bl%nvmax)*vresolution/clight)/(1d0 - real(Bl%nvmax)*vresolution/clight))
        ! wl13 = log10(lam_cont(ilam + 1)/(lam*f))/log10(lam_cont(ilam + 1)/lam_cont(ilam))
        ! ! see above wl11
        ! wl13 = min(1d0, max(0d0, wl13))
        ! wl23 = 1d0 - wl13

        ! flux_l1 = 0d0
        ! flux_l2 = 0d0
        ! do k = 1, ngrids
        !   do j = 1, npoints(k)
        !     PP => P(k, j)
        !     flux_l1 = flux_l1 + PP%flux_cont(ilam)*PP%A/real(ngrids)
        !     flux_l2 = flux_l2 + PP%flux_cont(ilam + 1)*PP%A/real(ngrids)
        !   end do
        ! end do
        ! PP => path2star

        ! flux_l1 = flux_l1 + PP%flux_cont(ilam)*PP%A
        ! flux_l2 = flux_l2 + PP%flux_cont(ilam + 1)*PP%A

        ! flux1 = flux_l1**wl11*flux_l2**wl21
        ! flux3 = flux_l1**wl13*flux_l2**wl23

        ! do iv = Bl%nvmin, Bl%nvmax
        !   fc = flux1 + (flux3 - flux1)*real(iv - Bl%nvmin)/real(Bl%nvmax - Bl%nvmin)
        !   flux(iv) = flux(iv) - flux2 + fc
        !   flux_cont(iv) = fc
        ! end do

        if (Bl%n == 1) then
          comment = trim(Mol(Bl%L(1)%imol)%name)//"  up: "//trim(int2string(Bl%L(1)%jup, '(i5)')) &
                    //"  low: "//trim(int2string(Bl%L(1)%jlow, '(i5)'))
        else
          count = 0
          do iv = 1, Bl%n
            count(Bl%L(iv)%imol) = count(Bl%L(iv)%imol) + 1
          end do
          comment = "blend of "
          do iv = 1, nmol
            if (count(iv) > 0) then
              comment = trim(comment)//trim(int2string(count(iv), '(i3)'))//" "//trim(Mol(iv)%name)
            end if
          end do

        end if

        Bl%F = 0d0
        lam_velo = lam*sqrt((1d0 + vresolution/clight)/(1d0 - vresolution/clight))
        dnu = dabs(clight*1d4*(1d0/lam_velo - 1d0/lam))
        lam_w = 0d0
        lam_w_min = 1d200
        lam_w_max = 0d0
        do iv = Bl%nvmin, Bl%nvmax
          lam_velo = lam*sqrt((1d0 + real(iv)*vresolution/clight)/(1d0 - real(iv)*vresolution/clight))
          if (lam_velo > lmin_next) then
            write (20, "(F15.8, 1pE16.8, 1pE16.8, 1pE16.8, 3X, A)") lam_velo, &
              flux(iv)*1e23/(distance*parsec)**2, &
              real(iv)*vresolution/1d5, &
              flux_cont(iv)*1e23/(distance*parsec)**2, &
              trim(comment)
            Bl%F = Bl%F + dnu*(flux(iv) - flux_cont(iv))
            lam_w = lam_w + lam_velo*dnu*(flux(iv) - flux_cont(iv))
            if (lam_velo < lam_w_min) lam_w_min = lam_velo
            if (lam_velo > lam_w_max) lam_w_max = lam_velo
          end if
        end do
        if (Bl%F /= 0d0) then     ! added PW, June 13, 2022
          lam_w = lam_w/Bl%F
        end if
        Bl%F = Bl%F*1d-3/(distance*parsec)**2
        write (21, *) lam_w, Bl%F, lam_w_min, lam_w_max, trim(comment)

        lmin_next = max(lmin_next, lam_velo)

      end if ! if(lam>lmin.and.lam<lmax)

      if (iblends < nblends) Bl => Bl%next

      call tellertje_time(iblends, nblends, nl, nltot, starttime)

    end do

    if (imagecube) then
      call output("")
      call output("Writing image cube to file imcube.fits ...")
      call writefitsfile(trim(imagecube_filename), imcube*1e23/(distance*parsec)**2, nlam_cube, npix, lmin, lmin, .true.)
      call output("")      
    end if

    ilam = ilam + 1
    do while (ilam <= nlamCont)
      if (lam_cont(ilam) < lmax) then
        flux0 = 0d0
        do j = 1, npoints(i)
          PP => P(i, j)
          flux0 = flux0 + PP%flux_cont(ilam)*PP%A
        end do
        flux0 = flux0 + path2star%flux_cont(ilam)*path2star%A
        write (20, "(F15.8, 1pE16.8)") lam_cont(ilam), flux0*1e23/(distance*parsec)**2
      end if
      ilam = ilam + 1
    end do

    close (unit=20)
    close (unit=21)

    !call cpu_time(stoptime)
    call SYSTEM_CLOCK(counts, count_rate, count_max)
    stoptime = DBLE(counts)/DBLE(count_rate)

    call output("Time used for the lines:"//trim(dbl2string(stoptime - starttime, '(f8.2)')) &
                //" s")
    call output("Time used per line:     "//trim(dbl2string((stoptime - starttime)/real(nlines), '(f8.2)')) &
                //" s")
    !   call output("Time used per line:     "//trim(dbl2string((stoptime-starttime)/real(nblends),'(f8.2)'))
    !   //" s")
    call clock_write(time, utime, "RaytraceLines:")

    return
  end subroutine RaytraceLines

  subroutine ContContrPath(p0, flux)
!   Compute the continuum flux along path p0 at the current line-center
!   wavelength (set by the last call to InterpolateLam). Integrates the
!   formal RT solution cell by cell using the pre-interpolated dust opacity
!   kext_l and source function (therm_l + scat_l). As a side effect, caches
!   exptau_dust(k), S_dust(k) and cont_contr(k) per path element so that
!   TraceFluxLines can reuse them without recomputing the dust RT.
! DOC by Claude
    use GlobalSetup
    use Constants
    implicit none

    type(Path), intent(inout) :: p0
    real(kind=8), intent(out) :: flux

    integer :: i, j, k
    real(kind=8) :: tau, fact, S, tau_dust, tau_d, tau_tot

    type(Cell), pointer :: CC

    fact = 1d0
    flux = 0d0
    tau_tot = 0d0

    do k = 1, p0%n
      i = p0%i(k)
      j = p0%j(k)
      if (i > 0 .and. i < nR .and. j > 0) then
        CC => C(i, j)
        tau_dust = CC%kext_l
        ! dust thermal source function
        S = CC%therm_l*tau_dust
        ! dust scattering source function
        S = S + CC%scat_l*tau_dust

        p0%S_dust(k) = S

        tau = tau_dust
        tau_d = tau*p0%d(k)
        if (tau_d > 1d-4) then
          p0%exptau_dust(k) = exp(-tau_d)
          p0%cont_contr(k) = S*(1d0 - p0%exptau_dust(k))/tau
        else
          p0%exptau_dust(k) = 1d0 - tau_d
          p0%cont_contr(k) = S*p0%d(k)
        end if
        flux = flux + p0%cont_contr(k)*fact

        fact = fact*p0%exptau_dust(k)
        tau_tot = tau_tot + tau_d
        if (tau_tot > tau_max) return
      end if
    end do

    return
  end subroutine ContContrPath

  subroutine TraceFluxLines(p0, flux, ii, vmult, imol_blend, v_blend, doit, nb, ib0)
    use GlobalSetup
    use Constants
    implicit none

    type(Path), intent(in) :: p0
    real(kind=8), intent(out) :: flux
    integer, intent(in) :: ii, vmult, imol_blend(nb)
    real(kind=8), intent(in) :: v_blend(nb)
    logical, intent(in) :: doit(nb)
    integer, intent(in) :: nb, ib0

    integer :: i, j, k, imol
    real(kind=8) :: tau, exptau, fact, prof, S, tau_gas, tau_dust, tau_d, tau_tot
    double precision :: v
    integer :: jj, ib
    type(Cell), pointer :: CC
    logical :: gas
    real(kind=8) :: rj

    fact = 1d0
    flux = 0d0
    tau_tot = 0d0

    v = real(ii) + ran1(idum) - 0.5d0

    do k = 1, p0%n
      i = p0%i(k)
      j = p0%j(k)
      if (i > 0 .and. i < nR .and. j > 0) then
        CC => C(i, j)

        tau = 0d0
        S = 0d0
        gas = .false.

        do ib = 1, nb
          if (doit(ib)) then
            imol = imol_blend(ib)
            rj = ((real(vmult)*p0%v(k) + v_blend(ib))*vres_mult/vresolution - v*vres_mult) &
                 *1d5/CC%line_width(imol)
            jj = int(rj)
            if (jj > -nvprofile .and. jj < nvprofile) then
              if (profile_nz(jj)) then
                prof = profile(jj)
                tau_gas = prof*CC%line_abs(ib + ib0 - 1)
                ! gas source function
                S = S + prof*CC%line_emis(ib + ib0 - 1)
                tau = tau + tau_gas
                gas = .true.
              end if
            end if
          end if
        end do

        if (gas) then
          tau_dust = CC%kext_l
          S = S + p0%S_dust(k)

          tau = tau + tau_dust

          tau_d = tau*p0%d(k)
          if (tau_d > 1d-4) then
            exptau = exp(-tau_d)
            flux = flux + S*(1d0 - exptau)*fact/tau
          else
            exptau = 1d0 - tau_d
            flux = flux + S*p0%d(k)*fact
          end if
        else
          tau_d = CC%kext_l*p0%d(k)
          exptau = p0%exptau_dust(k)
          flux = flux + p0%cont_contr(k)*fact
        end if

        fact = fact*exptau
        tau_tot = tau_tot + tau_d

        if (tau_tot > tau_max) return
      end if
    end do

    return
  end subroutine TraceFluxLines

  subroutine Trace2StarLines(p0, flux, ii, imol_blend, v_blend, nb)
    use GlobalSetup
    use Constants
    implicit none
    type(Path), intent(in) :: p0
    real(kind=8), intent(out) :: flux
    integer, intent(in) :: ii, imol_blend(nb)
    real(kind=8), intent(in) :: v_blend(nb)
    integer, intent(in) :: nb

    integer :: i, j, k, imol
    real(kind=8) :: tau, prof
    integer :: ib, jj
    type(Cell), pointer :: CC

    tau = 0d0

    do k = 1, p0%n
      i = p0%i(k)
      if (i == 0) exit
      j = p0%j(k)
      if (i /= 0 .and. i /= nR .and. j /= 0) then
        CC => C(i, j)
        tau = tau + CC%kext_l

        do ib = 1, nb
          imol = imol_blend(ib)
          jj = int(((p0%v(k) + v_blend(ib))*vres_mult/vresolution - real(ii)*vres_mult)*1d5/CC%line_width(imol))
          if (jj < -nvprofile) jj = -nvprofile
          if (jj > nvprofile) jj = nvprofile
          prof = profile(jj)
          tau = tau + prof*CC%line_abs(ib)
        end do
      end if
    end do

    flux = Fstar_l*exp(-tau)

    return
  end subroutine Trace2StarLines

  subroutine DetermineBlends(nv, maxblend, nvmin0, nvmax0)
    use GlobalSetup
    use Constants
    use InOut
    implicit none

    integer, intent(in) :: nv
    integer, intent(out) :: maxblend, nvmin0, nvmax0

    integer :: i, iv, j
    real(kind=8) :: maxvshift, maxmult, v, f
    type(Blend), pointer :: Bl
    real(kind=8) :: time,utime

    call clock(time,utime)
    call output("")
    call output("Determining blends...")        

    maxvshift = 2d0*real(nv)*vresolution
    maxmult = sqrt((1d0 + maxvshift/clight)/(1d0 - maxvshift/clight))

    Bl => Blends

    nblends = 0
    maxblend = 0
    nvmin0 = -nv
    nvmax0 = nv

    do i = 1, nlines
      Bl%n = 0
      Bl%nvmin = -nv
      Bl%nvmax = nv
      Bl%lam = Lines(i)%lam
      ! this just determines the number of lines in the blend.
      do j = 1, nlines
        if ((j > i .and. (Lines(j)%lam/Lines(i)%lam) < maxmult) .or. &
            (j < i .and. (Lines(i)%lam/Lines(j)%lam) < maxmult) .or. j == i) then
          Bl%n = Bl%n + 1
        end if
      end do
      allocate (Bl%L(Bl%n))
      allocate (Bl%v(Bl%n))
      if (Bl%n > maxblend) maxblend = Bl%n
      ! now actually assign the lines to the blend
      Bl%n = 0
      do j = 1, nlines
        if ((j > i .and. (Lines(j)%lam/Lines(i)%lam) < maxmult) .or. &
            (j < i .and. (Lines(i)%lam/Lines(j)%lam) < maxmult) .or. j == i) then
          Bl%n = Bl%n + 1
          Bl%L(Bl%n) = Lines(j)
          if (i == j) then
            Bl%v(Bl%n) = 0d0 ! so the current line i, is the center of the blend
          else
            f = (Lines(j)%lam/Lines(i)%lam)**2
            Bl%v(Bl%n) = clight*(f - 1d0)/(f + 1d0) ! this is the velocity fot line j on the vel grid of the current line i which sits a v=0
          end if
        end if
      end do
      allocate (Bl%ib0(-nv:nv))
      allocate (Bl%nb0(-nv:nv))
      do iv = -nv, nv
        do j = 1, Bl%n
          if (real(iv - nv)*vresolution <= Bl%v(j)) exit
        end do
        Bl%ib0(iv) = j
        do j = Bl%n, 1, -1
          if (real(iv + nv)*vresolution >= Bl%v(j)) exit
        end do
        Bl%nb0(iv) = j - Bl%ib0(iv) + 1
      end do

      v = -real(nv)*vresolution
      Bl%lmin = Bl%L(1)%lam*sqrt((1d0 + v/clight)/(1d0 - v/clight))
      v = real(nv)*vresolution
      Bl%lmax = Bl%L(1)%lam*sqrt((1d0 + v/clight)/(1d0 - v/clight))

      allocate (Bl%next)
      Bl => Bl%next
      nblends = nblends + 1
    end do

    call output("Number of lines: "//trim(int2string(nlines, '(i7)'))//"  Number of blends: "//trim(int2string(nblends, '(i7)')))
    call clock_write(time, utime, "DetermineBlends: ")

    return
  end subroutine DetermineBlends

  subroutine InterpolateLam(lam0, ilam)
    use GlobalSetup
    use Constants
    use InOut
    implicit none
    real(kind=8), intent(in) :: lam0
    integer, intent(in) :: ilam

    real(kind=8) :: wl1, wl2, w1, w2, x1, x2
    integer :: i, j

    if (lam0 < lam_cont(ilam) .or. lam0 > lam_cont(ilam + 1)) then
      call output("Error: InterpolateLam called with lam0 outside the range of the continuum grid points. This should not happen. Check the code.")
      stop
    end if
    w1 = (lam_cont(ilam + 1) - lam0)/(lam_cont(ilam + 1) - lam_cont(ilam))
    w2 = 1d0 - w1
    wl1 = log10(lam_cont(ilam + 1)/lam0)/log10(lam_cont(ilam + 1)/lam_cont(ilam))
    wl2 = 1d0 - wl1
    do i = 0, nR
      do j = 1, nTheta
        x1 = C(i, j)%kext(ilam)
        x2 = C(i, j)%kext(ilam + 1)
        C(i, j)%kext_l = (x1*w1) + (x2*w2)
        !   x1=C(i,j)%LRF(ilam)*C(i,j)%albedo(ilam)*C(i,j)%kext(ilam)
        !   x2=C(i,j)%LRF(ilam+1)*C(i,j)%albedo(ilam+1)*C(i,j)%kext(ilam+1)
        !   C(i,j)%scat_l=(x1*w1+x2*w2)/C(i,j)%kext_l
        !   x1=BB(ilam,C(i,j)%iT)
        !   x2=BB(ilam+1,C(i,j)%iT)
        !   C(i,j)%therm_l=w1*x1+w2*x2
        !   x1=(1d0-C(i,j)%albedo(ilam))*C(i,j)%kext(ilam)
        !   x2=(1d0-C(i,j)%albedo(ilam+1))*C(i,j)%kext(ilam+1)
        !   C(i,j)%therm_l=C(i,j)%therm_l*(w1*x1+w2*x2)/C(i,j)%kext_l

        ! Using the source function for now.
        x1 = C(i, j)%S(ilam)
        x2 = C(i, j)%S(ilam + 1)
        C(i, j)%therm_l = (w1*x1 + w2*x2)/2d0
        C(i, j)%scat_l = (w1*x1 + w2*x2)/2d0
      end do
    end do
    Fstar_l = (Fstar(ilam)**wl1)*(Fstar(ilam + 1)**wl2)

    return
  end subroutine InterpolateLam

  subroutine map_pixels_to_path(igrid, npath)
    use GlobalSetup
    use Constants
    use InOut
    implicit none
    integer, intent(in) :: igrid, npath
    real(kind=8) :: px(npath), py(npath), dist(npath), im_x, im_y
    real(kind=8) :: maxx, maxr, maxrxp, maxrxm !pixA
    type(Path), pointer :: p0
    real(kind=8) :: time, utime

    integer :: i, j, ipath, iminpath, maxnpixpath

    call clock(time, utime)

    ! If already allocated this was done already (same grid igrid is used again), no need to map things again
    if (allocated(P(igrid, 1)%im_ixy)) return

    !call output("Mapping pixels to path igrid,npath:"//trim(int2string(igrid, '(i7)'))//" "//trim(int2string(npath, '(i7)')))

    ! area of one pixel cm^2
    !pixA = (im_coord(2) - im_coord(1))**2

    ! the paths do not sample the whole image (just the disk)
    ! In theory this could be different (e.g. if not a disk)
    !maxr=maxval(abs(P(igrid,:)%y)) ! this should always be the max r, major axis
    !write(*,*) maxr/AU
    maxr = R(nr)
    !maxrxp=maxval(P(igrid,:)%x)
    !maxrxm=minval(P(igrid,:)%x)
    !maxr=max(maxrxp,maxrxm)

    P(igrid, :)%im_npix = 0
    px(1:npath) = P(igrid, 1:npath)%x
    py(1:npath) = P(igrid, 1:npath)%y

    maxnpixpath = npix*int(npix/10)
    do ipath = 1, npath
      ! FIXME: think about how large that can be
      allocate (P(igrid, ipath)%im_ixy(2, maxnpixpath))
    end do

    do i = 1, npix
      call tellertje(i,npix)
      im_x = im_coord(i)
      ! only care about the positive half
      do j = int(npix/2 + 1), npix
        im_y = im_coord(j)
        ! the Path grid does not cover a circle (due to inclination)
        ! along the x is the minor axis, which has different Rmax depending on sign
        if (sqrt((im_x)**2 + (im_y)**2) > maxr) cycle
        !if ((im_x<0).and.im_x<maxrxm) cycle FIXME: doesn't really work
        !if ((im_x>0).and.im_x>maxrxp) cycle
        dist = sqrt((px - im_x)**2 + (py - im_y)**2)
        ! don't need squareroot
        iminpath = minloc(dist, dim=1)
!                        if (i==int(npix/2)) then
!                                write(*,*) "out: ",i,j,im_x/AU,im_y/AU,dist(iminpath)/AU,iminpath,P(igrid,iminpath)%x/AU,P(igrid,iminpath)%y/AU
!                        end if

        P(igrid, iminpath)%im_npix = P(igrid, iminpath)%im_npix + 1

        if (P(igrid, iminpath)%im_npix > maxnpixpath) then
          call output('Problem found to many pixels for path')
          call output(dbl2string(P(igrid, iminpath)%x/AU) // ' ' // dbl2string(P(igrid, iminpath)%y/AU) // ' ' // int2string(int(P(igrid, iminpath)%im_npix,kind=4), '(i7)'))
          stop
        end if
        P(igrid, iminpath)%im_ixy(1, P(igrid, iminpath)%im_npix) = i
        P(igrid, iminpath)%im_ixy(2, P(igrid, iminpath)%im_npix) = j
      end do
    end do

    ! now there will be paths with no pixel assigned as in case for Paths being
    ! close to each other the pixel is only assigned to the closest one
    ! so simply go through the paths without pixels, and assign the closest pixel .

    do ipath = 1, npath
      p0 => P(igrid, ipath)
      if (p0%im_npix > 0) cycle
      ! find the closest pixel
      p0%im_npix = 1
      p0%im_ixy(1, 1) = minloc(abs(im_coord - p0%x), dim=1)
      p0%im_ixy(2, 1) = minloc(abs(im_coord - p0%y), dim=1)
    end do

    ! special treatment for the star
    ! assume here it is always just one pixel for the moment
    ! and only needs to be done once
    if (.not. allocated(path2star%im_ixy)) then
      allocate (path2star%im_ixy(2, 1))

      path2star%im_npix = 1
      path2star%im_ixy(1, 1) = minloc(abs(im_coord - path2star%x), dim=1)
      path2star%im_ixy(2, 1) = minloc(abs(im_coord - path2star%y), dim=1)
      !write (*, *) path2star%x, path2star%y, path2star%im_ixy(:, 1)
    end if

    call clock_write(time, utime, "map_pixels_to_path:")
    ! some log output and set the correction factor

    ! p0 => path2star
    ! write (99, *) p0%x/AU, p0%y/AU, p0%A/AU/AU, p0%im_ixy(:, 1), p0%im_npix, p0%im_npix*pixA/AU/AU

    ! do ipath = 1, npath
    !   p0 => P(igrid, ipath)
    !   write (99, *) p0%x/AU, p0%y/AU, p0%A/AU/AU, p0%im_ixy(:, 1), p0%im_npix, p0%im_npix*pixA/AU/AU
    ! end do

  end subroutine map_pixels_to_path

  subroutine AddImage(iv, lam_blend, flux0, p0, vmult, caller)
    ! Maps velocity channel iv (relative to lam_blend) to the global linear
    ! wavelength grid and accumulates flux into imcube.
    use GlobalSetup
    use Constants
    use InOut
    implicit none
    integer, intent(in) :: iv
    real(kind=8), intent(in) :: lam_blend
    real(kind=8), intent(in) :: flux0
    integer, intent(in) :: vmult
    type(Path), intent(in) :: p0
    character(len=*), intent(in) :: caller        ! just for logging, can be removed
    integer :: ix, iy, ipix, ivcube
    real(kind=8) :: fluxperpix, lam_velo

    ! the cube goes strictly from lmin (channel centers) to lmax, but parts of a bleand can be outside    
    ! compute absolute wavelength for this channel and map to global cube index
    lam_velo = lam_blend*sqrt((1d0 + real(iv)*vresolution/clight)/ &
                              (1d0 - real(iv)*vresolution/clight))
    if (lam_velo<lamb_cube(1) .or. lam_velo>lamb_cube(nlam_cube+1)) then
      ! this can happen for the channels at the edge of the blend, just ignore those.
      return
    end if

    ivcube = idnint((lam_velo - lam_cube(1))/dlam_cube) + 1
    ! Sanity check
    if (lam_velo<lamb_cube(ivcube).or.lam_velo>lamb_cube(ivcube+1)) then
      write(*,*) ivcube,lam_velo, lamb_cube(ivcube), lam_cube(ivcube),lamb_cube(ivcube+1)
      call output("Warning: velocity channel iv="//trim(int2string(iv, '(i5)'))// &
                  " at lam_velo="//trim(dbl2string(lam_velo, '(f15.8)'))// &
                  " is outside of the cube limits. Skipping this channel for the image cube. Caller: "//trim(caller))
      stop
    end if
    !write(*,*) caller," ", lam_blend,lam_velo,lmin,ivcube
    ! some parts of the line can be outside of the cube limits, we ignore those.
    if (ivcube < 1 .or. ivcube > nlam_cube) return

    
    fluxperpix = flux0/p0%im_npix    

    do ipix = 1, p0%im_npix

      ix = p0%im_ixy(1, ipix)
      iy = p0%im_ixy(2, ipix) 

      ! assume here x and y = 0 at the center of the image
      if (vmult < 0) then
        iy = npix - (iy - 1) ! for mirroring the whole thing
      end if

      imcube(iy, ix, ivcube) = imcube(iy, ix, ivcube) + fluxperpix ! P%y (Vertical) seems to be along the major axis - make it x
    end do

    ! if (ivcube==57) then
    !   !write(*,*) lamb_cube(57),lam_cube(57),lamb_cube(58)
    !   write(*,*) "caller: ", trim(caller), " iv: ", iv, " lam_blend: ", lam_blend, " lam_velo: ", lam_velo, " ivcube: ", ivcube, vmult,p0%im_npix,fluxperpix,p0%x,p0%y,ix,iy      
    ! end if  

    return
  end subroutine AddImage

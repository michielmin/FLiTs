subroutine write_imcube(filename, specunitwl)
  ! Write the image cube im, to a fits files with name filename.
  !
  ! Please note the flux values are stored in single precision to save disks space.  
  !
  ! TODO: activate the option for spectral units in frequency, but currently only wavelength is supported.
  !
  use GlobalSetup, only: distance, vresolution,nlam_cube,npix,dlam_cube,inc,delpix,imcube,lam_cube
  use Constants, only: pi, clight, parsec, AU
  use InOut
  IMPLICIT NONE
  character(len=*),intent(in) :: filename
  logical, intent(in) :: specunitwl
  real(kind=8) :: centerpix

  real(kind=8) :: degree, spixdegree
  character(len=500) :: VersionGIT

  integer status, unit, blocksize, bitpix, naxis
  integer group
  integer :: naxes(3)
  integer(kind=8) :: nelements, fpixel
  logical truefalse
  real(kind=8) :: refFrqHz, delHz, restFrqHz,reflam

  inquire (file=filename, exist=truefalse)
  if (truefalse) then
    call output("FITS file already exists, overwriting")
    open (unit=90, file=filename)
    close (unit=90, status='delete')
  end if
  
  status = 0
  ! Get an unused Logical Unit Number to use to create the FITS file
  call ftgiou(unit, status)
  ! create the new empty FITS file
  blocksize = 1
  call ftinit(unit, filename, blocksize, status)

  ! initialize parameters about the FITS image (IMDIM x IMDIM 64-bit reals)
  bitpix = -32 ! single precision should be enough
  naxes(1) = npix
  naxes(2) = npix
  if (nlam_cube > 1) then
    naxis = 3
    naxes(3) = nlam_cube
  else
    naxis = 2
    naxes(3) = 1
  end if

  ! write the required header keywords
  call FTPHPS(unit,bitpix, naxis, naxes, status)  

  degree = pi/180.0           ! one degree in rad
  !----- constant size of all pixels in [sr] -----
  ! dist is from Parameter.in (is already converted to cm)  
  spixdegree = delpix/(distance*parsec)/degree
  ! this is the reference wavelength, which is the first channel in the cube. This is used for the spectral axis keywords in the FITS header.
  reflam=lam_cube(1) 
  ! also put the restFrquency in Hz, for reference
  restFrqHz = clight/(reflam/1.e4)

  ! put some coordinate system
  centerpix = npix/2 + 1.0 ! require odd number of pixels
  call ftpkys(unit, 'ctype1', 'RA---SIN', '', status)
  call ftpkys(unit, 'cunit1', 'deg     ', '', status)
  call ftpkyd(unit, 'crval1', 68.0, 12, 'dummy value', status)
  call ftpkyd(unit, 'cdelt1', spixdegree, 12, '', status)
  call ftpkyd(unit, 'crota1', 0.0, 12, '', status)
  call ftpkyd(unit, 'crpix1', centerpix, 12, 'center pixel', status)

  call ftpkys(unit, 'ctype2', 'DEC--SIN', '', status)
  call ftpkys(unit, 'cunit2', 'deg     ', '', status)
  call ftpkyd(unit, 'crval2', 24.0, 12, 'dummy value', status)
  call ftpkyd(unit, 'cdelt2', spixdegree, 12, '', status)
  call ftpkyd(unit, 'crota2', 0.0, 12, '', status)
  call ftpkyd(unit, 'crpix2', centerpix, 12, 'center pixel', status)
  call ftpkys(unit, 'RADESYS', 'ICRS', '', status)
      
  if (specunitwl) then
    call ftpkys(unit, 'CTYPE3', 'WAVE', '[micron]', status)
    call ftpkys(unit, 'cunit3', 'um      ', '', status)
    call ftpkyd(unit, 'crval3', reflam, 12, 'Reference wl', status)
    call ftpkyd(unit, 'CDELT3', dlam_cube, 12, 'd_wl (channel width)', status)
    call ftpkyd(unit, 'crpix3', 1.0, 12, 'Reference channel', status)
  else
    ! FIXME: not really supported at the momement
    refFrqHz = clight/(reflam/1.e4)
    delHz = -restFrqHz*vresolution/clight
    !----- spectral axis -----
    ! Frequency/Spectral coordinate
    ! the first channel (velo) point is the reference
    call ftpkys(unit, 'CTYPE3', 'FREQ', '[Hz]', status)
    call ftpkys(unit, 'cunit3', 'Hz      ', '', status)
    call ftpkyd(unit, 'crval3', refFrqHz, 12, 'Reference FREQ', status)
    call ftpkyd(unit, 'CDELT3', delHz, 12, 'd_Hz (channel width)', status)
    call ftpkyd(unit, 'crpix3', 1.0, 12, 'Reference channel', status)
  end if

  call ftpkys(unit, 'SPECSYS', 'LSRK', 'Spectral reference frame', status)
  call ftpkys(unit, 'timesys', 'UTC     ', '', status)  
  call ftpkys(unit, 'BUNIT', 'Jy/pixel', 'Flux density per pixel', status)  
  ! Restfrequency for the first spectral channel 
  call ftpkyd(unit, 'RESTFREQ', restFrqHz, 12, '', status) 

  call ftpkyd(unit, 'incl', inc, 7, 'Inclination [deg]', status)
  call ftpkyd(unit, 'dist', distance, 7, 'Distance [pc]', status)
  call ftpkys(unit, 'origin', 'FLiTs/ProDiMo model', '', status)
  call ftpkys(unit, 'version', trim(VersionGIT()), 'FLiTs version', status)
  
  if (status /= 0) then
    call output("Error writing header to FITS file: "//trim(filename)//". Status: "//int2string(status))
  end if

  ! write the array to the FITS file
  group = 1
  fpixel = 1_8
  nelements = int(naxes(1),kind=8)*int(naxes(2),kind=8)*int(naxes(3),kind=8)
  ! write and convert to proper units
  call ftppre(unit, group, fpixel, nelements, imcube(:,:,:)*real(1e23/(distance*parsec)**2,kind=4), status)
  if (status /= 0) then
    call output("Error writing data to FITS file: "//trim(filename)//". Status: "//int2string(status))    
  end if

  ! close the file and free the unit number
  call ftclos(unit, status)
  call ftfiou(unit, status)

  return
end subroutine write_imcube

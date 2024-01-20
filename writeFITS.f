	subroutine writefitsfile(filename,im,nlam,n)
	use GlobalSetup, only : Rout,distance,lmin,vresolution
	use Constants, only : pi,clight,parsec,AU
	IMPLICIT NONE
	character*500 filename
	integer,intent(in) :: nlam,n
	real centerpix
	real*8 im(n,n,nlam)
	real*8 degree,delpix,spixdegree

	integer status,unit,blocksize,bitpix,naxis,naxes(3)
	integer i,j,group,fpixel,nelements
	logical simple,extend,truefalse
	real refFrqHz,delHz

	inquire(file=filename,exist=truefalse)
	if(truefalse) then
		call output("FITS file already exists, overwriting")
		open(unit=90,file=filename)
		close(unit=90,status='delete')
	endif

	status=0
C   Get an unused Logical Unit Number to use to create the FITS file
	call ftgiou(unit,status)
C   create the new empty FITS file
	blocksize=1
	call ftinit(unit,filename,blocksize,status)

C     initialize parameters about the FITS image (IMDIM x IMDIM 64-bit reals)
	simple=.true.
	!bitpix=-64
	bitpix=-32 ! single precision should be enough
	naxes(1)=n
	naxes(2)=n
	if(nlam.gt.1) then
		naxis=3
		naxes(3)=nlam
	else
		naxis=2
		naxes(3)=1
	endif
	extend=.true.

C     write the required header keywords
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	degree = pi/180.0           ! one degree in rad	

	!----- constant size of all pixels in [sr] -----
	! dist is from Parameter.in (is already converted to cm)
	delpix=2.d0*Rout*AU/real(n)	! this is in cm
	spixdegree  = delpix/(distance*parsec)/degree



	! put some coordinate system 
	centerpix=n/2+1.0 ! FIXME: requires odd number of pixesl
	call ftpkys(unit,'ctype1','RA---SIN','',status)
	call ftpkys(unit,'cunit1','deg     ','',status)
	call ftpkyd(unit,'crval1',68.0,12,'dummy value',status)
	call ftpkyd(unit,'cdelt1',spixdegree,12,'',status)
	call ftpkyd(unit,'crota1',0.0,12,'',status)
	call ftpkyd(unit,'crpix1',centerpix,12,'center pixel',status)

	call ftpkys(unit,'ctype2','DEC--SIN','',status)
	call ftpkys(unit,'cunit2','deg     ','',status)
	call ftpkyd(unit,'crval2',24.0,12,'dummy value',status)
	call ftpkyd(unit,'cdelt2',spixdegree,12,'',status)
	call ftpkyd(unit,'crota2',0.0,12,'',status)
	call ftpkyd(unit,'crpix2',centerpix,12,'center pixel',status)
	call ftpkys(unit,'RADESYS','ICRS','',status)

	! put in some grid for the velocity 
	! reference is lmin

	! the the first point (which should be lmin as the reference
	refFrqHz=clight/(lmin/1.e4)
	delHz=-refFrqHz*vresolution/clight

    ! Frequency/Spectral coordinate
    ! the first channel (velo) point is the reference
	call ftpkys(unit,'CTYPE3','FREQ','[Hz]',status)
	call ftpkys(unit,'cunit3','Hz      ','',status)
	call ftpkyd(unit,'crval3',refFrqHz,12,'Reference FREQ',status)
	call ftpkyd(unit,'CDELT3',delHz,12,'d_Hz (channel width)',status)
	call ftpkyd(unit,'crpix3',1.0,12,'Reference channel',status)
	! this is for velocity spectral units.
	call ftpkys(unit,'SPECSYS','LSRK','Spectral reference frame',status)
	call ftpkyj(unit,'VELREF',257,'Radio velocity in CASA',status)
	call ftpkys(unit,'timesys','UTC     ','',status)
	call ftpkys(unit,'cellscal','CONSTANT','',status)
	call ftpkys(unit,'origin','FLiTs/ProDiMo model','',status)
	! FIXME: does not make so much sense here if we have not a single line spectrum 
	! but for now focus on a single line spectrum 
	call ftpkyd(unit,'RESTFREQ',refFrqHz,12,'',status) ! this means first spectral axis point would give v=0


C     write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)*naxes(3)
      call ftpprd(unit,group,fpixel,nelements,im,status)

C     close the file and free the unit number
      call ftclos(unit, status)
      call ftfiou(unit, status)


	return
	end



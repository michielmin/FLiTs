	subroutine writefitsfile(filename,im,nlam,n)
	use GlobalSetup, only: distance
	use Constants, only: parsec
	IMPLICIT NONE
	character*500 filename
	integer n,nlam
	real*8 im(n,n,nlam)

      integer status,unit,blocksize,bitpix,naxis,naxes(3)
      integer i,j,group,fpixel,nelements
      logical simple,extend,truefalse

	inquire(file=filename,exist=truefalse)
	if(truefalse) then
		call output("FITS file already exists, overwriting")
		open(unit=90,file=filename)
		close(unit=90,status='delete')
	endif

      status=0
C     Get an unused Logical Unit Number to use to create the FITS file
      call ftgiou(unit,status)
C     create the new empty FITS file
      blocksize=1
      call ftinit(unit,filename,blocksize,status)

C     initialize parameters about the FITS image (IMDIM x IMDIM 64-bit reals)
      simple=.true.
      bitpix=-64
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

C     write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)*naxes(3)
      call ftpprd(unit,group,fpixel,nelements,im*1e23/(distance*parsec)**2,status)

C     close the file and free the unit number
      call ftclos(unit, status)
      call ftfiou(unit, status)


	return
	end



c=========================================================================================
c The program FLiTs (Fast Line Tracer(s)) is a line raytracer especially designed for 
c protoplanetary disks. It reads in the temperature and density structure from somewhere
c else and computes the level populations of the lines either by LTE, or reads them in
c from ProDiMo.
c The program uses a single, unformatted input file which is the only argument to the 
c program needed. All arguments can be changed from the command line as well by using the
c -s keyword=value syntax.
c So in principle the program is run like:
c FLiTs input.dat -s keyword=value -s keyword2=value2
c That's it, have fun!
c=========================================================================================
	program FLiTs
	IMPLICIT NONE

	call output("=================================================================="
	call output("Let's get the show on the road!!")
	call output("=================================================================="

c initialization
	call Initialize()
c read in density, temperature, opacity, and local radiation field
	call ReadStructure()
c read in line data: which lines, Einstein coefficients of all lines
	call ReadLineData()
c setup/readin level populations
	call SetupLevels()
c prepare the remaining things in the structure for the raytracing
	call PrepareStructure()
	
c the action starts!!
c first setup the paths
	call SetupPaths()
c now do the continuum raytracing at the reduced spectral resolution
	call RaytraceContinuum()


c and now the real interesting part!!
c do the line raytracing
	call RaytraceLines()

c finally write all the output files
	call OutputFiles()
	
c well that's it. we seem to be done!
c have a good day
	call output("=================================================================="
	call output("Success!!")
	call output("Have a nice day!")
	call output("=================================================================="

	end
	

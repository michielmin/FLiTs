	subroutine VersionDateTime(string)
	character*500 version,string
	parameter(version=trim(__DATE__)//' '//trim(__TIME__))

	string=version
	
	return
	end
	
	
	character*500 function VersionGIT()
#include "gitversion.h"

	VersionGIT=gitversion
	
	return
	end
	

	subroutine VersionDateTime(string)
  character(len=500) :: version, string
  parameter(version=trim(__DATE__)//' '//trim(__TIME__))

  string = version

  return
end subroutine VersionDateTime


character(len=500) function VersionGIT()
#include "gitversion.h"

  VersionGIT = gitversion

  return
end function VersionGIT
	

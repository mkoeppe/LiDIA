//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/lidia_file.h"
#include	<sys/types.h>
#include	<sys/stat.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool
lidia_file::file_exists (const char * fn)
{
	struct stat	statbuf;

	if (stat(fn, &statbuf) != 0) {
		return false;
	}

	return true;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

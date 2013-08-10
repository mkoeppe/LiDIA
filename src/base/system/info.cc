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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/LiDIA.h"
#include	<cstring>
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



int qo_info = 0;
bool qo_special = false;



void lidia_internal_xinfo(const char *v, lidia_size_t a, lidia_size_t b)
{
	// GET ENV Variable
	static char *LIDIA_INFO_MODE = NULL;

	if (LIDIA_INFO_MODE == NULL) {
		LIDIA_INFO_MODE = getenv("LIDIA_INFO_MODE");
		if (LIDIA_INFO_MODE == NULL) {
			LIDIA_INFO_MODE = new char[256];
			strcpy(LIDIA_INFO_MODE, "NONE");
		}
	}
	if (strcmp(LIDIA_INFO_MODE, "MODE_1") == 0) {
		char buffer[256];  // This should really be replaced by a std::string and sstream!
		sprintf(buffer, "%s:                    (%d/%d)                 ", v, a, b);
		std::cout << buffer << std::flush;
		lidia_size_t len = strlen(buffer);
		for (lidia_size_t j = 0; j < len; j++)
			std::cout << "\b" << std::flush;
	}
	else
		if (strcmp(LIDIA_INFO_MODE, "MODE_2") == 0) {
			char buffer[256];  // This should really be replaced by a std::string and sstream!
			sprintf(buffer, "%s:                    (%d/%d)              ", v, a, b);
			std::cout << buffer << std::endl;
		}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

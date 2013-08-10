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
#include	"LiDIA/path.h"
#include	"LiDIA/LiDIA.h"

#include	<cstring>
#include	<cstdlib>
#include	<fstream>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifndef LIDIA_DEBUGLIST_NAME
# define LIDIA_DEBUGLIST_NAME	"/dev/null"
#endif



bool check_debug_list(const char *f, const char *m)
{

	static char (*class_list)[256] = NULL;
	static char (*function_list)[256] = NULL;

	// GET ENV Variable
	static lidia_size_t ANZAHL = 0;
	static char *LIDIA_DEBUG_NAME = NULL;

	if (LIDIA_DEBUG_NAME == NULL) {
		LIDIA_DEBUG_NAME = getenv("LIDIA_DEBUG_NAME");
		if (LIDIA_DEBUG_NAME == NULL) {
			LIDIA_DEBUG_NAME = new char[256];
			strcpy(LIDIA_DEBUG_NAME, LIDIA_DEBUGLIST_NAME);
		}
	}
	if (class_list == NULL) {
		std::ifstream in(LIDIA_DEBUG_NAME);
		if (in.fail()) {
			return true;
		}

		in >> ANZAHL;

		class_list = new char[ANZAHL][256];
		function_list = new char[ANZAHL][256];

		for (register lidia_size_t i = 0; i < ANZAHL; i++) {
			in >> class_list[i];
			in >> function_list[i];

			//std::cout << class_list[i] << " - " << function_list[i] << std::endl;
		}
		in.close();
	}

	for (register lidia_size_t i = 0; i < ANZAHL; i++)
		if (strstr(f, class_list[i]) != NULL && (strstr(m, function_list[i]) != NULL || !strcmp("*", function_list[i])))
			return false;
	return true;
}



void default_debug_handler(const char *f, const char *m)
{
	std::cout << "\n debug_handler";
	std::cout << "::" << f;
	std::cout << "::" << m;
	std::cout << "\n";
	std::cout.flush();
}



debug_handler_ptr lidia_debug_handler = default_debug_handler;



debug_handler_ptr set_debug_handler(debug_handler_ptr new_handler)
{
	debug_handler_ptr old_handler = lidia_debug_handler;
	lidia_debug_handler = new_handler;
	return old_handler;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

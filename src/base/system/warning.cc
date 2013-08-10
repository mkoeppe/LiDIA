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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/LiDIA.h"
#include	<iostream>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void default_warning_handler(char *f, char *m)
{
	std::cout << "\n warning_handler";
	std::cout << "::" << f;
	std::cout << "::" << m;
	std::cout << "\n";
	std::cout.flush();
}



warning_handler_ptr lidia_warning_handler = default_warning_handler;



warning_handler_ptr set_warning_handler(warning_handler_ptr new_handler)
{
	warning_handler_ptr old_handler = lidia_warning_handler;
	lidia_warning_handler = new_handler;
	return old_handler;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

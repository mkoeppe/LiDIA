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
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void default_memory_handler(char *f, char *m)
{
	std::cout << "\n memory_handler";
	std::cout << "::" << f;
	std::cout << "::" << m;
	std::cout << "memory exhausted\n";
	std::cout.flush();
	abort();
}



memory_handler_ptr lidia_memory_handler = default_memory_handler;



memory_handler_ptr set_memory_handler(memory_handler_ptr new_handler)
{
	memory_handler_ptr old_handler = lidia_memory_handler;
	lidia_memory_handler = new_handler;
	return old_handler;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

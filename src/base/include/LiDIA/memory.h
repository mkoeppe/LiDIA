// -*- C++ -*-
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
//	$Id: memory.h,v 2.3 2002/06/24 09:41:30 lidiaadm Exp $
//
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_MEMORY_H_GUARD_
#define LIDIA_MEMORY_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



typedef void (*memory_handler_ptr)(char *, char *);

extern memory_handler_ptr lidia_memory_handler;
void default_memory_handler(char *, char *);
memory_handler_ptr set_memory_handler(memory_handler_ptr);

#ifdef LIDIA_MEMORY
#define memory_handler(vec, name, message) if ((vec) == 0) \
				     lidia_memory_handler(name, message)
#else
#define memory_handler(vec, name, message) { }
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MEMORY_H_GUARD_

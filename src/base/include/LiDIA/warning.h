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
//	$Id: warning.h,v 2.3 2002/06/24 09:41:34 lidiaadm Exp $
//
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_WARNING_H_GUARD_
#define LIDIA_WARNING_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



typedef void (*warning_handler_ptr)(char *, char *);

extern warning_handler_ptr lidia_warning_handler;
void default_warning_handler(char *, char *);
warning_handler_ptr set_warning_handler(warning_handler_ptr);

#ifdef LIDIA_WARNINGS
#define warning_handler(name, message) lidia_warning_handler(name, message)
#else
#define warning_handler(name, message) { }
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_WARNING_H_GUARD_

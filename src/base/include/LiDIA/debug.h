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
//	$Id: debug.h,v 2.3 2002/06/24 09:41:25 lidiaadm Exp $
//
//	Author	: Thomas Papanikolaou (TP), Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_DEBUG_H_GUARD_
#define LIDIA_DEBUG_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



typedef void (*debug_handler_ptr)(const char *, const char *);

extern debug_handler_ptr lidia_debug_handler;
void default_debug_handler(const char *, const char *);
debug_handler_ptr set_debug_handler(debug_handler_ptr);

bool check_debug_list(const char *, const char *);


#ifdef LIDIA_DEBUG
#define debug_handler(name, message) {if (check_debug_list(name, message))\
  lidia_debug_handler(name, message); }
#define debug_handler_l(name, message, level) {if (LIDIA_DEBUG <= level)\
  {if (check_debug_list(name, message))\
    {lidia_debug_handler(name, message); }}}

#define debug_handler_c(name, message, level, code) {if (LIDIA_DEBUG <= level)\
  if (check_debug_list(name, message))\
    {lidia_debug_handler(name, message); code; }}
#else
#define debug_handler(name, message) { }
#define debug_handler_l(name, message, level) { }
#define debug_handler_c(name, message, level, code) { }
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_DEBUG_H_GUARD_

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
//	$Id$
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_INFO_H_GUARD_
#define LIDIA_INFO_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



extern int qo_info; // whether to print timings or not (quadratic order)
extern bool qo_special;

void lidia_internal_xinfo(const char *, lidia_size_t, lidia_size_t);

#ifdef LIDIA_INFO
# define lidia_info_handler(a) {a}
# define lidia_xinfo_handler(a, b, c) lidia_internal_xinfo(a, b, c);

# define lidia_qo_info_handler(a) {if (!qo_special && (qo_info > 1)) {a}}
# define lidia_qo_xinfo_handler(a, b, c) {if (!qo_special && (qo_info > 1)) {lidia_internal_xinfo(a, b, c); }}
#else

# define lidia_info_handler(a) {}
# define lidia_xinfo_handler(a, b, c) {}

# define lidia_qo_info_handler(a) {}
# define lidia_qo_xinfo_handler(a, b, c) {}
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_INFO_H_GUARD_

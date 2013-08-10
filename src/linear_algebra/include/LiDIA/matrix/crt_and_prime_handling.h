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


#ifndef LIDIA_CRT_AND_PRIME_HANDLING_H_GUARD_
#define LIDIA_CRT_AND_PRIME_HANDLING_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// CRT
//

bigint chinrest(const bigint *values, const bigint *prim);

void chinrest(bigint & RES, const bigint * values, const bigint * prim);



//
// prime handling
//

bigint * get_primes(const bigint & C, const bigint & m, const bool SW_COPY = false);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_CRT_AND_PRIME_H_GUARD_ANDLING_INL

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
//	Author	: Damian Weber (DW)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_DLP_H_GUARD_
#define LIDIA_DLP_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
# include	"LiDIA/rational_factorization.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



bigint search_generator (const bigint &p, const bigint &g0);
bigint search_generator (const bigint &p);

#ifndef HEADBANGER
bigint search_generator (const bigint &p, const bigint &g0,
			 rational_factorization &f);

bigint pohlig_hellman_shanks (const bigint &g, const bigint &a,
			      const bigint &p, const bigint &q);
bigint discrete_log (const bigint &a, const bigint &b,
		     const bigint &p,
		     const rational_factorization &f);
#endif

extern bigint dl (const bigint &a, const bigint &b, const bigint &p,
		  int verbose = 0);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_DLP_H_GUARD_

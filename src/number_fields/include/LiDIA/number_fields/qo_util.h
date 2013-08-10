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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_QO_UTIL_H_GUARD_
#define LIDIA_QO_UTIL_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_PAIR_H_GUARD_
# include	"LiDIA/pair.h"
#endif
#ifndef LIDIA_MATRIX_H_GUARD_
# include	"LiDIA/matrix.h"
#endif
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// pretty runtime output function
void MyTime(long t);

// utility function for group structure computations
void decode_vector(matrix< bigint > & Bmat, const bigint & Bjj,
                   long r, long q, base_vector< long > & Rvec,
                   base_vector< long > & Qvec, long nR, long nQ);


bigint int_key(const int & G);
bigint long_key(const long & G);
bigint bigint_key(const bigint & G);
bigint pair_bigint_key(const pair< bigint, bigint > & G);
bigint bigrational_key(const bigrational & G);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QO_UTIL_H_GUARD_

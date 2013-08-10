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
//	Author	: Nigel Smart (NS)
//                Adaption of John Cremona's code
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_EC_ARITH_H_GUARD_
#define LIDIA_EC_ARITH_H_GUARD_


#ifndef LIDIA_BIGCOMPLEX_H_GUARD_
# include	"LiDIA/bigcomplex.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_NMBRTHRY_FUNCTIONS_H_GUARD_
# include	"LiDIA/nmbrthry_functions.h"
#endif
#ifndef LIDIA_MODULAR_OPERATIONS_INL_GUARD_
# include	"LiDIA/modular_operations.inl"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// Complex agm algorithm (should perhaps have real one as well)
bigcomplex cagm(const bigcomplex& a, const bigcomplex& b);


// Routines more at home in elliptic curves
// Why would anyone else want to use them ?

base_vector< bigint > int_roots_cubic(const bigint& a, const bigint& b,
				      const bigint& c);
bigcomplex* solve_cubic(const bigcomplex& c1, const bigcomplex& c2,
                        const bigcomplex& c3);
bigcomplex* solve_quartic(const bigcomplex& a, const bigcomplex& b, const bigcomplex& c, const bigcomplex& d);
bigcomplex* solve_real_quartic(const bigfloat& a, const bigfloat& b, const bigfloat& c, const bigfloat& d, const bigfloat& e);
base_vector< bigint > square_divs(const bigint& N);
bigcomplex normalize(bigcomplex& w1, bigcomplex& w2);
void getc4c6(const bigcomplex& w1, const bigcomplex& w2,
             bigcomplex& c4, bigcomplex &c6);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_EC_ARITH_H_GUARD_

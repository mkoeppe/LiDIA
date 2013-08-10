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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SINGLE_FACTOR_CC_GUARD_
#define LIDIA_SINGLE_FACTOR_CC_GUARD_



#ifndef LIDIA_SINGLE_FACTOR_H_GUARD_
# include	"LiDIA/single_factor.h"
#endif
#ifndef LIDIA_FACTORIZATION_H_GUARD_
# include	"LiDIA/factorization.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//*********************************************************************
//				class single_factor< T >
//*********************************************************************

template< class T >
lidia_size_t
ord_divide (const single_factor< T > &a, single_factor< T > &b)
{
	if (a.is_one())
		lidia_error_handler("single_factor< T >", "ord_divide::1st argument mustn't be 1");

	single_factor< T > tmp;
	lidia_size_t e = 0;
	gcd(tmp, a, b);
	while (tmp == a) {
		e++;
		divide(b, b, a);
		gcd(tmp, a, b);
	}
	return e;
}



template< class T >
factorization< T >
single_factor< T >::factor () const
	//standard factorization algorithm,
	//called by factorization < T >::factor_all_components(), here: dummy version
{
	lidia_error_handler("single_factor< T >", "factor(void)::not implemented");
	factorization< T > F;
	F.assign(*this);
	return F;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SINGLE_FACTOR_CC_GUARD_

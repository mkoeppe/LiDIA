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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigrational.h"
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/bigcomplex.h"
#include	"LiDIA/polynomial.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



polynomial< bigint >::operator base_polynomial< bigcomplex > () const
{
	base_polynomial< bigcomplex > x;
	x.set_degree(deg);
	for (lidia_size_t i = 0; i <= deg; i++)
		x[i] = bigcomplex(bigfloat(coeff[i]));
	return x;
}



polynomial< bigrational >::operator base_polynomial< bigcomplex > () const
{
	base_polynomial< bigcomplex > x;
	x.set_degree(deg);
	for (lidia_size_t i = 0; i <= deg; i++)
		x[i] = bigcomplex(bigfloat(coeff[i].numerator(),
					   coeff[i].denominator()));
	return x;
}



polynomial< bigfloat >::operator base_polynomial< bigcomplex > () const
{
	base_polynomial< bigcomplex > x;
	x.set_degree(deg);
	for (lidia_size_t i = 0; i <= deg; i++)
		x[i] = bigcomplex(coeff[i]);
	return x;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

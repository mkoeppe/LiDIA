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
#include	"LiDIA/bigfloat.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void
besselj (bigfloat & y, int n, const bigfloat & x)
{
	bigfloat e(1UL), f(1UL), s;
	long i, m = 0, ex, t = bigfloat::binary_precision, z;

	//
	// use sign identities to make n positiv
	//

	if (n < 0) {
		n = -n;
		if (n % 2)
			m = 1;
	}

	//
	// calculate n!
	//

	for (i = 2; i <= n; i++)
		multiply(f, f, i);
	z = t;

	//
	// increase the precision if needed
	//

	ex = x.e + x.prec;
	if (ex > z)
		bigfloat::binary_precision = z = ex;

	//
	// calculate the series
	//

	power(y, x, n);
	y.e -= n;
	divide(y, y, f);
	ex = z + y.e;
	if (ex < 0)
		ex = -ex;
	ex += z;
	bigfloat::binary_precision = ex + bits_per_base_digit
		- (ex % bits_per_base_digit);
	f.assign(y);
	y.assign_one();
	e.assign_one();
	square(s, x);
	s.negate();
	s.e -= 2;

	i = 1;
	for (;;) {
		multiply(e, e, s);
		divide(e, e, i);
		divide(e, e, n + i);
		if (e.is_approx_zero())
			break;
		add(y, y, e);
		i++;
	}

	//
	// restore the precision; final multiplication; restore
	// the sign and delete local variables
	//

	bigfloat::binary_precision = t;
	multiply(y, y, f);
	if (m)
		y.negate();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

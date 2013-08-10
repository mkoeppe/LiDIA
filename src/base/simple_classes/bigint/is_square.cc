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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool
bigint::is_square (bigint & r) const
{
	// Modified for LiDIA from the PARI/GP function
	//
	// int carrecomplet(GEN x, GEN *pt) of arith1.c
	//
	// The following static arrays store the squares
	// mod 64, 63, 65 and 11. For example 0 and 1 are
	// squares mod 63 but 2 is not.
	static int squares_mod64[] =
	{1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
	 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
	 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
	 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0};

	static int squares_mod63[] =
	{1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
	 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0,
	 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
	 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};

	static int squares_mod65[] =
	{1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0,
	 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0,
	 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
	 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0,
	 1};

	static int squares_mod11[] =
	{1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0};

	static bigint r2;
	unsigned long t = 0;

	// Special cases
	switch (sign ()) {
	case -1:
		// A negative x cannot be a square
		return false;
	case 0:
		// The square root of zero is easy to calculate;-)
		r.assign_zero ();
		return true;
	case 1:
		// Here x is positive so try the first (and cheepest)
		// square test. Notice that x % 64 == x & 63
		t = least_significant_digit ();
		if (!squares_mod64[t & 63])
			return false;
	}

	// If x can be represented in a machine word then we can test
	// more efficiently.
	if (length () == 1) {
		int r = squares_mod63[t % 63] + squares_mod65[t % 65] +
			squares_mod11[t % 11];
		if (r != 3)
			return false;
	}
	else {
		if (!squares_mod63[remainder (*this, 63)])
			return false;
		if (!squares_mod65[remainder (*this, 65)])
			return false;
		if (!squares_mod11[remainder (*this, 11)])
			return false;
	}
	// We have lost: x is a square mod 64, 63, 65 and 11 so have to
	// take the sqrt. We calculate r = floor(sqrt(x)) and check
	// r^2 == x
	sqrt (r, *this);
	square (r2, r);
	return (r2.compare (*this) == 0);
}



bool
bigint::is_square () const
{
	bigint r;

	return is_square(r);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

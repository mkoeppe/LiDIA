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
//	Author	: Markus Maurer (MM), Stefan Neis (SN)
//		  Thomas Papanikolaou (TP), Tobias Hahn (TH)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_POWER_FUNCTIONS_H_GUARD_
#define LIDIA_POWER_FUNCTIONS_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif



//
// To use these template functions, the type T must
// provide the following functions:
//
//	T.is_one();
//	T.assign_one ();
//
//	operator = (const T &);
//
//      square  (T&, const T&, const T&);
//      multiply(T&, const T&, const T&);
//





#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// Function: power_left_to_right
//
// Method: c = a^b by repeated squaring,
//         traversing the bits from left to right
//
// Special cases:
//
//   b < 0  ->lidia_error_handler is called
//


template< class T > void
lidia_power_left_to_right (T & c, const T & a, unsigned long b);

#ifndef HAVE_HP_UX_CC
template< class T > void
lidia_power_left_to_right (T & c, const T & a, bigint b);

#else
template< class T > inline void
lidia_power_left_to_right (T & c, const T & a, bigint b)
{
	if (b.is_negative())
		lidia_error_handler("lidia_power_left_to_right(T, const T&, bigint)",
				    "Negative exponent.");

	else if (b.is_zero() || a.is_one())
		c.assign_one();

	else if (b.is_one())
		c = a;

	else {
		T multiplier;

		multiplier = a;
		c = a;

		// Note: The bitlength is at least 2, since b == 1
		// is handled above. The first bit of the bits
		// [0, ..., bit_length()-1] is used up by
		// the assignment c = a, so the next bit to handle
		// will be the one with index bit_length()-2.

		unsigned int i;

		for (i = b.bit_length()-2; i > 0; i--) {
			square(c, c);
			if (b.bit(i))
				multiply(c, c, multiplier);
		}

		// i == 0

		square(c, c);
		if (b.is_odd())
			multiply(c, c, multiplier);
	}
}
#endif




//
// The following function implements the left-to-right exponentiation too.
// The only difference is, that this function obtains an additional
// parameter, which is a pointer to a multiplication function that is used
// in the inner loop of the algorithm.
//


//
// Function: power_left_to_right
//
// Method: c = a^b by repeated squaring,
//         traversing the bits from left to right
//
// Special cases:
//
//   b < 0  ->lidia_error_handler is called
//


template< class T > void
lidia_power_left_to_right (T & c,
			   const T & a,
			   unsigned long b,
			   void (*fast_mul) (T&, const T&, const T&));



template< class T > void
lidia_power_left_to_right (T & c,
			   const T & a,
			   bigint b,
			   void (*fast_mul) (T&, const T&, const T&));


//
//
//
//
// right-to-left algorithms
//
//
//
//


//
// Function: lidia_power_right_to_left
//
// Method: c = a^b by repeated squaring,
//         traversing the bits from right to left
//
// Special cases:
//
//   b < 0  ->lidia_error_handler is called
//
// Note that this is slower than lidia_power_left_to_right,
// if the results are not bounded in size (e.g. for bigints),
// whereas it might be slightly faster, if all
// powers are small (e.g. for bigmods).
//

template< class T > void
lidia_power_right_to_left (T & c, const T & a, unsigned long b);

template< class T > void
lidia_power_right_to_left (T & c, const T & a, bigint b);

//
// Function: lidia_power_sliding_window
//
// Method: c = a^b by repeated squaring and multipliciation
//         with precomputed powers of a. If the window size k
//         is one, this algorithm is equivalent to
//         lidia_power_left_to_right.
//
// Input:  a - the base
//         b - the exponent
//         k - the window size
//
// Output: c = a^b
//
// Special cases:
//
//   b < 0  ->lidia_error_handler is called
//
// Note that this algorithm precomputes 2^(k-1) powers of a
//

template< class T > void
lidia_power_sliding_window  (T & c, const T & a, const bigint &b,
			     unsigned long k);


//
//
//
//
// interface versions
//
//
//
//


template< class T >
inline void
lidia_power (T & c, const T & a, const bigint & b)
{
	lidia_power_left_to_right(c, a, b);
}



template< class T >
inline void
lidia_power (T & c, const T & a, unsigned long b)
{
	lidia_power_left_to_right(c, a, b);
}



template< class T >
inline void
lidia_power (T & c,
	     const T & a,
	     const bigint & b,
	     void (*fast_mul) (T&, const T&, const T&))
{
	lidia_power_left_to_right(c, a, b, fast_mul);
}



template< class T >
inline void
lidia_power (T & c,
	     const T & a,
	     unsigned long b,
	     void (*fast_mul) (T&, const T&, const T&))
{
	lidia_power_left_to_right(c, a, b, fast_mul);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/power_functions.cc"
#endif



#endif	// LIDIA_POWER_FUNCTIONS_H_GUARD_

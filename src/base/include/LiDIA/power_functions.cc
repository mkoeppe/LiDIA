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


#ifndef LIDIA_POWER_FUNCTIONS_CC_GUARD_
#define LIDIA_POWER_FUNCTIONS_CC_GUARD_



#ifndef LIDIA_POWER_FUNCTIONS_H_GUARD_
# include	"LiDIA/power_functions.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
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
lidia_power_left_to_right (T & c, const T & a, unsigned long b)
{
	if (b == 0 || a.is_one())
		c.assign_one();

	else if (b == 1)
		c = a;

	else {
		unsigned long exponent;
		T multiplier;

		exponent = b;

		multiplier = a;
		c = a;

		// Note: The bitlength is at least 2, since b == 1
		// is handled above. The first bit of the bits
		// [0, ..., bit_length()-1] is used up by
		// the assignment c = a, so the next bit to handle
		// will be the one with index bit_length()-2.

		lidia_size_t length = 0;
		while (b != 0) {
			length++;
			b >>= 1;
		}
		length -= 2;

		// shift_left(i, 1, length);
		unsigned long i = 1UL << length;

		while (i != 0) {
			square(c, c);
			if (exponent & i)
				multiply(c, c, multiplier);
			i >>= 1;
		}
	}
}



#ifndef HAVE_HP_UX_CC
template< class T > void
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
			   void (*fast_mul) (T&, const T&, const T&))
{
	if ((b == 0) || a.is_one())
		c.assign_one();

	else if (b == 1)
		c = a;

	else {
		unsigned long exponent;
		T multiplier;

		exponent = b;

		multiplier = a;
		c = a;

		// Note: The bitlength is at least 2, since b == 1
		// is handled above. The first bit of the bits
		// [0, ..., bit_length()-1] is used up by
		// the assignment c = a, so the next bit to handle
		// will be the one with index bit_length()-2.

		lidia_size_t length = 0;
		while (b != 0) {
			length++;
			b >>= 1;
		}
		length -= 2;

		// shift_left(i, 1, length);
		unsigned long i = 1UL << length;

		while (i != 0) {
			square(c, c);
			if (exponent&i)
				fast_mul(c, c, multiplier);
			i >>= 1;
		}
	}
}



template< class T > void
lidia_power_left_to_right (T & c,
			   const T & a,
			   bigint b,
			   void (*fast_mul) (T&, const T&, const T&))
{
	if (b.is_negative())
		lidia_error_handler("lidia_power_left_to_right(T, const T&, E, void (*fast_mul) (T&, const T&, const T&))",
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
				fast_mul(c, c, multiplier);
		}

		// i == 0

		square(c, c);
		if (b.is_odd())
			fast_mul(c, c, multiplier);
	}
}



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
lidia_power_right_to_left (T & c, const T & a, unsigned long b)
{
	if ((b == 0) || a.is_one())
		c.assign_one();

	else {
		unsigned long exponent;
		T multiplier;

		exponent = b;
		multiplier = a;

		c.assign_one();

		while (exponent != 1) {
			if (exponent & 1)
				multiply(c, c, multiplier);
			square(multiplier, multiplier);
			exponent >>= 1;
		}

		multiply(c, c, multiplier);
	}
}



template< class T > void
lidia_power_right_to_left (T & c, const T & a, bigint b)
{
	if (b.is_negative())
		lidia_error_handler("lidia_power_right_to_left(T, const T&, E)",
				    "Negative exponent.");

	else if (b.is_zero() || a.is_one())
		c.assign_one();

	else {
		T multiplier;

		multiplier = a;
		c.assign_one();

		while (!b.is_one()) {
			if (b.is_odd())
				multiply(c, c, multiplier);
			square(multiplier, multiplier);
			shift_right (b, b, 1);
		}

		multiply(c, c, multiplier);
	}
}



//
// Calculates c = a^b by using a sliding window technique
// of window size k.
//
template< class T > void
lidia_power_sliding_window (T & c, const T & a, const bigint & b,
			    unsigned long k)
{
	if (b.is_negative())
		lidia_error_handler("lidia_power_sliding_window",
				    "Negative exponent.");

	else if (k< 1 || k > sizeof (long)*8)
		lidia_error_handler("lidia_power_sliding_window",
				    "k not in [1, sizeof (long)*8]");

	else if (b.is_zero() || a.is_one())
		c.assign_one();

	else {
		// Precompute a, a^2 and all odd powers of a up to a^(2^k - 1)

		// a2 = a^2
		T a2;
		square (a2, a);

		// Array to hold odd powers.
		T *odd = NULL;

		// Allocate memory for odd powers.
		// There are 2^(k-1) - 1 precomputations left.
		unsigned long no_powers = 1 << (k-1);
		odd = new T [no_powers];

		// Calculate first entry.
		odd [0] = a;

		// INVARIANT: For all j in [0, i-1]. odd[j] == a^(2j + 1)
		for (unsigned long ii = 1; ii < no_powers; ii++) {
			multiply (odd [ii], a2, odd [ii-1]);
		}

		// Initialize loop.
		c.assign_one ();
		long bitlength = b.bit_length ();
		long i = bitlength - 1;

		// INVARIANT: c = a^(b >> (i + 1))
		//            i < 0 =  > i == -1
		while (i >= 0) {
			if (b.bit (i) == 0) {
				// In this case just square the result and move one
				// position to the right.
				square (c, c);
				i--;
			} else {
				// Get maximum window that contains an odd number.
				// Note that we might not have k bits left.
				long l = (static_cast<unsigned long>(i) >= k ? i - k + 1 : 0);
				// INVARIANT: For each j in [max (i-k+1, 0), l-1].  b.bit (j) == 0
				//            l <= i
				while (b.bit (l) == 0)
					l++;

				// Get all bits in [l, i]
				bigint mask = 1;
				mask = mask << (i - l + 1);
				mask --;
				mask = mask << l;
				bigint bigIndex = (b & mask) >> l;
				// This is one point where k > sizeof (long)*8 would be a problem.
				long index;
				bigIndex.longify (index);
				// So far, index is the exponent we need. To look it up we need
				// to convert this into an index:
				index = (index + 1)/2 - 1;

				// "Left shift" the result.
				long j = 0;
				// INVARIANT: c = c0^(2^j)
				while (j < i-l+1) {
					j++;
					square (c, c);
				}

				// Multiply with a^index
				multiply (c, c, odd[index]);
				i = l - 1;
			}
		}
		delete [] odd;
	}
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_POWER_FUNCTIONS_CC_GUARD_

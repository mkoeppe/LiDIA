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
//	Author	: Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/quadratic_form.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



base_vector< quadratic_form >
compute_class_group(const bigint & Delta)
{
	//
	// We check whether Delta is a negative integer congruent
	// 0 or 1 modulo 4.
	//
	
	bigint abs_Delta;
	abs_Delta.assign( Delta );
	abs_Delta.negate();

	if( Delta.is_positive() ||
	    ( abs_Delta.bit(0) == 0 && abs_Delta.bit(1) == 1) ||
	    ( abs_Delta.bit(0) == 1 && abs_Delta.bit(1) == 0))
		lidia_error_handler("compute_class_group(bigint & Delta)",
				    "Delta has to be a negative integer congruent 0 or 1 modulo 4.");

	//
	// A temporary variable of type bigint.
	//
	bigint tmp(Delta);

	//
	// Integers 0 and 1 as constant objects of type "bigint".
	//
	const bigint bigint_zero(0);
	const bigint bigint_one(1);

	//
	// The quadratic forms of discriminant Delta are stored
	// in an instance of the class base_vector. The current
	// length of the instance is stored in the variable "length".
	//
	base_vector< quadratic_form > reduced_forms(0, EXPAND);
	lidia_size_t length = 0;

	//
	// sqrt(|Delta|/3) is an upper bound of a and therefore for b either
	//
	bigint upper_bound_of_b;
	tmp.abs();
	divide(tmp, tmp, bigint(3));
	sqrt(upper_bound_of_b, tmp);
	// We can use the function void sqrt( bigint &, bigint & ); for
	// a proof see ...

	//
	// A reduced form is represented by ( a, b, candidate_of_c ).
	//
	bigint a(1);
	bigint candidate_of_c(0);
	bigint b = abs_Delta.bit(0); // b is even if and only if Delta is even.

	//
	// Assign the main form to the first entry of "reduced_forms".
	//
	if (b == bigint_zero) {
		divide(candidate_of_c, Delta, bigint(4));
		candidate_of_c.negate();
		reduced_forms[0].assign(bigint_one, bigint_zero, candidate_of_c);
	}
	else {
		subtract(candidate_of_c, bigint_one, Delta);
		divide(candidate_of_c, candidate_of_c, bigint(4));
		reduced_forms[0].assign(bigint_one, bigint_one, candidate_of_c);
	}

	//
	// candidate_of_c_times_a <- ( b^2 - Delta )/4
	//
	bigint candidate_of_c_times_a(candidate_of_c);

	//
	// Increase a by 1: a <- a + 1.
	//
	inc(a);

	//
	// The loop on b: Test all possible values of b.
	//
	while (b <= upper_bound_of_b) {
		//
		// The loop on a: Test all possible values of a.
		//
		while (a*a <= candidate_of_c_times_a) {
			//
			// Test whether a divides candidate_of_c_times_a
			//
			remainder(tmp, candidate_of_c_times_a, a);
			if (tmp == bigint_zero) {
				//
				// If ( a, b, candidate_of_c ) is reduced
				// then gcd( a, b, candidate_of_c ) = 1 follows.
				//
				divide(candidate_of_c, candidate_of_c_times_a, a);
				if (gcd(gcd(a, b), candidate_of_c) == bigint_one) {
					//
					// Concatenate the forms to the array "reduced_forms".
					//
					length++;
					reduced_forms[length].assign(a, b, candidate_of_c);

					if (b != 0 && a != b && a != candidate_of_c) {
						length++;
						reduced_forms[length].assign(a, -b, candidate_of_c);
					}
				}
			}
			inc(a);
		}

		//
		// Adapt candidate_of_c_times_a.
		//
		candidate_of_c_times_a += b + 1;

		//
		// b <- b + 2,   a <- b.
		//
		inc(b);
		inc(b);
		a = b;
	}

	return(reduced_forms);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

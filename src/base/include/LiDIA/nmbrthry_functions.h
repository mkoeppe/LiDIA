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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_NMBRTHRY_FUNCTIONS_H_GUARD_
#define LIDIA_NMBRTHRY_FUNCTIONS_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
# include	"LiDIA/rational_factorization.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// file: divisors.cc

void h_divisors     (sort_vector< bigint > & div, rational_factorization & f, int flag = 0);
void h_all_divisors (sort_vector< bigint > & div, rational_factorization & f, int flag = 0);

sort_vector< bigint > divisors (rational_factorization & f);
sort_vector< bigint > divisors (const bigint & N);

sort_vector< bigint > all_divisors (rational_factorization & f);
sort_vector< bigint > all_divisors (const bigint & N);

sort_vector< bigint > square_free_divisors (rational_factorization & f);
sort_vector< bigint > square_free_divisors (const bigint & N);

sort_vector< bigint > all_square_free_divisors (rational_factorization & f);
sort_vector< bigint > all_square_free_divisors (const bigint & N);

sort_vector< bigint > square_divides_n_divisors (rational_factorization & f);
sort_vector< bigint > square_divides_n_divisors (const bigint & N);

sort_vector< bigint > all_square_divides_n_divisors (rational_factorization & f);
sort_vector< bigint > all_square_divides_n_divisors (const bigint & N);


bool is_divisor (const bigint & factor, const bigint & number);
int valuation(const bigint & factor, const bigint & number);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_NMBRTHRY_FUNCTIONS_H_GUARD_

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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_MV_FRACT_H_GUARD_
#define LIDIA_MV_FRACT_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_GF2N_H_GUARD_
# include	"LiDIA/gf2n.h"
#endif
#ifndef LIDIA_MV_POLY_H_GUARD_
# include	"LiDIA/mv_poly.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class mv_fract
{
	//
	// this class represents multivariate rational funktion
	// the numerator and denominator are of type mv_poly
	//
	//

	mv_poly numerator;
	mv_poly denominator;

public:

	//
	// constructors and destructor
	//

	mv_fract();
	mv_fract(const mv_poly & n);
	mv_fract(const mv_poly & n, const mv_poly & d);
	mv_fract(const mv_fract & p);
	~mv_fract();

	//
	// assignments
	//

	friend void swap (mv_fract & a, mv_fract & b);

	mv_fract& operator = (const mv_fract & p);
	mv_fract& operator = (const mv_poly & t);

	void assign_zero()
	{
		numerator.assign_zero();
		denominator.assign_one();
	}

	void assign_one()
	{
		numerator.assign_one();
		denominator.assign_one();
	}

	void assign(const mv_poly & n, const mv_poly & d)
	{
		numerator = n;
		denominator = d;
	}

	void assign(const mv_poly & n)
	{
		numerator = n;
		denominator.assign_one();
	}

	//
	// access
	//
	const mv_poly & get_numerator() const
	{
		return numerator;
	}
	const mv_poly & get_denominator() const
	{
		return denominator;
	}

	//
	// comparisons
	//

	friend bool operator == (const mv_fract & a, const mv_fract & b);

	// other comparisons

	bool is_zero() const;
	bool is_one() const;

	//
	// functions
	//


	// If the numerator of q contains linear term with variable X_k
	// the function solves for X_k. k is set to the variable index, q is set
	// to the result (note: it becomes polynomial !!), and true is returned.
	//
	// Otherwise the return value is false, and q will is not changed.

	friend bool solve_rat_one(mv_fract & q, lidia_size_t & k);


	// If the numerator of q contains linear term with variable X_k and the
	// the denominator does not have X_k, the function solves the numerator for
	// X_k. q is set to the result  (note: it becomes polynomial !!), and
	// true is returned.
	// Otherwise the return value is false, and q will is not changed.

	friend bool solve_rat_fixed(mv_fract & q, lidia_size_t k);


	// substitute all X_k in q with the rational function p (that does not
	// contain X_k itself), return result in q.

	friend void substitute(mv_fract & q, const mv_fract & p, lidia_size_t k);


	// variable X_k in t is substituted by bit 'bit'

	friend void evaluate_one(mv_fract & t, lidia_size_t k, unsigned char bit);


	// check whether either numerator or denominator of *this have variable X_k
	// If so, return true, otherwise false

	bool has_var_k(lidia_size_t k) const;

	// clean both numerator and denominator from common factor, constant
	// denominators are changed to 1.


	void clean();

	// split the numerator of t = t1 * X_k + t2, set t = t1/t2
	// if t does not have variable X_k, return error.

	friend void split_x_k(mv_fract & t, lidia_size_t k);


	//
	//  output
	//

	friend std::ostream & operator << (std::ostream & out, const mv_fract & a);
};



inline void swap (mv_fract & a, mv_fract & b)
{
	swap(a.numerator, b.numerator);
	swap(a.denominator, b.denominator);
}



inline bool operator != (const mv_fract & a, const mv_fract & b)
{
	return !(a == b);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MV_FRACT_H_GUARD_

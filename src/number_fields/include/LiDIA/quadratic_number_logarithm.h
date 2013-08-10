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


#ifndef LIDIA_QUADRATIC_NUMBER_LOGARITHM_H_GUARD_
#define LIDIA_QUADRATIC_NUMBER_LOGARITHM_H_GUARD_


#ifndef LIDIA_QUADRATIC_IDEAL_H_GUARD_
# include	"LiDIA/quadratic_ideal.h"
#endif
#ifndef LIDIA_XBIGFLOAT_H_GUARD_
# include	"LiDIA/xbigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// A quadratic number logarithm is a tupel L = (A, B, a, k, s, t), where A is a
// quadratic ideal of the quadratic order O, B is a reduced quadratic ideal
// of O, a an xbigfloat, k >= 4, and s = +-1, such that there exists a
// number alpha in the field of fractions of O with
//
//    1) A / alpha = B
//    2) |a - Ln alpha| < 2^{-k}, if t = 0,
//       |a - ln alpha| < 2^{-k}, if t = 1,
//    3) s = sign(alpha), (1, if alpha >= 0, and -1, otherwise).
//
// The tupel L represents alpha uniquely, i.e., if there is element gamma
// in the field of fractions of O satisfying 1)-3), then gamma = alpha.
//

class quadratic_number_logarithm
{
private:

	quadratic_ideal A;
	quadratic_ideal B;
	xbigfloat       a;
	long            k;
	int             s;
	int             t;

	//
	//  constructors / destructor
	//

public:
	quadratic_number_logarithm ();
	quadratic_number_logarithm (const quadratic_number_logarithm & L);
	~quadratic_number_logarithm ();

	//
	//  access
	//

	bigint get_discriminant () const;
	const quadratic_order & which_order() const;
	const quadratic_order & get_order() const;

	const quadratic_ideal & get_A () const;
	const quadratic_ideal & get_B () const;

	const xbigfloat & get_log_approximation () const;
	long get_log_accuracy () const;
	int get_log_type () const;

	int get_sign() const;


	//
	//  assignments
	//

	quadratic_number_logarithm & operator = (const quadratic_number_logarithm & L);

	void assign (const quadratic_number_logarithm & L);

	void assign (const quadratic_ideal & pA,
		     const quadratic_ideal & pB,
		     const xbigfloat       & pa,
		     long                    pk,
		     int                     ps,
		     int                     pt = 0);

	friend void swap (quadratic_number_logarithm & x,
			  quadratic_number_logarithm & y);

	//
	//  arithmetic
	//

	void invert();


	//
	//  input / output
	//

	friend std::istream &
	operator >> (std::istream & in, quadratic_number_logarithm &);

	friend std::ostream &
	operator << (std::ostream & out, const quadratic_number_logarithm &);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QUADRATIC_NUMBER_LOGARITHM_H_GUARD_

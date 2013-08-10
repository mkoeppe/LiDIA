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
//	Author	: Volker Mueller (VM), Markus Maurer (MM), Andrea Rau (AR)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_EC_DOMAIN_PARAMETERS_P1363_H_GUARD_
#define LIDIA_EC_DOMAIN_PARAMETERS_P1363_H_GUARD_



#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif
#ifndef LIDIA_POINT_H_GUARD_
# include	"LiDIA/point.h"
#endif
#ifndef LIDIA_GF_ELEMENT_H_GUARD_
# include	"LiDIA/gf_element.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_TIMER_H_GUARD_
# include	"LiDIA/timer.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



const int MUCH_INFO = 2;
const int LITTLE_INFO = 1;
const int NO_INFO = 0;



class EC_domain_parameters_P1363
{
	//
	// q   Size of the finite field.
	// a, b Coefficients of the elliptic curve E.
	// r   The prime divisor of #E.
	// k   The cofactor #E / r.
	// x   x-coordinate of point on E over GF(q) of order r.
	// y   x-coordinate of point on E over GF(q) of order r.
	//

private:
	static lidia_size_t defaultBitsize_;
	static lidia_size_t defaultPercentage_;

	gf_element	                a, b;
	bigint               		k, r, q;
	point< gf_element >	 	G;
	bool initialized; 		// true, iff generate_parameters
	// has been called.
public:
	static int GF2N;
	static int GFP;

	//
	// constructor / destructor
	//
	EC_domain_parameters_P1363();
	~EC_domain_parameters_P1363();

	//
	// Access
	//
	// Error, if initialized == false
	//
	const bigint & get_q () const;
	const gf_element & get_a () const;
	const gf_element & get_b () const;
	const bigint & get_k () const;
	const bigint & get_r () const;
	const point< gf_element > & get_G () const;

	lidia_size_t default_bitsize();
	lidia_size_t default_percentage();

	//
	// Assignments
	//
	void assign(const EC_domain_parameters_P1363 & I);
	EC_domain_parameters_P1363 & operator = (const EC_domain_parameters_P1363 & I);

	//
	// High level functions
	//
	void generate_parameters(int field, int bitsize_factor, int info);

	void generate_parameters(int field, int info);

	void generate_parameters(int field, //GF2N oder GFP
				 int bitsize_factor,
				 int percentage,
				 int info);

private:
	bool is_strong_curve(rational_factorization & rf_order,
			     const bigint & order,
			     const bigint & co_factor,
			     const bigint & p,
			     int info) const;

	void get_twist_coeff(gf_element & new_a, gf_element & new_b);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_EC_DOMAIN_PARAMETERS_P1363_H_GUARD_

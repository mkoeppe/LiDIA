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
//	Author	: Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


// This class holds points on elliptic curves as rational
// functions modulo some gf2n_poly_modulus. This is nice for
// baby-giantstep algorithms and for rational function versions
// of the eigenvalue search algorithm.


#ifndef LIDIA_WECO2_RAT_FUNCTION_H_GUARD_
#define LIDIA_WECO2_RAT_FUNCTION_H_GUARD_



#ifndef LIDIA_GF2N_RATIONAL_FUNCTION_H_GUARD_
# include	"LiDIA/gf2n_rational_function.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class weco2_rat_function
{
public:
	typedef gf2n ff_element;
	typedef gf2n_polynomial ff_pol;
	typedef gf2n_rational_function  ff_rat;
	typedef gf2n_poly_modulus  ff_pol_m;

private:

	static ff_element  *a6; // coefficient of curve
	static ff_pol_m  modulus; // global modulus for rational functions
	static ff_pol      curve; // X^3 + a6

	ff_rat   x; // x-coordinate as polynomial in X
	ff_rat   y0, y1; // y-coordinate is given as Y * y1 + y0
	bool zero; // Zero-Point Yes/No


 public:

	weco2_rat_function();
	weco2_rat_function(const ff_rat & xx, const ff_rat & yy0);
	weco2_rat_function(const weco2_rat_function &);
	weco2_rat_function(const ff_rat & xx, const ff_rat & yy0, const ff_rat & yy1);
	~weco2_rat_function();

	static void initialize (const ff_element & a6, const ff_pol_m & pm);

	bool on_curve() const;

	void kill()
	{
		x.kill();
		y0.kill();
		y1.kill();
		zero = true;
	}

	// ---------- COMPARISONS ----------

	friend bool operator == (const weco2_rat_function & P,
				 const weco2_rat_function & Q);
	friend bool operator != (const weco2_rat_function & P,
				 const weco2_rat_function & Q)
	{
		return ! (P == Q);
	}

	bool is_zero() const
	{
		return zero;
	}


	// ---------- ARITHMETICAL OPERATIONS ----------

	friend void negate (weco2_rat_function & Q, const weco2_rat_function & P);
	// Q = -P

	friend void multiply_by_2 (weco2_rat_function & Q, const weco2_rat_function & P);
	// Q = 2*P

	friend void add_PQ (weco2_rat_function & R, const weco2_rat_function & Q,
			    const weco2_rat_function & P);
	// R = P + Q, P != +- Q, P != 0, Q != 0

	friend void add (weco2_rat_function & R, const weco2_rat_function & Q,
			 const weco2_rat_function & P);
	// R = P + Q

	friend void subtract (weco2_rat_function & R, const weco2_rat_function & Q,
			      const weco2_rat_function & P);
	// R = P - Q

	friend void multiply (weco2_rat_function & Q, lidia_size_t n,
			      const weco2_rat_function & P);
	// Q = n * P
	friend void multiply2 (weco2_rat_function & Q, lidia_size_t n,
			       const weco2_rat_function & P);
	// Q = n * P, alternative algorithm, comparison has to be done


	friend void multiply (weco2_rat_function & R1, weco2_rat_function & R2,
			      lidia_size_t n1, lidia_size_t n2,
			      const weco2_rat_function & P);
	// R1 = n1 * P, R2 = n2 * P, parallel multiplication


	// ---------- ASSIGNMENTS ----------

	weco2_rat_function & operator = (const weco2_rat_function& wep);
	void assign_xy ();
	void assign (const ff_rat &xx, const ff_rat &yy0);
	void assign (const ff_rat &xx, const ff_rat &yy0, const ff_rat & yy1);
	void assign (const weco2_rat_function & w);
	void assign_zero()
	{
		zero = true;
	}

	// ---------- ACCESS ----------

	const ff_rat & get_x() const
	{
		return x;
	}

	const ff_rat & get_y0() const
	{
		return y0;
	}

	const ff_rat & get_y1() const
	{
		return y1;
	}

	//------------- Output ---------

	friend std::ostream & operator << (std::ostream & o, const weco2_rat_function &);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_WECO2_RAT_FUCTION_H_GUARD_

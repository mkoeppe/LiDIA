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
//	Author	: Markus Maurer (MM), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


// This class holds points on elliptic curves as rational
// functions modulo some Fp_poly_modulus. This is nice for 
// baby-giantstep algorithms and for rational function versions
// of the eigenvalue search algorithm. 


#ifndef LIDIA_WEP_RAT_FUNCTION_H_GUARD_
#define LIDIA_WEP_RAT_FUNCTION_H_GUARD_



#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif
#ifndef LIDIA_FP_RATIONAL_FUNCTION_H_GUARD_
# include	"LiDIA/Fp_rational_function.h"
#endif
#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include	"LiDIA/Fp_polynomial.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// Here we use projective coordinates of the form (x : y : z) such
// that this point corresponds to the affine point (x / z^2, y/ z^3).
// This form is chosen since it allows fast additions formulas which are
// taken from the IEEE P1363 Annex A.


class wep_rat_function
{
public:		
	typedef bigmod ff_element; 
	typedef Fp_polynomial ff_pol; 
	typedef Fp_rational_function ff_rat; 
	typedef Fp_poly_modulus  ff_pol_m; 

private:

	static ff_element  a; // coefficients of curve
	static ff_element  b; 
	static ff_pol_m  modulus; // global modulus for rational functions

	ff_pol   x; // x-coordinate as polynomial in X
	ff_pol   y; // y-coordinate is given as Y * y
	ff_pol   z; // z-coordinate 

public:

	wep_rat_function(); 
	wep_rat_function(const ff_pol & xx, const ff_pol & yy); 
	~wep_rat_function(); 
 
	static void initialize (const ff_element & aa, const ff_element & bb, 
				const ff_pol_m & pm); 

	bool on_curve() const; 

private:
 
	friend void mult_by_curve(ff_pol & f, const ff_pol & ff); 

public:

	// ---------- COMPARISONS ----------

	friend bool operator == (const wep_rat_function & P, 
				 const wep_rat_function & Q); 
	friend bool operator != (const wep_rat_function & P, 
				 const wep_rat_function & Q)
	{
		return ! (P == Q);
	}

	bool is_zero() const
	{
		return z.is_zero();
	}


	// ---------- ARITHMETICAL OPERATIONS ----------

	friend void negate (wep_rat_function & Q, const wep_rat_function & P); 
	// Q = -P
 
	friend void multiply_by_2 (wep_rat_function & Q, const wep_rat_function & P); 
	// Q = 2*P
   
	friend void add_PQ (wep_rat_function & R, const wep_rat_function & Q, 
			    const wep_rat_function & P); 
	// R = P + Q, P != +- Q, P != 0, Q != 0
   
	friend void add (wep_rat_function & R, const wep_rat_function & Q, 
			 const wep_rat_function & P); 
	// R = P + Q
	
        friend void add_P_XY (wep_rat_function & R, const wep_rat_function &P);
	// R = P + (X,Y) 
 
	friend void subtract (wep_rat_function & R, const wep_rat_function & Q, 
			      const wep_rat_function & P); 
	// R = P - Q

	friend void multiply (wep_rat_function & Q, lidia_size_t n, 
			      const wep_rat_function & P); 
	// Q = n * P


	friend void multiply (wep_rat_function & R1, wep_rat_function & R2, 
			      lidia_size_t n1, lidia_size_t n2, 
			      const wep_rat_function & P); 
	// R1 = n1 * P, R2 = n2 * P, parallel multiplication


	// ---------- ASSIGNMENTS ----------
  
	wep_rat_function & operator = (const wep_rat_function& wep); 
	void assign_xy (); 
	void assign (const ff_pol &xx, const ff_pol &yy, const ff_pol & zz); 
	void assign (const ff_pol &xx, const ff_pol &yy); 
	void assign (const wep_rat_function & w); 
	void assign_zero()
	{
		x.assign_one();
		y.assign_one();
		z.assign_zero();
	}

	// ---------- ACCESS ----------
 
	ff_rat get_x() const
	{ 
		ff_pol h; 
		square(h, z, wep_rat_function::modulus); 
		return ff_rat(x, h); 
	}

	ff_rat get_y() const 
	{ 
		ff_pol h; 
		square(h, z, wep_rat_function::modulus); 
		multiply(h, h, z, wep_rat_function::modulus); 
		return ff_rat(y, h); 
	}

	friend std::ostream & operator << (std::ostream &, const wep_rat_function &); 

}; 



inline
void mult_by_curve(wep_rat_function::ff_pol & f,
		   const wep_rat_function::ff_pol & ff)
{
	wep_rat_function::ff_pol h;

	h.set_modulus(wep_rat_function::a.modulus());
	multiply_by_x_mod(h, f, ff);
	multiply_by_x_mod(h, h, ff);
	add(h, h, wep_rat_function::a.mantissa() * f);
	multiply_by_x_mod(h, h, ff);
	add(f, h, wep_rat_function::b.mantissa() * f);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_WEP_RAT_FUNCTION_H_GUARD_

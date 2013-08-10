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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_FP_POLY_MODULUS_H_GUARD_
#define LIDIA_FP_POLY_MODULUS_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include	"LiDIA/Fp_polynomial.h"
#endif
#ifndef LIDIA_FP_POLYNOMIAL_FFT_H_GUARD_
# include	"LiDIA/finite_fields/Fp_polynomial_fft.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#ifndef HEADBANGER


class Fp_ratfunc;
class Fp_poly_multiplier;



//******************************************************************
//
//        Modular Arithmetic with Pre-conditioning
//
//******************************************************************

// If you need to do a lot of arithmetic modulo a fixed f,
// build Fp_poly_modulus F for f.  This pre-computes information about f
// that speeds up the computation a great deal.
// f should be monic, and deg(f) > 0.

class Fp_poly_modulus
{

	friend void add      (Fp_ratfunc&, const Fp_ratfunc&, const Fp_ratfunc&);
	friend void subtract (Fp_ratfunc&, const Fp_ratfunc&, const Fp_ratfunc&);

	Fp_polynomial f; // modulus, a monic polynomial
	lidia_size_t n; // n = deg(f)
	bool use_FFT; // flag indicating whether FFT should be used.
	lidia_size_t crov; // crossover point

	lidia_size_t k; // least k s/t 2^k >= deg_f
	lidia_size_t l; // least l s/t 2^l >= 2deg_f-3
	fft_rep FRep; // 2^k point rep of f
	// H = rev((rev(f))^{-1} remainder x^{deg_f-1})
	fft_rep HRep; // 2^l point rep of H
	fft_data fd; // used for FRep and HRep


	void rem21(Fp_polynomial& x, const Fp_polynomial& a) const;
	// x = a % f
	// deg(a) <= 2(n-1), where n = f.degree()


public:

	const Fp_polynomial& modulus () const
	{
		debug_handler("Fp_poly_modulus", "modulus (void)");
		return f;
	}

	Fp_poly_modulus() : n(-1), use_FFT(false)
	{
		debug_handler("Fp_poly_modulus", "Fp_poly_modulus (void)");
	}

	Fp_poly_modulus(const Fp_poly_modulus & P)
	{
		debug_handler("Fp_poly_modulus", "Fp_poly_modulus (Fp_poly_modulus&)");
		build(P.modulus());
	}

	Fp_poly_modulus(const Fp_polynomial& ff) // MM
	{
		debug_handler("Fp_poly_modulus", "Fp_poly_modulus (Fp_polynomial&)");
		build(ff);
	}

	~Fp_poly_modulus()
	{
		debug_handler("Fp_poly_modulus", "~Fp_poly_modulus (void)");
	}

	void build (const Fp_polynomial& ff); // MM

	Fp_poly_modulus & operator = (const Fp_poly_modulus & P);


	// -------------------- friends --------------------


	friend void remainder (Fp_polynomial& x, const Fp_polynomial& a,
			       const Fp_poly_modulus& F);
    	// x = a % f, no restrictions on deg(a); makes repeated calls to rem21


	friend void multiply (Fp_polynomial& x, const Fp_polynomial& a,
			      const Fp_polynomial& b, const Fp_poly_modulus& F);
	// x = (a * b) % f
	// deg(a), deg(b) < f.degree()

	friend void multiply (Fp_polynomial& x, const Fp_polynomial& a,
			      const Fp_poly_multiplier& B, const Fp_poly_modulus& F);
	// x = (a * b) % f

	friend void square (Fp_polynomial& x, const Fp_polynomial& a,
			    const Fp_poly_modulus& F);
	// x = a^2 % f
	// deg(a) < f.degree()


	//fractions.[ch]
	friend void remainder (Fp_polynomial& x, fft_rep& R1,
			       const Fp_poly_modulus& F, modular_fft_rep& R2, Fp_polynomial& P1);
	friend void add_frac (Fp_polynomial& x, Fp_polynomial& y,
			      const Fp_polynomial& a, const Fp_polynomial& b,
			      const Fp_polynomial& c, const Fp_polynomial& d,
			      const Fp_poly_modulus& F);
	friend void subtract_frac (Fp_polynomial& x, Fp_polynomial& y,
				   const Fp_polynomial& a, const Fp_polynomial& b,
				   const Fp_polynomial& c, const Fp_polynomial& d,
				   const Fp_poly_modulus& F);
	friend bool eq_frac (const Fp_polynomial& a, const Fp_polynomial& b,
			     const Fp_polynomial& c, const Fp_polynomial& d,
			     const Fp_poly_modulus& F);
	//compose.[ch]
	friend void update_map (base_vector< bigint > & x, const base_vector< bigint > & a,
				const Fp_poly_multiplier& B, const Fp_poly_modulus& F);


	friend class Fp_poly_multiplier;
};

void remainder (Fp_polynomial& x, const Fp_polynomial& a,
		const Fp_poly_modulus& F);
    	// x = a % f, no restrictions on deg(a); makes repeated calls to rem21


void multiply (Fp_polynomial& x, const Fp_polynomial& a,
	       const Fp_polynomial& b, const Fp_poly_modulus& F);
	// x = (a * b) % f
	// deg(a), deg(b) < f.degree()

void multiply (Fp_polynomial& x, const Fp_polynomial& a,
	       const Fp_poly_multiplier& B, const Fp_poly_modulus& F);
	// x = (a * b) % f

void square (Fp_polynomial& x, const Fp_polynomial& a,
	     const Fp_poly_modulus& F);
	// x = a^2 % f
	// deg(a) < f.degree()


	//fractions.[ch]
void remainder (Fp_polynomial& x, fft_rep& R1,
		const Fp_poly_modulus& F, modular_fft_rep& R2, Fp_polynomial& P1);
void add_frac (Fp_polynomial& x, Fp_polynomial& y,
	       const Fp_polynomial& a, const Fp_polynomial& b,
	       const Fp_polynomial& c, const Fp_polynomial& d,
	       const Fp_poly_modulus& F);
void subtract_frac (Fp_polynomial& x, Fp_polynomial& y,
		    const Fp_polynomial& a, const Fp_polynomial& b,
		    const Fp_polynomial& c, const Fp_polynomial& d,
		    const Fp_poly_modulus& F);
bool eq_frac (const Fp_polynomial& a, const Fp_polynomial& b,
	      const Fp_polynomial& c, const Fp_polynomial& d,
	      const Fp_poly_modulus& F);
	//compose.[ch]
void update_map (base_vector< bigint > & x, const base_vector< bigint > & a,
		 const Fp_poly_multiplier& B, const Fp_poly_modulus& F);



typedef Fp_poly_modulus poly_modulus;


void power (Fp_polynomial& x, const Fp_polynomial& a, const bigint& e,
	    const Fp_poly_modulus& F);
// x = a^e % f

void power_x (Fp_polynomial& x, const bigint&  e, const Fp_poly_modulus& F);
// x = X^e % f

void power_x_plus_a (Fp_polynomial& x, const bigint& a, const bigint& e,
		     const Fp_poly_modulus& F);
// x = (X + a)^e % f


#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_FP_POLY_MODULUS_H_GUARD_

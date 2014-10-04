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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================

// Description  : poly_moduluss over gf2n, basic operations


#ifndef LIDIA_GF2N_POLY_MODULUS_H_GUARD_
#define LIDIA_GF2N_POLY_MODULUS_H_GUARD_


#ifndef LIDIA_GF2N_POLYNOMIAL_H_GUARD_
# include	"LiDIA/gf2n_polynomial.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class gf2n_poly_modulus
{
	friend class gf2n_polynomial;

private:

	gf2n_polynomial mod; // modulus
	gf2n_polynomial *tab; // table[i] = recip_{i+deg(m)}(mod)
	int size_tab; // tab[0, ..., size-1]
	int used_table_entries;

	static gf2n_polynomial help_pol;


public:

	gf2n_poly_modulus()
	{
		used_table_entries = 0;
		tab = NULL;
	}

	gf2n_poly_modulus (const gf2n_polynomial & p)
	{
		used_table_entries = 0;
		tab = NULL;
		build(p);
	}

	gf2n_poly_modulus(const gf2n_poly_modulus & p);

	~gf2n_poly_modulus()
	{
		if (tab != NULL)
			delete[] tab;

		// resize the stack for Karatzuba multiplications

		delete[] gf2n_polynomial::gf2n_p_tmp;
		gf2n_polynomial::gf2n_p_tmp = NULL;
		gf2n_polynomial::gf2n_p_top = gf2n_polynomial::gf2n_p_tmpsize = 0;
	}


	//*** basic functions  **********************************

	void build(const gf2n_polynomial &);

	const gf2n_polynomial & modulus () const
	{
		return mod;
	}

	gf2n_poly_modulus & operator = (const gf2n_poly_modulus & a)
	{
		assign(a);
		return *this;
	}

	void assign (const gf2n_poly_modulus & a);


	//********* functions which need reduction ***************


	friend void plain_rem(class gf2n_polynomial &, const class gf2n_polynomial &,
			      class gf2n_poly_modulus &);

	friend void remainder(class gf2n_polynomial & r,
			      const class gf2n_polynomial & a,
			      gf2n_poly_modulus & fp);

	friend void kara_rem (class gf2n_polynomial &, const class gf2n_polynomial &,
			      class gf2n_poly_modulus &);

	friend void multiply(gf2n_polynomial& x, const gf2n_polynomial& a,
			     const gf2n_polynomial& b, gf2n_poly_modulus & F);
	// x = (a * b) % f
	// deg(a), deg(b) < f.degree()


	friend void square(gf2n_polynomial& x, const gf2n_polynomial& a,
			   gf2n_poly_modulus & F);
	// x = a^2 % f
	// deg(a) < f.degree()


	friend void shift_left(gf2n_polynomial &erg, const gf2n_polynomial &a,
			       gf2n_poly_modulus &F, unsigned int d);
	friend void shift_left(gf2n_polynomial &erg, const gf2n_polynomial &a,
			       gf2n_poly_modulus &F)
	{
		shift_left(erg, a, F, 1);
	}


	// compute a^i mod F.mod for i = 0, ..., d-1, return table
	// Note that deallocation has to be done by user !!

	friend gf2n_polynomial *power_table(const gf2n_polynomial & a,
					    gf2n_poly_modulus & F, int d);



	friend void Xq(gf2n_polynomial &, gf2n_poly_modulus &,
		       unsigned int d); // X^(2^d)

	friend void power(gf2n_polynomial &, const gf2n_polynomial &,
			  const bigint &, gf2n_poly_modulus &);


	// computes g = f(X) with t[i] = X^i mod modulus (reduced)

	friend void compose(gf2n_polynomial & g, const gf2n_polynomial & f,
			    const gf2n_polynomial & X,
			    const gf2n_polynomial * t);

	// computes g = f(X) mod F.mod, reductions with gf2n_poly_modulus

	friend void compose(gf2n_polynomial & g, const gf2n_polynomial & f,
			    const gf2n_polynomial & X,
			    gf2n_poly_modulus &);


	// Xq is X^q mod F.mod. The function computes the degree of all the
	// irreducible factors of F.mod (which is assumed to factor (d ... d))
	// and returns it.
	// d is assumed a multiple of the common degree,
	// if d == -1, use degree of F.mod as multiple

	friend lidia_size_t compute_degree(const gf2n_polynomial & Xq,
					   gf2n_poly_modulus & F,
					   lidia_size_t d);

};

lidia_size_t compute_degree(const gf2n_polynomial & Xq,
			    gf2n_poly_modulus & F,
			    lidia_size_t d = -1);


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_GF2N_POLY_MODULUS__H_GUARD_

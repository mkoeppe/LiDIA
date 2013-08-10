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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/Fp_poly_multiplier.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//******************************************************************
//
//			class poly_argument
//
//*******************************************************************

lidia_size_t poly_argument::poly_arg_bound = 0;
bool poly_argument::change_enable = true;

void poly_argument::set_poly_arg_bound(lidia_size_t n)
{
	//this value can only be changed if no poly_arguments have been
	//declared yet
	debug_handler("poly_argument", "set_poly_arg_bound(lidia_size_t)");
	if (change_enable == true) {
		poly_arg_bound = n;
		change_enable = false;
	}
	else
		lidia_error_handler("poly_argument", "set_poly_arg_bound"
				    "(lidia_size_t)::poly_arg_bound can only be set before "
				    "declaring any poly_arguments");
}



// The routine poly_argument::build (see below) which is implicitly called
// by the various compose and update_map routines builds a table
// of polynomials.  If poly_arg_bound > 0, then only up to
// about poly_arg_bound bigmod's are allocated.  Setting this too
// low may lead to slower execution.
// If poly_arg_bound <= 0, then it is ignored, and space is allocated
// so as to maximize speed.
// Initially, poly_arg_bound = 0.

// If a single h is going to be used with many g's
// then you should build a poly_argument for h,
// and then use the compose routine below.
// poly_argument::build computes and stores h, h^2, ..., h^m mod f.
// After this pre-computation, composing a polynomial of degree
// roughly n with h takes n/m multiplies mod f, plus n^2
// scalar multiplies.
// Thus, increasing m increases the space requirement and the pre-computation
// time, but reduces the composition time.
// If poly_arg_bound > 0, a table of size less than m may be built.

// m must be > 0, otherwise an error is raised

void poly_argument::build(const Fp_polynomial& h, const Fp_poly_modulus& F, lidia_size_t m)
{
	debug_handler("poly_argument", "build (Fp_polynomial&, Fp_poly_modulus&, lidia_size_t)");

	const Fp_polynomial &f = F.modulus();
	h.comp_modulus(f, "poly_argument::build");

	if (m <= 0) {
		lidia_error_handler("poly_argument", "build (Fp_polynomial&, "
				    "Fp_poly_modulus&, lidia_size_t)::bad args");
		return;
	}


	if (poly_arg_bound > 0) {
		m = comparator< lidia_size_t >::min(m, poly_arg_bound/f.degree());
		m = comparator< lidia_size_t >::max(m, 1);
	}

	Fp_poly_multiplier M;
	if (h.degree() < f.degree())
		M.build(h, F);
	else {
		Fp_polynomial tmp;
		remainder(tmp, h, F);
		M.build(tmp, F);
	}

	if (H.capacity() < m+1)
		H.set_capacity(m+1);
	H.set_size(m+1);

	H[0].set_modulus(h);
	(H[0]).assign_one();
	H[1].assign(h);

	lidia_size_t i;
	for (i = 2; i <= m; i++)
		multiply(H[i], H[i-1], M, F); //Fp_poly_multiplier
}



//***************************************************************
//
//                      compose
//
//****************************************************************

void poly_argument::compose(Fp_polynomial& x, const Fp_polynomial& g,
			    const Fp_poly_modulus& F) const
{
	debug_handler("poly_argument", "compose(Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	g.comp_modulus(F.modulus(), "poly_argument::compose");

	if (g.degree() <= 0) {
		x.assign(g);
		return;
	}

	Fp_polynomial s, t;
	lidia_size_t n = F.modulus().degree();

	bigint *scratch = new bigint[n];
	memory_handler(scratch, "poly_argument",
		       "compose :: Error in memory allocation");

	lidia_size_t m = H.size() - 1;
	lidia_size_t l = ((g.degree()+m)/m) - 1;

	Fp_poly_multiplier M(H[m], F);

	compose_InnerProd(t, g, l*m, l*m + m - 1, H, n, scratch);

	lidia_size_t i;
	for (i = l-1; i >= 0; i--) {
		compose_InnerProd(s, g, i*m, i*m + m - 1, H, n, scratch);
		multiply(t, t, M, F); //Fp_poly_multiplier
		add(t, t, s);
	}

	delete[] scratch;
	x.assign(t);
}



void poly_argument::compose_InnerProd(Fp_polynomial& x,
				      const Fp_polynomial &v, lidia_size_t low, lidia_size_t high,
				      const base_vector< Fp_polynomial > & H, lidia_size_t n, bigint * t)
	// only used in compose
	// assumes t.size() >= n
{
	debug_handler("poly_argument", "compose_InnerProduct(Fp_polynomial&, Fp_polynomial&, lidia_size_t, lidia_size_t, base_vector< Fp_polynomial > &, lidia_size_t, bigint*)");

	if (H.size() == 0) {
		lidia_error_handler("poly_argument",
				    "compose_InnerProduct(...)::wrong argument size");
		return;
	}

	bigint s; //static
	const bigint &p = H[0].modulus();
	lidia_size_t i, j;

	for (j = 0; j < n; j++)
		t[j].assign_zero();

	high = comparator< lidia_size_t >::min(high, v.degree());
	for (i = low; i <= high; i++) {
		const bigint* h = H[i-low].coeff;
		lidia_size_t m = H[i-low].c_length;

		const bigint & w = v.coeff[i];

		for (j = 0; j < m; j++) {
			multiply(s, w, h[j]);
			add(t[j], t[j], s);
		}
	}

	x.set_modulus(p);
	x.set_degree(n-1);
	for (j = 0; j < n; j++)
		Remainder(x.coeff[j], t[j], p);
	x.remove_leading_zeros();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/Fp_poly_modulus.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void append_irred_factor(factorization< Fp_polynomial > &F,
			 const Fp_polynomial & f, lidia_size_t i)
{
	single_factor< Fp_polynomial > tmp(f);
	if (f.degree() > 0)
		tmp.set_prime_flag(decomposable_object::prime);
	F.append(tmp, i);
}



//
// Task:	compute the factorization of f
//
// Conditions:	f must be a monic polynomial of degree 2,
//		input may not alias output

void
factor_quadratic_pol(factorization< Fp_polynomial > & F, const Fp_polynomial &f)
{
	debug_handler("Fp_polynomial", "factor_quadratic_pol(factorization< Fp_polynomial > &, Fp_polynomial&)");

	if (f.degree() != 2 || !f.is_monic()) {
		lidia_error_handler("Fp_polynomial", "factor_quadratic_pol()::"
				    "input must be a monic polynomial of degree 2");
		return; // LC
	}

	bool verbose = single_factor< Fp_polynomial >::verbose();
	const bigint &p = f.modulus();

	if (verbose) std::cout << "factoring quadratic polynomial" << std::endl;

	F.kill();
	Fp_polynomial g;

	if (f.const_term().is_zero()) {
		shift_right(g, f, 1);
		append_irred_factor(F, g);
		g.assign_x();
		append_irred_factor(F, g);
		return;
	}

	bigint t;
	square(t, f[1]);
	subtract(t, t, f[0]*4); // t = b^2 - 4a

	if (t.is_zero() || jacobi(t, p) == 1 || (p == 2 && f[1].is_zero())) {
		if (verbose) std::cout << " splitting f in 2 linear factors " << std::endl;
		t.assign(find_root(f));
		t.negate();

		g.set_modulus(p);
		g.assign_x();
		g.set_coefficient(t, 0);
		append_irred_factor(F, g);

		divide(g, f, g);
		append_irred_factor(F, g);
		F.normalize();
	}
	else {
		if (verbose) std::cout << " f is irreducible with degree 2 " << std::endl;
		append_irred_factor(F, f);
	}
}



// Literature:	Th. Pfahler, Diplomarbeit, Universitaet des Saarlandes, 1996
//
// Task:	computes the factorization of f
//
// Conditions:	f is assumed to be monic, non-zero, square-free, f(0)!=0
//		f may not alias a component of Fact
//
// Algorithm:	uses some heuristics to determine which algorithm is
//		probably the fastest
//

static
void sf_factor(factorization< Fp_polynomial > &Fact, const Fp_polynomial& f)
{
	debug_handler("Fp_polynomial", "sf_factor(factorization< Fp_polynomial > &, Fp_polynomial&)");

	bool verbose = single_factor< Fp_polynomial >::verbose();
	lidia_size_t n = f.degree();
	Fact.kill();

	// test for binomial/trinomial
	const bigint &p = f.modulus();
	lidia_size_t i, non_zeros = 2; //first non-zero coefficients: f(n), f(0)

	for (i = 1; i < n; i++)
		if (!f[i].is_zero()) {
			non_zeros++;
			if (non_zeros > 3) break;
		}

	switch(non_zeros) {
	case(2) : Fact.assign(factor_binomial(f));
		if (Fact.no_of_components() > 1
		    || Fact.no_of_prime_components() == 1) return;

		Fact.kill(); // otherwise, we couldn't split f
		break; // and have to use the normal routines

	case(3) : // this is a very simple test
		// if  f = x^p - x - a,  a!=0  =>  f irred.
		if (n == p && f[1] == p-1) {
			if (verbose) std::cerr << "found irred. trinomial" << std::endl;
			append_irred_factor(Fact, f);
			return;
		}
		break; // otherwise, use normal routines
	}

	// unfortunately, we must work harder for the rest...

	my_timer t;
	Fp_polynomial x_to_the_p, x;
	x.set_modulus(f);
	x.assign_x();
	Fp_poly_modulus F(f);

	if (verbose) t.start("computing x^p...");
	power_x(x_to_the_p, p, F);
	if (verbose) t.stop();



	// splitting linear factors

	if (verbose) t.start("splitting linear factors");
	Fp_polynomial d1, g1, ff;
	subtract(g1, x_to_the_p, x);
	gcd(d1, g1, f); // d1 = gcd(x^p-x, f)
	if (!d1.is_one()) {
		if (verbose)
			std::cerr << " root_edf, number of factors = " << d1.degree() << std::endl;

		root_edf(Fact, d1);

		divide(ff, f, d1);
		if (ff.degree() < 1) {
			if (verbose) t.stop(); // no factors left !
			return;
		}

		F.build(ff);
		remainder(x_to_the_p, x_to_the_p, F);
	}
	else
		ff.assign(f); // f has no factors of degree 1
	if (verbose) t.stop();



	factorization< Fp_polynomial > Fact2;
#if 0
	// splitting factors of degree 2
	Fp_polynomial d2, g2;
	if (verbose) t.start("splitting factors of degree 2");
	power_compose(g2, x_to_the_p, 2, F);
	subtract(g2, g2, x);
	gcd(d2, g2, ff); // d2 = gcd(x^{p^2}-x, f)
	if (!d2.is_one()) {
		if (verbose) std::cerr << " edf, degree 2, number of factors = "
				       << d2.degree()/2 << std::endl;
		F.build(d2);
		Fp_polynomial tmp;
		remainder(tmp, x_to_the_p, F);

		edf(Fact2, F, tmp, 2);
		multiply(Fact, Fact, Fact2);

		divide(ff, ff, d2);
		if (ff.degree() < 1) {
			if (verbose) t.stop(); // no factors left !
			return;
		}
		F.build(ff);
		remainder(x_to_the_p, x_to_the_p, F);
	}
	if (verbose) t.stop();
#endif



	// general case

	if (verbose) t.start("general case");

	lidia_size_t bl = p.bit_length(); // considering space complexity
	if (n*bl< 100000 || n/bl > 20)
		sf_berlekamp_work(Fact2, x_to_the_p, F);
	else {
		// maybe we should set poly_argument::poly_arg_bound
		sf_can_zass_work(Fact2, x_to_the_p, F);
	}

	multiply(Fact, Fact, Fact2);

	if (verbose) t.stop();
}



//***********************************************************************
//
//			    "interface"
//
//***********************************************************************

void
factor(factorization< Fp_polynomial > &F, const Fp_polynomial& f)
{
	debug_handler("Fp_polynomial", "factor(factorization< Fp_polynomial > &, Fp_polynomial&)");

	factor_generic(F, f, sf_factor);
}



factorization< Fp_polynomial >
single_factor< Fp_polynomial >::factor() const
{
	debug_handler("single_factor< Fp_polynomial >", "factor()");

	factorization< Fp_polynomial > F;
	LiDIA::factor(F, rep);
	return F;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

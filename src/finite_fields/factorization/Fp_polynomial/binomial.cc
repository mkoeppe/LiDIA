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
#include	"LiDIA/rational_factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// IMPORTANT:
//
// Everything in this file is heavily based on
// Lidl, Niederreiter: Introduction to finite fields and their applications
// Cambridge University Press, 1986
//
// In this file, this book is referred to as 'LN'.
//




//
// Task:	returns multiplicative order of a in Z/pZ
// 		returns -1 if factorization failed
//
// Conditions:	p must be prime and positive

static
bigint order_mod(const bigint &aa, const bigint &p, int verbose = 0)
{
	debug_handler("binomial.c", "order_mod(bigint&, bigint&, int)");

	bigint a(aa);
	if (a.is_negative() || a >= p)
		remainder(a, a, p);

	if (a.is_zero())
		return bigint(0);
	if (a.is_one())
		return bigint(1);

	bigint m = p-1;
	rational_factorization f(m);

	my_timer t;

	if (verbose)
		t.start("factoring p-1 : ");
	if (m < 1000000) {
		int mm;
		m.intify(mm);
		mm >>= 1; // since p-1 is (almost always) divisible by 2...
		if (mm < 5)
			mm = 5;
		f.trialdiv(mm);
	}
	else
		f.factor();
	if (verbose)
		t.stop();

	if (!f.is_prime_factorization()) {
		lidia_warning_handler("binomial.c", "order_mod(bigint&, bigint&)::"
				      "could not factor p-1");
		return bigint(-1);
	}

	bigint p_i, pow;
	lidia_size_t i, e_i, end = f.no_of_comp();
	for (i = 0; i < end; i++) {
		p_i = f.base(i);
		e_i = f.exponent(i);
		power_mod(pow, a, (m/p_i), p);
		while ((e_i > 0) && (pow.is_one())) {
			m = m/p_i;
			e_i--;
			power_mod(pow, a, (m/p_i), p);
		}
	}
	return m;
}




//***********************************************************************
//
//	special algorithms for factoring binomials
//	- case1_factor_binomial, case2_factor_binomial,
//	case3_factor_binomial, case4_factor_binomial
//	are "sub"cases
//	- factor_binomial
//	is the main function
//	GENERAL ASSUMPTION: f must be a monic binomial, f(0) != 0
//	========
//
//***********************************************************************


//
// Literature:	LN, Theorem 3.75
//
// Task:	returns a factor of f = x^t - a (mod p)
//
// Conditions:	r must be a prime factor of t, dividing d = (p-1)/ord(a)
//		=> a = b^r  for some b,   x^t - a = x^(t1*r) - b^r
//
// Algorithm:	also correct for GF(p^n)
//

static
Fp_polynomial
case1_factor_binomial(const bigint & a, lidia_size_t t, lidia_size_t r, const bigint & p)
{
	debug_handler("binomial.c", "case1_factor_binomial(bigint&, lidia_size_t, lidia_size_t, bigint&)");

	lidia_size_t t1 = t / r;
	bigint b;

	//compute b so that a = b^r
	bool verbose = single_factor< Fp_polynomial >::verbose();

	if (verbose)
		std::cerr << "Compute " << r << "-th root..." << std::endl;
	Fp_polynomial g;
	g.set_modulus(p);
	g.set_coefficient(r);
	g.set_coefficient(-a, 0);
	b.assign(find_root(g));
	if (verbose)
		std::cerr << "done." << std::endl;


	// x^t1 - b  is a divisor of f
	Fp_polynomial factor;
	factor.set_modulus(g);
	factor.set_coefficient(t1);
	if (r % 2 != 0)
		b.negate(); //if r is even, -b is also an r-th root
	factor.set_coefficient(b, 0); //was actually: factor[0] = -b
	return factor;
}



//
// Literature:	LN, Theorem 3.75
//
// Task:	returns a factor of f = x^t - a (mod p)
//
// Conditions:	r must be a prime factor of t, not dividing e = ord(a)
//		=> r1*r = 1 mod (p-1)  for some r1,
//		   x^t - a = x^(t1*r) - a^(r1*r)
//

static
Fp_polynomial
case2_factor_binomial(const bigint & a, lidia_size_t t, lidia_size_t r, const bigint & p)
{
	debug_handler("binomial.c", "case2_factor_binomial(bigint&, lidia_size_t, lidia_size_t, bigint&)");

	lidia_size_t t1 = t / r;
	bigint r1;

	InvMod(r1, r, p - 1);
	bigint b;
	power_mod(b, a, r1, p); //b = a^r1  mod p

//x^t1 - a^r1  is a divisor of f
	Fp_polynomial factor;
	factor.set_modulus(p);
	factor.set_coefficient(t1);
	factor.set_coefficient(-b, 0);
	return factor;
}



//
// Literature:	LN, Theorem 3.75
//
// Task:	returns a factor of f = x^t - a (mod p)
//
// Conditions:	all prime factors of t divide e = ord(a) but not d = (p-1)/e
//		t == 0 mod 4, p != 1 mod 4
//		e divides (p-1) ==> p odd ==> p == 3 mod 4
//		==> e = 2 mod 4
//
// Algorithm:	also correct for GF(p^n)
//

static
Fp_polynomial
case3_factor_binomial(const bigint & a, lidia_size_t t, const bigint & e, const bigint & p)
{
	debug_handler("binomial.c", "case3_factor_binomial(bigint&, lidia_size_t, bigint&, bigint&)");

	bigint d = (e >> 1) + 1;
	bigint d_half = d >> 1; //d_half = d/2
	bigint c, tmp;

	power_mod(c, a, d_half, p);
	InvMod(tmp, 2, p);
	MulMod(c, c, tmp, p); //c = a^d_half * 2^(-1)

	add(tmp, p, 1);
	shift_right(tmp, tmp, 2); //tmp = (p+1)/4

	power_mod(c, c, tmp, p); //c = (a^d_half * 2^(-1)) ^ ((p+1)/4)
	//a^d = 4c^4   since a^2 == a^(p+1)

	Fp_polynomial factor;
	factor.set_modulus(p);
	lidia_size_t t2 = t >> 2; //t2 = t/4;

	//x^t - a = x^t + a^(e/2 + 1)
	//        = x&t + a^d
	//        = x^(4*t2) + 4c^4
	//        = (x^(2*t2) + 2*c*x^t2 + 2*c^2) * (x^(2*t2) - 2*c*x^t2 + 2c^2)
	//        = factor * g

	factor.set_coefficient(2 * t2);
	square(tmp, c);

	tmp.multiply_by_2(); //tmp = c^2
	factor.set_coefficient(tmp, 0);

	c.multiply_by_2(); //c = 2c
	factor.set_coefficient(c, t2);

	return factor;
}



//
// Literature:	LN, Theorem 3.76
//
// Task:	returns the factorization of f = x^t - a (mod p)
//
// Conditions:	e = ord(a), p == 3 mod 4,
//		all prime factors of t divide e, but not d = (p-1)/e
//		p = 2^A * u - 1, u odd,  2^A must divide t
//
// Algorithm:	also correct for GF(p^n)
//

static
factorization< Fp_polynomial >
case4_factor_binomial(const bigint & a, lidia_size_t t, lidia_size_t A, const bigint & e, const bigint & p)
{
	debug_handler("binomial.c", "case4_factor_binomial(bigint&, lidia_size_t, bigint&, bigint&)");

	bool verbose = single_factor< Fp_polynomial >::verbose();

	lidia_size_t B = 1 << (A - 1);
	lidia_size_t v = t >> (A - 1); //t = B * v

	if (verbose) {
		std::cerr << "f factors as a product of " << B << " monic irreducible ";
		std::cerr << "polynomials of degree " << v << std::endl;
	}

	Fp_polynomial g;
	g.set_modulus(p);
	lidia_size_t i;
	bigint c;

	if (verbose)
		std::cerr << "Computing special polynomial" << std::endl;

	//build g
	c.assign_one();
	// c(0) = 1
	// c(i) = ( (B-i-1)! B )/( i! (B-2i)! )
	// c(i) = c(i-1) * ( (B-2i+2)(B-2i+1) )/( (i)(B-i) )

	g.set_coefficient(B);
	for (i = 1; i <= (B / 2); i++) {
		c *= (B - 2 * i + 2) * (B - 2 * i + 1);
		c /= i * (B - i);
		g.set_coefficient(c, B - 2 * i);
	}

	//find roots of g
	if (verbose)
		std::cerr << "Finding roots" << std::endl;
	base_vector< bigint > roots;
	rec_find_roots(roots, g);

	//build complete factorization
	if (verbose)
		std::cerr << "Build complete factorization" << std::endl;
	bigint r;
	bigint tmp, D;

	shift_right(D, e, 1);
	inc(D); //D = (e/2)+1, D is even !

	//compute r s.t.
	// 2Br == D  mod (p-1)
	//we know that gcd(2B,p-1) = gcd(B,p-1) = 2

	xgcd_left(tmp, 2 * B, p - 1); //2 == tmp*(2B)  mod (p-1)

	MulMod(r, tmp, D, p - 1); //have 2r == D  mod (p-1)
	if ((r & 1) != 0)
		lidia_error_handler("factor_binomial", "d is not divisible by 2 !!!");
	r = r >> 1;

	bigint b;
	power_mod(b, a, r, p); //b = a^r

	MulMod(tmp, b, b, p); //square_mod


	//factor = x^v - b*roots[i]*x^(v/2) - b^2
	if (verbose)
		std::cerr << "Add factors..." << std::endl;

	factorization< Fp_polynomial > F;
	Fp_polynomial factor;
	factor.set_modulus(p);
	factor.set_coefficient(v);
	tmp.negate();
	factor.set_coefficient(tmp, 0);
	b.negate();

	for (i = 0; i < B; i++) {
		MulMod(tmp, b, roots[i], p);
		factor.set_coefficient(tmp, v / 2);
		append_irred_factor(F, factor);
	}
	return F;
}



//***********************************************************************
//
//			the main routine
//
//***********************************************************************

//
// Task:	returns a factorization of the binomial f
//
// Conditions:	f must be a monic binomial with f(0)!=0
//
// Algorithm:	also correct for GF(p^n)
//

factorization< Fp_polynomial >
factor_binomial(const Fp_polynomial & f, int FACTOR_P_MINUS_ONE)
{
	debug_handler("binomial.c", "factor_binomial(Fp_polynomial&, int)");

	bool verbose = single_factor< Fp_polynomial >::verbose();
	factorization< Fp_polynomial > F;

	if (!f.is_binomial() || f.const_term().is_zero()) {
		if (verbose)
			std::cerr << "factor_binomial::input must be a binomial with nonzero "
				"constant term" << std::endl;
		F.assign(f);
		return F;
	}

	if (verbose)
		std::cerr << "factor_binomial " << f << std::endl;

	const bigint & p = f.modulus();
	bigint p_minus_one(p);
	dec(p_minus_one);

	lidia_size_t i;
	lidia_size_t deg = f.degree();

	switch(deg) {
	case(1) :
		if (verbose)
			std::cerr << "polynomial has degree 1 =  > irreducible" << std::endl;
		append_irred_factor(F, f);
		return F;
	case(2) :
		factor_quadratic_pol(F, f);
		return F;
	}

	//now, deg > 2

	Fp_polynomial g, h;
	g.set_modulus(p);

	if (deg == p_minus_one && f.const_term() == p_minus_one) {
		//f = (x^(p-1) - 1)
		if (verbose)
			std::cerr << "f = (x^(" << p << "-1) - 1) mod " << p << std::endl;
		g.assign_x();
		lidia_size_t small_p;
		p.sizetify(small_p);
		for (i = 1; i < small_p; i++) {
			g.set_coefficient(i, 0);
			append_irred_factor(F, g);
		}
		return F;
	}


	lidia_size_t r;
	bigint a = -f.const_term();

	//f = x^t - a, t>=2, a!=0

	fac_vec t_fac(deg); //t_fac contains the factorization of deg
	int ii, num_factors = t_fac.number_of_factors();

	if (f.const_term() == p_minus_one &&
	    !(num_factors == 1 && t_fac[0].a == 1)) {
		//f = x^t - 1, t not prime
		g.assign_one();
		g.negate();
		F.append(f);
		for (ii = 0; ii < num_factors; ii++) {
			if (verbose)
				std::cerr << "Refinement with x^"
					  << deg/t_fac[ii].q << " - 1" << std::endl;

			g[deg/t_fac[ii].q].assign_one();
			F.refine(g);
			g[deg/t_fac[ii].q].assign_zero();
		}
		return F;
	}


	if (p.length() > 5 && !FACTOR_P_MINUS_ONE) {
//	if (verbose) std::cerr << "avoid factorization of (modulus-1)" << std::endl;
		F.append(f);
		return F;
	}

	bigint e, d;
	e = order_mod(a, p); //order of a in GF(p)
	if (e == -1) {
		if (verbose)
			std::cerr << "factorization of (modulus-1) failed" << std::endl;
		F.append(f);
		return F;
	}


	divide(d, p - 1, e);
	long rem;

	if (verbose) {
		std::cerr << "Computing prime factors of the degree and of the order ";
		std::cerr << "of the const term..." << std::endl;
		std::cerr << "deg = " << deg << "   e = " << e << std::endl;
	}


	lidia_size_t count = 0;
	bool hit;
	for (ii = num_factors-1; ii >= 0; ii--) {
		hit = false;
		r = t_fac[ii].q;

		remainder(rem, d, r);
		if (rem == 0) {
			//if (! r does not divide (p-1)/e)
			hit = true;
			if (verbose) {
				std::cerr << "prime factor " << r << " of " << deg << " divides ";
				std::cerr << "(p-1)/e = " << d << std::endl;
				std::cerr << " =  > split into 2 factors" << std::endl;
			}
			g = case1_factor_binomial(a, deg, r, p);
			if (count == 0) {
				divide(h, f, g);
				F.append(g);
				F.append(h);
				count++;
			}
			else
				F.refine(g);

//	    if (verbose)
//	    	std::cerr << "cannot split f (cannot take roots in Fp)" << std::endl;
//	    F.append(f);
		}

		if (!hit) {
			remainder(rem, e, r);
			if (rem != 0) {
				//if (! r divides e)
				if (verbose) {
					std::cerr << "prime factor " << r << " of " << deg;
					std::cerr << " does not divide the order " << e << std::endl;
					std::cerr << " =  > split into 2 factors" << std::endl;
				}
				g = case2_factor_binomial(a, deg, r, p);
				if (count == 0) {
					divide(h, f, g);
					F.append(g);
					F.append(h);
					count++;
				}
				else
					F.refine(g);
			}
		}
	}

	if (count != 0) {
		F.refine();
		return F;
	}

	//now, all prime factors of t divide e, but not d
	if (verbose)
		std::cerr << "testing second condition..." << std::endl;

	if (deg % 4 != 0) {
		if (verbose)
			std::cerr << "deg != 0 mod 4 =  > f is irreducible" << std::endl;
		append_irred_factor(F, f);
		return F;
	}

	//now, deg == 0 mod 4
	if (p % 4 == 1) {
		if (verbose)
			std::cerr << "deg % 4 == 0  and  modulus % 4 == 1 =  >"
				"f is irreducible" << std::endl;

		append_irred_factor(F, f);
		return F;
	}

	//now, deg == 0 mod 4 and p == 3 mod 4
	if (verbose)
		std::cerr << "deg % 4 == 0  and  modulus % 4 == 3" << std::endl;
	bigint tmp(p);
	inc(tmp);
	lidia_size_t A = 0;
	while (tmp.is_even()) {
		A++;
		tmp.divide_by_2();
	}

	if (verbose) {
		std::cerr << "deg = " << deg << "   p = " << p << "   A = " << A << std::endl;
		std::cerr << "p = 2^A*u - 1" << std::endl;
	}

	if (deg % (1 << A) == 0) {
		//if (deg is divisible by 2^A)
		if (verbose)
			std::cerr << "deg is divisible by 2^A =  > compute canonical "
				"factorization" << std::endl;

		F = case4_factor_binomial(a, deg, A, e, p);
	}
	else {
		if (verbose)
			std::cerr << "deg is not divisible by 2^A =  > split into "
				"2 factors" << std::endl;

		g = case3_factor_binomial(a, deg, e, p);
		F.append(g);
		g.set_coefficient(-g[deg / 4], deg / 4);
		//negate coeff. of x^(deg/4), see case3_factor_binomial
		F.append(g);
	}
	return F;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

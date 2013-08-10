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
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/Fp_poly_modulus.h"
#include	"LiDIA/Fp_poly_multiplier.h"




#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool plain_eq_frac(const Fp_polynomial& a, const Fp_polynomial& b,
		   const Fp_polynomial& c, const Fp_polynomial& d, const Fp_poly_modulus& F)
{
	debug_handler("fractions.c", "plain_eq_frac(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	Fp_polynomial t1, t2;

	multiply(t1, a, d);
	multiply(t2, b, c);
	subtract(t1, t1, t2);
	remainder(t2, t1, F);

	return t2.is_zero();
}



void plain_add_frac(Fp_polynomial& x, Fp_polynomial& y, const Fp_polynomial& a,
		    const Fp_polynomial& b, const Fp_polynomial& c, const Fp_polynomial& d,
		    const Fp_poly_modulus& F)
{
	debug_handler("fractions.c", "plain_add_frac(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	Fp_polynomial t1, t2;

	multiply(t1, a, d);
	multiply(t2, b, c);
	add(t1, t1, t2);
	remainder(x, t1, F);
	multiply(y, b, d, F); //Fp_poly_modulus
}



void plain_subtract_frac(Fp_polynomial& x, Fp_polynomial& y,
			 const Fp_polynomial& a, const Fp_polynomial& b,
			 const Fp_polynomial& c, const Fp_polynomial& d, const Fp_poly_modulus& F)
{
	debug_handler("fractions.c", "plain_subtract_frac(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	Fp_polynomial t1, t2;

	multiply(t1, a, d);
	multiply(t2, b, c);
	subtract(t1, t1, t2);
	remainder(x, t1, F);
	multiply(y, b, d, F); //Fp_poly_modulus
}



void remainder(Fp_polynomial& x, fft_rep& R1, const Fp_poly_modulus& F,
	       modular_fft_rep& R2, Fp_polynomial& P1)
{
	debug_handler("fractions.c", "remainder(Fp_polynomial&, fft_rep&, Fp_poly_modulus&, modular_fft_rep&, Fp_polynomial&)");

	lidia_size_t n = F.n;
	lidia_size_t d = 2*n-2;
	lidia_size_t index;

	modular_fft_rep R3;
	R3.init(F.fd);
	R3.set_size(F.k);

	R1.from_fft_rep(P1, n, d); // save R1 for future use

	R2.set_size(F.l);
	for (index = 0; index < F.fd.number_of_primes(); index++) {
		R2.to_modular_fft_rep(P1, index);
		multiply(R2, F.HRep, R2, index);
		R2.from_modular_fft_rep(n-2, 2*n-4, index);
	}
	R2.get_result(P1, n-2, 2*n-4);

	R2.set_size(F.k);
	for (index = 0; index < F.fd.number_of_primes(); index++) {
		R2.to_modular_fft_rep(P1, index);
		multiply(R2, F.FRep, R2, index);
		reduce(R3, R1, F.k, index);
		subtract(R2, R3, R2, index);
		R2.from_modular_fft_rep(0, n-1, index);
	}
	R2.get_result(x, 0, n-1);
}



void subtract_frac(Fp_polynomial& x, Fp_polynomial& y, const Fp_polynomial& a,
		   const Fp_polynomial& b, const Fp_polynomial& c,
		   const Fp_polynomial& d, const Fp_poly_modulus& F)
{
	debug_handler("fractions.c", "subtract_frac(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	lidia_size_t m, n, k, index;

	if (!F.use_FFT) {
		plain_subtract_frac(x, y, a, b, c, d, F);
		return;
	}

	n = F.n;
	m = 2*n-2;
	k = next_power_of_two(m+1);


	fft_rep R_a, R_b;
	R_a.init(F.fd);
	R_b.init(F.fd);
	R_a.set_size(k);
	R_b.set_size(k);

	modular_fft_rep R1(F.fd), R2(F.fd), R3(F.fd), R4(F.fd);

	for (index = 0; index < F.fd.number_of_primes(); index++) {
		R1.to_modular_fft_rep(b, index);
		R2.to_modular_fft_rep(d, index);
		R3.to_modular_fft_rep(a, index);
		R4.to_modular_fft_rep(c, index);

		multiply(R3, R3, R2, index);
		multiply(R4, R4, R1, index);

		subtract(R3, R3, R4, index); // numerator
		multiply(R4, R1, R2, index); // denominator

		reduce(R_a, R3, k, index); // = copy
		reduce(R_b, R3, k, index); // = copy
	}

	Fp_polynomial P1;
	P1.set_max_degree(n-1);

	remainder(x, R_a, F, R1, P1);
	remainder(y, R_b, F, R1, P1);
}



void add_frac(Fp_polynomial& x, Fp_polynomial& y, const Fp_polynomial& a,
	      const Fp_polynomial& b, const Fp_polynomial& c,
	      const Fp_polynomial& d, const Fp_poly_modulus& F)
{
	debug_handler("fractions.c", "AddFrac(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	lidia_size_t m, n, k, index;

	if (!F.use_FFT) {
		plain_add_frac(x, y, a, b, c, d, F);
		return;
	}

	n = F.n;
	m = 2*n-2;
	k = next_power_of_two(m+1);


	fft_rep R_a, R_b;
	R_a.init(F.fd);
	R_b.init(F.fd);
	R_a.set_size(k);
	R_b.set_size(k);
	modular_fft_rep R1(F.fd), R2(F.fd), R3(F.fd), R4(F.fd);

	for (index = 0; index < F.fd.number_of_primes(); index++) {
		R1.to_modular_fft_rep(b, index);
		R2.to_modular_fft_rep(d, index);
		R3.to_modular_fft_rep(a, index);
		R4.to_modular_fft_rep(c, index);

		multiply(R3, R3, R2, index);
		multiply(R4, R4, R1, index);

		add(R3, R3, R4, index); // numerator
		multiply(R4, R1, R2, index); // denominator

		reduce(R_a, R3, k, index); // = copy
		reduce(R_b, R3, k, index); // = copy
	}

	Fp_polynomial P1;
	P1.set_max_degree(n-1);

	remainder(x, R_a, F, R1, P1);
	remainder(y, R_b, F, R1, P1);
}



bool eq_frac(const Fp_polynomial& a, const Fp_polynomial& b,
	     const Fp_polynomial& c, const Fp_polynomial& d, const Fp_poly_modulus& F)
{
	debug_handler("fractions.c", "eq_frac(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	lidia_size_t m, n, k;

	if (!F.use_FFT)
		return plain_eq_frac(a, b, c, d, F);

	n = F.n;
	m = 2*n-2;
	k = next_power_of_two(m+1);


	fft_rep R_(F.fd);
	R_.set_size(k);
	modular_fft_rep R1(F.fd), R2(F.fd), R3(F.fd);

	lidia_size_t index;
	for (index = 0; index < F.fd.number_of_primes(); index++) {
		R1.to_modular_fft_rep(a, index);
		R2.to_modular_fft_rep(d, index);
		multiply(R1, R1, R2, index);

		R2.to_modular_fft_rep(b, index);
		R3.to_modular_fft_rep(c, index);
		multiply(R2, R2, R3, index);

		subtract(R1, R1, R2, index);
		reduce(R_, R1, k, index); // = copy
	}

	Fp_polynomial P1;
	P1.set_max_degree(n-1);
	remainder(P1, R_, F, R2, P1);

	return P1.is_zero();
}



//***********************************************************************
//
//			    class eq_frac_info
//
//***********************************************************************

void eq_frac_info::build(const Fp_polynomial& a, const Fp_polynomial& b,
			 const Fp_poly_modulus& F)
{
	debug_handler("eq_frac_info", "build(Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	a.comp_modulus(b, "eq_frac_info::build");
	a.comp_modulus(F.modulus(), "eq_frac_info::build");

	modulus = a.modulus();

	lidia_size_t n = F.modulus().degree();
	base_vector< bigint > L(n, n);
	lidia_size_t i;

	for (i = 0; i < n; i++)
		L[i].assign(randomize(a.modulus()));

	Fp_poly_multiplier M;

	M.build(a, F);
	update_map(L1, L, M, F);

	M.build(b, F);
	update_map(L2, L, M, F);
}



bool eq_frac_info::eq_frac(const Fp_polynomial& c, const Fp_polynomial& d) const
{
	debug_handler("eq_frac_info", "eq_frac(Fp_polynomial&, Fp_polynomial&)");
	c.comp_modulus(d, "eq_frac_info::eq_frac");
	if (c.modulus() != modulus) {
		lidia_error_handler("eq_frac_info", "eq_frac(...)::different moduli");
		return false;
	}

	bigint x, y;

	inner_product(x, L1, d);
	inner_product(y, L2, c);

	return x == y;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

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
#include	"LiDIA/Fp_poly_modulus.h"
#include	"LiDIA/Fp_poly_multiplier.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//******************************************************************
//
//			  class Fp_poly_modulus
//
//******************************************************************

// If you need to do a lot of arithmetic modulo a fixed f,
// build Fp_poly_modulus F for f.  This pre-computes information about f
// that speeds up the computation a great deal.
// f should be monic, and deg(f) > 0.

void Fp_poly_modulus::build(const Fp_polynomial& pol)
{
	debug_handler("Fp_poly_modulus", "build(Fp_polynomial&)");

	//if ( f == pol )
	//	return; 	//initializing with the same polynomial again

	bool re_init = (f.degree() > 0);

	f.assign(pol);
	n = pol.degree();

	if (n <= 0 || !pol.is_monic()) {
		lidia_error_handler("Fp_poly_modulus", "build(Fp_polynomial&)::"
				    "modulus must be monic and have degree >= 1");
		return;
	}

	const bigint &p = pol.modulus();
	crov = Fp_polynomial::crossovers.fftmul_crossover(p);
	if (n <= crov) {
		use_FFT = false;
		return;
	}

	use_FFT = true;

	k = next_power_of_two(n);
	l = next_power_of_two(2*n - 3);

	lidia_size_t l_safe = next_power_of_two(2*n - 1);
	//in functions multiply,square it must be guaranteed that both fft_reps
	//use an fft_table of at least this size

	if (re_init) {
		fft_data fd2(p, l_safe);
		fd = fd2;
	}
	else
		fd.init(p, l_safe);

	FRep.init(fd);
	FRep.set_size(k);
	FRep.to_fft_rep(pol);

	Fp_polynomial P1, P2;

	copy_reverse(P1, pol, 0, n);
	invert(P2, P1, n-1);
	copy_reverse(P1, P2, 0, n-2);

	HRep.init(fd);
	HRep.set_size(l);
	HRep.to_fft_rep(P1);
}



Fp_poly_modulus & Fp_poly_modulus::operator = (const Fp_poly_modulus & P)
{
	debug_handler("Fp_poly_modulus", "operator = (Fp_poly_modulus&)");
	build(P.modulus());
	return *this;
}



#if 0
static
Fp_polynomial f_star(const Fp_polynomial &f)
{
	Fp_polynomial P1, P2, a;
	copy_reverse(P1, f, 0, f.degree());
	invert(P2, P1, f.degree());
	copy_reverse(a, P2, 0, f.degree()-2);
	return a;
}
#endif



// x = a % f
// deg(a) <= 2(n-1), where n = f.degree()
void Fp_poly_modulus::rem21(Fp_polynomial& x, const Fp_polynomial& a) const
{
	debug_handler("Fp_poly_modulus", "rem21(Fp_polynomial&, Fp_polynomial&)");
	lidia_size_t i, index, K, da, ds;

	da = a.degree();

	if (da > 2*n-2) {
		lidia_error_handler("Fp_poly_modulus", "rem21(Fp_polynomial&, "
				    "Fp_polynomial&)::bad_args");
		return;
	}

	if (da < n) {
		x.assign(a);
		return;
	}

	if (!use_FFT || (da - n) <= crov) {
		plain_rem(x, a, f);
		return;
	}

	Fp_polynomial P1;
	P1.set_max_degree(n-1);

	lidia_size_t num_primes = fd.number_of_primes();
	if (num_primes == 1) {
		fft_rep Ra(fd);
		Ra.set_size(l);
		Ra.to_fft_rep(a, n, 2*(n-1));
		multiply(Ra, HRep, Ra);
		Ra.from_fft_rep(P1, n-2, 2*n-4);

		Ra.set_size(k);
		Ra.to_fft_rep(P1);
		multiply(Ra, FRep, Ra);
		Ra.from_fft_rep(P1, 0, n-1);
	}
	else {
		modular_fft_rep R1(fd); // more space efficient version

		R1.set_size(l);
		for (index = 0; index < num_primes; index++) {
			R1.to_modular_fft_rep(a, n, 2*(n-1), index);
			multiply(R1, HRep, R1, index);
			R1.from_modular_fft_rep(n-2, 2*n-4, index);
		}
		R1.get_result(P1, n-2, 2*n-4);

		R1.set_size(k);
		for (index = 0; index < num_primes; index++) {
			R1.to_modular_fft_rep(P1, index);
			multiply(R1, FRep, R1, index);
			R1.from_modular_fft_rep(0, n-1, index);
		}
		R1.get_result(P1, 0, n-1);
	}


	ds = P1.degree();
	K = 1 << k;

	x.set_modulus(P1);
	x.set_degree(n-1);
	const bigint* aa = a.coeff;
	const bigint* ss = P1.coeff;
	bigint* xx = x.coeff;
	const bigint & p = P1.modulus();

	for (i = 0; i < n; i++) {
		if (i <= ds)
			SubMod(xx[i], aa[i], ss[i], p);
		else
			xx[i].assign(aa[i]);

		if (i + K <= da)
			AddMod(xx[i], xx[i], aa[i+K], p);
	}

	x.remove_leading_zeros();
}



// x = a % f, no restrictions on deg(a);  makes repeated calls to rem21
void remainder(Fp_polynomial& x, const Fp_polynomial& a, const Fp_poly_modulus& F)
{
	debug_handler("Fp_poly_modulus", "remainder(Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	a.comp_modulus(F.f, "Fp_poly_modulus::remainder");

	lidia_size_t da = a.degree();
	lidia_size_t n = F.n;

	if (da <= 2*n-2) {
		F.rem21(x, a);
		return;
	}
	else
		if (!F.use_FFT) {
			plain_rem(x, a, F.f);
			return;
		}

	Fp_polynomial buf;
	buf.set_modulus(F.f);
	buf.set_max_degree(2*n-2);

	lidia_size_t a_len = da+1;

	while (a_len > 0) {
		lidia_size_t old_buf_len = buf.degree() + 1;
		lidia_size_t amt = comparator< lidia_size_t >::min(2*n - 1 - old_buf_len, a_len);

		buf.set_degree(old_buf_len+amt-1);

		lidia_size_t i;

		for (i = old_buf_len+amt-1; i >= amt; i--)
			buf[i].assign(buf[i-amt]);

		for (i = amt-1; i >= 0; i--)
			buf[i].assign(a[a_len-amt+i]);
		buf.remove_leading_zeros();

		F.rem21(buf, buf);
		a_len -= amt;
	}
	x.assign(buf);
}



// x = (a * b) % f
// deg(a), deg(b) < deg_f
void multiply(Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& b,
	      const Fp_poly_modulus& F)
{
	debug_handler("Fp_poly_modulus", "multiply(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	a.comp_modulus(b, "Fp_poly_modulus::multiply");
	a.comp_modulus(F.modulus(), "Fp_poly_modulus::multiply");

	lidia_size_t  da, db, d, n, k, index;

	da = a.degree();
	db = b.degree();
	n = F.n;

	if (da >= n || db >= n) {
		lidia_error_handler("Fp_poly_modulus", "multiply(Fp_polynomial&, "
				    "Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)::degree of "
				    "Fp_polynomials must be < degree of Fp_poly_modulus");
		return;
	}

	if (!F.use_FFT || da + db - n <= F.crov) {
		Fp_polynomial P1;
		plain_mul(P1, a, b);
		plain_rem(x, P1, F.f);
		return;
	}

	d = da + db + 1;
	k = next_power_of_two(d);
	k = comparator< lidia_size_t >::max(k, F.k);
	//here, k <= next_power_of_two( 2*n - 1 )


	Fp_polynomial P1;
	P1.set_max_degree(n - 1);

	lidia_size_t num_primes = F.fd.number_of_primes();
	if (num_primes == 1)	//TOM: probably isn't worth the distinction...
	{
		fft_rep Ra, Rb;
		Ra.init(F.fd);
		Rb.init(F.fd);
		Ra.set_size(k);
		Rb.set_size(k);

		Ra.to_fft_rep(a);
		Rb.to_fft_rep(b);
		multiply(Ra, Ra, Rb);
		Ra.from_fft_rep(P1, n, d-1);

		Rb.set_size(F.k);
		reduce(Rb, Ra, F.k);

		Ra.set_size(F.l);
		Ra.to_fft_rep(P1);
		multiply(Ra, F.HRep, Ra);
		Ra.from_fft_rep(P1, n-2, 2*n-4);

		Ra.set_size(F.k);
		Ra.to_fft_rep(P1);
		multiply(Ra, F.FRep, Ra);

		subtract(Rb, Rb, Ra);
		Rb.from_fft_rep(x, 0, n-1);
	}
	else {
		fft_rep R1;
		R1.init(F.fd);
		R1.set_size(F.k);
		modular_fft_rep Ra(F.fd); //allocate enough space
		//for convol. of deg. k, F.k, F.l
		modular_fft_rep Rb;
		Rb.init(F.fd);
		Ra.set_size(k);
		Rb.set_size(k);
		for (index = 0; index < num_primes; index++) {
			Ra.to_modular_fft_rep(a, index);
			Rb.to_modular_fft_rep(b, index);
			multiply(Ra, Ra, Rb, index);
			reduce(R1, Ra, F.k, index); // store results
			Ra.from_modular_fft_rep(n, d-1, index);
		}
		Ra.get_result(P1, n, d-1);

		Ra.set_size(F.l);
		for (index = 0; index < num_primes; index++) {
			Ra.to_modular_fft_rep(P1, index);
			multiply(Ra, F.HRep, Ra, index);
			Ra.from_modular_fft_rep(n-2, 2*n-4, index);
		}
		Ra.get_result(P1, n-2, 2*n-4);

		Ra.set_size(F.k);
		Rb.set_size(F.k);
		for (index = 0; index < num_primes; index++) {
			Rb.to_modular_fft_rep(P1, index);
			multiply(Rb, F.FRep, Rb, index);

			subtract(Ra, R1, Rb, index);
			Ra.from_modular_fft_rep(0, n-1, index);
		}
		Ra.get_result(x, 0, n-1);
	}
}



// x = a^2 % f			a.degree() < f.degree()
void square(Fp_polynomial& x, const Fp_polynomial& a, const Fp_poly_modulus& F)
{
	debug_handler("Fp_poly_modulus", "square(Fp_polynomial&, Fp_polynomial&, Fp_poly_modulus&)");

	a.comp_modulus(F.modulus(), "Fp_poly_modulus::square");

	lidia_size_t  da, d, n, k, index;

	da = a.degree();
	n = F.n;

	if (da >= n) {
		lidia_error_handler("Fp_poly_modulus", "square(Fp_polynomial&, "
				    "Fp_polynomial&, Fp_poly_modulus&)::degree of Fp_polynomial "
				    "must be < degree of Fp_poly_modulus");
		return;
	}

	if (!F.use_FFT || 2*da - n <= F.crov) {
		Fp_polynomial P1;
		plain_sqr(P1, a);
		plain_rem(x, P1, F.f);
		return;
	}


	d = 2*da + 1;

	k = next_power_of_two(d);
	k = comparator< lidia_size_t >::max(k, F.k);
	//here, k <= next_power_of_two( 2*n - 1 )

	Fp_polynomial P1;
	P1.set_max_degree(n - 1);

	modular_fft_rep Ra(F.fd);
	modular_fft_rep Rb;

	fft_rep R1;
	R1.init(F.fd);
	R1.set_size(F.k);
	Ra.set_size(k);
	for (index = 0; index < F.fd.number_of_primes(); index++) {
		Ra.to_modular_fft_rep(a, index);
		multiply(Ra, Ra, Ra, index);
		reduce(R1, Ra, F.k, index); //save results fo future use
		Ra.from_modular_fft_rep(n, d-1, index);
	}
	Ra.get_result(P1, n, d-1);

	Ra.set_size(F.l);
	for (index = 0; index < F.fd.number_of_primes(); index++) {
		Ra.to_modular_fft_rep(P1, index);
		multiply(Ra, F.HRep, Ra, index);
		Ra.from_modular_fft_rep(n-2, 2*n-4, index);
	}
	Ra.get_result(P1, n-2, 2*n-4);

	Ra.set_size(F.k);
	Rb.init(F.fd);
	Rb.set_size(F.k);

	for (index = 0; index < F.fd.number_of_primes(); index++) {
		Rb.to_modular_fft_rep(P1, index);
		multiply(Rb, F.FRep, Rb, index);

		subtract(Ra, R1, Rb, index);
		Ra.from_modular_fft_rep(0, n-1, index);
	}
	Ra.get_result(x, 0, n-1);
}



// x = a^e % f
void power(Fp_polynomial& h, const Fp_polynomial& g, const bigint& e,
	   const Fp_poly_modulus& F)
{
	debug_handler("Fp_poly_modulus", "power(Fp_polynomial&, Fp_polynomial&, bigint&, Fp_poly_modulus&)");

	const Fp_polynomial &f = F.modulus();
	g.comp_modulus(f, "Fp_poly_modulus::power");

	if (e.is_negative()) {
		lidia_error_handler("Fp_poly_modulus", "power(Fp_polynomial&, "
				    "Fp_polynomial&, bigint&, Fp_poly_modulus&)::exponent must "
				    "be positive");
		return;
	}

	lidia_size_t i, n = e.bit_length();

	Fp_poly_multiplier G(g, F);

	h.set_modulus(f);
	h.set_max_degree(f.degree() - 1);
	h.assign_one();

	for (i = n - 1; i >= 0; i--) {
		square(h, h, F);
		if (e.bit(i))
			multiply(h, h, G, F); //Fp_poly_multiplier
	}
}



// x = X^e % f
void power_x(Fp_polynomial& h, const bigint& e, const Fp_poly_modulus& F)
{
	debug_handler("Fp_poly_modulus", "power_x(Fp_polynomial&, bigint&, Fp_poly_modulus&)");

	if (e.is_negative()) {
		lidia_error_handler("Fp_poly_modulus", "power_x(Fp_polynomial&, "
				    "bigint&, Fp_poly_modulus&)::exponent must be positive");
		return;
	}

	lidia_size_t i, n = e.bit_length();

	const Fp_polynomial &f = F.modulus();
	h.set_modulus(f);
	h.set_max_degree(f.degree() - 1);
	h.assign_one();

	for (i = n - 1; i >= 0; i--) {
		square(h, h, F);
		if (e.bit(i))
			multiply_by_x_mod(h, h, f);
	}
}



// x = (X + a)^e % f
void power_x_plus_a(Fp_polynomial& h, const bigint& a, const bigint& e,
		    const Fp_poly_modulus& F)
{
	debug_handler("Fp_poly_modulus", "power_x_plus_a(Fp_polynomial&, bigint&, bigint&, Fp_poly_modulus&)");

	if (e.is_negative()) {
		lidia_error_handler("Fp_poly_modulus", "power_x_plus_a(Fp_polynomial&, "
				    "bigint&, bigint&, Fp_poly_modulus&)::exponent must be > 0");
		return;
	}

	const Fp_polynomial &f = F.modulus();
	lidia_size_t n = f.degree();

	Fp_polynomial t1, t2;
	t1.set_max_degree(n-1);
	t2.set_max_degree(n-1);

	bigint la;
	Remainder(la, a, f.modulus()); //allows input to alias output

	lidia_size_t i, m = e.bit_length();

	h.set_modulus(f);
	h.set_max_degree(n-1);
	h.assign_one();

	for (i = m - 1; i >= 0; i--) {
		square(h, h, F);
		if (e.bit(i)) {
			multiply_by_x_mod(t1, h, f);
			multiply_by_scalar(t2, h, la);
			add(h, t1, t2);
		}
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

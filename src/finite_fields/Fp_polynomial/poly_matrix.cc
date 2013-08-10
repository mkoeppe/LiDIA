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
#include	"LiDIA/finite_fields/Fp_polynomial_fft.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***************************************************************
//                                    auxiliary functions
//	all assume that all input parameters have correct modulus
//***************************************************************

#if 0
//TPf: speedups in plain_div_rem make plain_div_rem_vec obsolete...

void plain_div_rem_vec(Fp_polynomial& q, Fp_polynomial& r, const Fp_polynomial& a, const Fp_polynomial& b, bigint *x)
	//special version of plain_div_rem, uses memory offered by vector x
	//about 10% faster
	//assumes x.size() > a.degree()+1
	//cmp. plain_rem_vec in file gcd.cc
{
	debug_handler("Fp_polynomial", "plain_div_rem_vec (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, bigint *)");

	lidia_size_t da, db, dq, i, j;
	const bigint *bp;
	const bigint &p = a.modulus();

	bigint LCInv, t, s;

	da = a.degree();
	db = b.degree();
	dq = da - db;

	if (db < 0) {
		lidia_error_handler("Fp_polynomial", "plain_div_rem_vec(...)::"
				    "division by zero");
		return;
	}

	if (dq < 0) {
		r.assign(a);
		q.set_modulus(a); //assigns zero
		return;
	}

	Fp_polynomial lb;
	if (&q == &b) {
		lb.assign(b);
		bp = lb.coeff;
	}
	else
		bp = b.coeff;

	bool lc_is_one;
	if (bp[db].is_one())
		lc_is_one = true;
	else {
		lc_is_one = false;
		InvMod(LCInv, bp[db], p);
	}


	for (i = 0; i <= da; i++)
		x[i].assign(a.coeff[i]);

	q.set_modulus(a);
	q.set_degree(dq);
	bigint *qp = q.coeff;

	for (i = dq; i >= 0; i--) {
		Remainder(t, x[i+db], p);
		if (!lc_is_one)
			MulMod(t, t, LCInv, p);
		qp[i].assign(t);

		negate(t, t);

		for (j = db-1; j >= 0; j--) {
			multiply(s, t, bp[j]);
			add(x[i+j], x[i+j], s);
		}
	}

	r.set_modulus(a);
	r.set_degree(db-1);
	for (i = 0; i < db; i++)
		Remainder(r.coeff[i], x[i], p);
	r.remove_leading_zeros();
}
#endif



//***************************************************************
//
//			  class poly_matrix
//
//****************************************************************

poly_matrix::poly_matrix(const bigint& p)
{
	debug_handler("poly_matrix", "poly_matrix(bigint&)");
	elts[0][0].set_modulus(p);
	elts[0][1].set_modulus(elts[0][0]);
	elts[1][0].set_modulus(elts[0][0]);
	elts[1][1].set_modulus(elts[0][0]);
}



poly_matrix & poly_matrix::operator = (const poly_matrix& M)
{
	debug_handler("poly_matrix", "operator = (poly_matrix&)");
	elts[0][0].assign(M.elts[0][0]);
	elts[0][1].assign(M.elts[0][1]);
	elts[1][0].assign(M.elts[1][0]);
	elts[1][1].assign(M.elts[1][1]);

	return *this;
}



void poly_matrix::kill()
{
	debug_handler("poly_matrix", "kill (void)");
	elts[0][0].kill();
	elts[0][1].kill();
	elts[1][0].kill();
	elts[1][1].kill();
}



void poly_matrix::multiply_and_kill(poly_matrix& B, poly_matrix& C)
	// (*this) = B*C, B and C are destroyed
{
	debug_handler("poly_matrix", "multiply_and_kill (poly_matrix&, poly_matrix&)");
	lidia_size_t db = (B(1, 1)).degree();
	lidia_size_t dc = (C(1, 1)).degree();
	lidia_size_t da = db + dc;

	lidia_size_t k = next_power_of_two(da+1);
	const bigint& p = B(0, 0).modulus();
	fft_data F(p, k);
	modular_fft_rep B00(F), B01(F), B10(F), B11(F),
		C0(F), C1(F),
		A00(F), A10(F);

	lidia_size_t index;
	for (index = 0; index < F.number_of_primes(); index++) {
		B00.to_modular_fft_rep(B(0, 0), index);
		B01.to_modular_fft_rep(B(0, 1), index);
		B10.to_modular_fft_rep(B(1, 0), index);
		B11.to_modular_fft_rep(B(1, 1), index);

		C0.to_modular_fft_rep(C(0, 0), index);
		C1.to_modular_fft_rep(C(1, 0), index);

		add_mul(A00, B00, C0, B01, C1, index);
		A00.from_modular_fft_rep(0, da, index);

		add_mul(A10, B10, C0, B11, C1, index);
		A10.from_modular_fft_rep(0, da, index);

		C0.to_modular_fft_rep(C(0, 1), index);
		C1.to_modular_fft_rep(C(1, 1), index);

		add_mul(B01, B00, C0, B01, C1, index);
		B01.from_modular_fft_rep(0, da, index);

		add_mul(B11, B10, C0, B11, C1, index);
		B11.from_modular_fft_rep(0, da, index);
	}
#if 0
	Fp_polynomial tmp1, tmp2, tmp3, tmp4;
	LiDIA::multiply(tmp1, B(0, 0), C(0, 0));
	LiDIA::multiply(tmp2, B(0, 1), C(1, 0));
	LiDIA::multiply(tmp3, B(1, 0), C(0, 0));
	LiDIA::multiply(tmp4, B(1, 1), C(1, 0));
	add(elts[0][0], tmp1, tmp2);
	add(elts[1][0], tmp3, tmp4);
	LiDIA::multiply(tmp1, B(0, 0), C(0, 1));
	LiDIA::multiply(tmp2, B(0, 1), C(1, 1));
	LiDIA::multiply(tmp3, B(1, 0), C(0, 1));
	LiDIA::multiply(tmp4, B(1, 1), C(1, 1));
	add(elts[0][1], tmp1, tmp2);
	add(elts[1][1], tmp3, tmp4);
#endif

	B.kill();
	C.kill();

	A00.get_result(elts[0][0], 0, da);
	A10.get_result(elts[1][0], 0, da);
	B01.get_result(elts[0][1], 0, da);
	B11.get_result(elts[1][1], 0, da);
}



void poly_matrix::multiply(Fp_polynomial& U, Fp_polynomial& V) const
	// (U, V)^T = M*(U, V)^T
{
	debug_handler("poly_matrix", "multiply (Fp_polynomial&, Fp_polynomial&)");

	const bigint& p = elts[0][0].modulus();

	if (p.is_zero()) {
		lidia_error_handler("poly_matrix", "multiply(Fp_polynomial&, "
				    "Fp_polynomial&)::modulus was not set");
		return;
	}

	lidia_size_t d = U.degree() - (elts[1][1]).degree();
	lidia_size_t k = next_power_of_two(d - 1);

	// When the GCD algorithm is run on polynomials of degree n, n-1,
	// where n is a power of two, then d-1 is likely to be a power of two.
	// It would be more natural to set k = next_power_of_two(d+1), but this
	// would be much less efficient in this case.

	lidia_size_t n = (1 << k);
	lidia_size_t xx;
	//bigint a0, a1, b0, b1, c0, d0, u0, u1, v0, v1, nu0, nu1, nv0; //static bigmod
	bigint nu0, nu1, nv0; //static bigmod
	bigint t1, t2; //static bigint

	if (n == d-1)
		xx = 1;
	else
		if (n == d)
			xx = 2;
		else
			xx = 3;

	switch (xx) {
	case 1:
#if 0
		  (elts[0][0]).get_coefficient(a0, 0);
		  (elts[0][0]).get_coefficient(a1, 1);
		  (elts[0][1]).get_coefficient(b0, 0);
		  (elts[0][1]).get_coefficient(b1, 1);
		  (elts[1][0]).get_coefficient(c0, 0);
		  (elts[1][1]).get_coefficient(d0, 0);
		  U.get_coefficient(u0, 0);
		  U.get_coefficient(u1, 1);
		  V.get_coefficient(v0, 0);
		  V.get_coefficient(v1, 1);

		  LiDIA::multiply(t1, a0, u0);
		  LiDIA::multiply(t2, b0, v0);
		  add(t1, t1, t2);
		  Remainder(nu0, t1, p);

		  LiDIA::multiply(t1, a1, u0);
		  LiDIA::multiply(t2, a0, u1);
		  add(t1, t1, t2);
		  LiDIA::multiply(t2, b1, v0);
		  add(t1, t1, t2);
		  LiDIA::multiply(t2, b0, v1);
		  add(t1, t1, t2);
		  Remainder(nu1, t1, p);

		  LiDIA::multiply(t1, c0, u0);
		  LiDIA::multiply(t2, d0, v0);
		  add (t1, t1, t2);
		  Remainder(nv0, t1, p);
#endif
		  LiDIA::multiply(t1, (elts[0][0])[0], U[0]);
		  LiDIA::multiply(t2, (elts[0][1])[0], V[0]);
		  add(t1, t1, t2);
		  Remainder(nu0, t1, p);

		  LiDIA::multiply(t1, (elts[0][0])[1], U[0]);
		  LiDIA::multiply(t2, (elts[0][0])[0], U[1]);
		  add(t1, t1, t2);
		  LiDIA::multiply(t2, (elts[0][1])[1], V[0]);
		  add(t1, t1, t2);
		  LiDIA::multiply(t2, (elts[0][1])[0], V[1]);
		  add(t1, t1, t2);
		  Remainder(nu1, t1, p);

		  LiDIA::multiply(t1, (elts[1][0])[0], U[0]);
		  LiDIA::multiply(t2, (elts[1][1])[0], V[0]);
		  add (t1, t1, t2);
		  Remainder(nv0, t1, p);

		  break;

	case 2:

#if 0
		  (elts[0][0]).get_coefficient(a0, 0);
		  (elts[0][1]).get_coefficient(b0, 0);
		  U.get_coefficient(u0, 0);
		  V.get_coefficient(v0, 0);

		  LiDIA::multiply(t1, a0, u0);
		  LiDIA::multiply(t2, b0, v0);
		  add(t1, t1, t2);
		  Remainder(nu0, t1, p);
#endif
		  LiDIA::multiply(t1, (elts[0][0])[0], U[0]);
		  LiDIA::multiply(t2, (elts[0][1])[0], V[0]);
		  add(t1, t1, t2);
		  Remainder(nu0, t1, p);

		  break;

	case 3:
		break;
	}

	fft_data F(p, k);
	modular_fft_rep RU(F), RV(F), R1(F), R2(F), R3(F), R4(F), T1(F), T2(F);

	lidia_size_t index;
	for (index = 0; index < F.number_of_primes(); index++) {
		RU.to_modular_fft_rep(U, 0, U.degree(), index);
		RV.to_modular_fft_rep(V, 0, U.degree(), index);

		R1.to_modular_fft_rep(elts[0][0], index);
		R2.to_modular_fft_rep(elts[0][1], index);

		add_mul(R1, R1, RU, R2, RV, index);

		R1.from_modular_fft_rep(0, d, index);

		R3.to_modular_fft_rep(elts[1][0], index);
		R4.to_modular_fft_rep(elts[1][1], index);

		add_mul(R3, R3, RU, R4, RV, index);
		R3.from_modular_fft_rep(0, d-1, index);
	}
	R1.get_result(U, 0, d);
	R3.get_result(V, 0, d-1);


#if 0
	Fp_polynomial tmp1, tmp2, tmp3, tmp4;
	LiDIA::multiply(tmp1, U, elts[0][0]);
	LiDIA::multiply(tmp2, V, elts[0][1]);

	LiDIA::multiply(tmp3, U, elts[1][0]);
	LiDIA::multiply(tmp4, V, elts[1][1]);

	add(U, tmp1, tmp2);
	add(V, tmp3, tmp4);

#endif


	// now fix-up results

	switch (xx) {
	case 1:
		SubMod(t1, U[0], nu0, p);
		U.set_coefficient(t1, d-1);
		U.set_coefficient(nu0, 0);

		SubMod(t1, U[1], nu1, p);
		U.set_coefficient(t1, d);
		U.set_coefficient(nu1, 1);

		SubMod(t1, V[0], nv0, p);
		V.set_coefficient(t1, d-1);
		V.set_coefficient(nv0, 0);

		break;

	case 2:
		SubMod(t1, U[0], nu0, p);
		U.set_coefficient(t1, d);
		U.set_coefficient(nu0, 0);

		break;
	}
}



// deg(U) > deg(V),   1 <= d_red <= deg(U)+1.
// This computes a 2 x 2 polynomial matrix M_out such that
//    M_out * (U, V)^T = (U', V')^T,
// where U', V' are consecutive polynomials in the Euclidean remainder
// sequence of U, V, and V' is the polynomial of highest degree
// satisfying deg(V') <= deg(U) - d_red.
void poly_matrix::half_gcd(const Fp_polynomial& U, const Fp_polynomial& V, lidia_size_t d_red)
{
	debug_handler("poly_matrix", "half_gcd (Fp_polynomial&, Fp_polynomial&, lidia_size_t)");

	if (V.is_zero() || V.degree() <= U.degree() - d_red) {
		(elts[0][0]).assign_one(); (elts[0][1]).assign_zero();
		(elts[1][0]).assign_zero(); (elts[1][1]).assign_one();
		return;
	}

	lidia_size_t n = U.degree() - 2*d_red + 2;
	if (n < 0)
		n = 0;


	Fp_polynomial U1, V1;
	shift_right(U1, U, n);
	shift_right(V1, V, n);

	const bigint &p = U.modulus();
	if (d_red < Fp_polynomial::crossovers.halfgcd_crossover(p)) {
		iter_half_gcd(U1, V1, d_red);
		return;
	}

	lidia_size_t d1 = (d_red + 1)/2;
	if (d1 < 1)
		d1 = 1;
	if (d1 >= d_red)
		d1 = d_red - 1;
	poly_matrix M1(p);


	M1.half_gcd(U1, V1, d1);
	M1.multiply(U1, V1);
	lidia_size_t d2 = V1.degree() - U.degree() + n + d_red;

	if (V1.is_zero() || d2 <= 0) {
		(*this) = M1;
		return;
	}


	Fp_polynomial Q;
	poly_matrix M2(p);

	div_rem(Q, U1, U1, V1);
	swap(U1, V1);

	M2.half_gcd(U1, V1, d2);

	Fp_polynomial t;
	t.set_max_degree((M1(1, 1)).degree()+Q.degree());


	LiDIA::multiply(t, Q, M1(1, 0));
	subtract(t, M1(0, 0), t);
	swap(M1(0, 0), M1(1, 0));
	swap(M1(1, 0), t);

	t.assign_zero();
	t.set_max_degree((M1(1, 1)).degree()+Q.degree());

	LiDIA::multiply(t, Q, M1(1, 1));
	subtract(t, M1(0, 1), t);
	swap(M1(0, 1), M1(1, 1));
	swap(M1(1, 1), t);

	multiply_and_kill(M2, M1);
}



// same as above, except that U is replaced by U', and V by V'
void poly_matrix::xhalf_gcd(Fp_polynomial& U, Fp_polynomial& V, lidia_size_t d_red)
{
	debug_handler("poly_matrix", "xhalf_gcd (Fp_polynomial&, Fp_polynomial&, lidia_size_t)");
	if (V.is_zero() || V.degree() <= U.degree() - d_red) {
		(elts[0][0]).assign_one(); (elts[0][1]).assign_zero();
		(elts[1][0]).assign_zero(); (elts[1][1]).assign_one();
		return;
	}

	lidia_size_t du = U.degree();

	const bigint &p = U.modulus();
	if (d_red < Fp_polynomial::crossovers.halfgcd_crossover(p)) {
		iter_half_gcd(U, V, d_red);
		return;
	}

	lidia_size_t d1 = (d_red + 1)/2;
	if (d1 < 1)
		d1 = 1;
	if (d1 >= d_red)
		d1 = d_red - 1;

	poly_matrix M1(p);

	M1.half_gcd(U, V, d1);
	M1.multiply(U, V);

	lidia_size_t d2 = V.degree() - du + d_red;

	if (V.is_zero() || d2 <= 0) {
		(*this) = M1;
		return;
	}

	Fp_polynomial Q;
	poly_matrix M2(p);

	div_rem(Q, U, U, V);
	swap(U, V);

	M2.xhalf_gcd(U, V, d2);

	Fp_polynomial t;
	t.set_max_degree((M1(1, 1)).degree()+Q.degree());

	LiDIA::multiply(t, Q, M1(1, 0));
	subtract(t, M1(0, 0), t);
	swap(M1(0, 0), M1(1, 0));
	swap(M1(1, 0), t);

	t.assign_zero();
	t.set_max_degree((M1(1, 1)).degree()+Q.degree());

	LiDIA::multiply(t, Q, M1(1, 1));
	subtract(t, M1(0, 1), t);
	swap(M1(0, 1), M1(1, 1));
	swap(M1(1, 1), t);

	multiply_and_kill(M2, M1);
}



void poly_matrix::iter_half_gcd(Fp_polynomial& U, Fp_polynomial& V, lidia_size_t d_red)
{
	debug_handler("poly_matrix", "iter_half_gcd (Fp_polynomial&, Fp_polynomial&, lidia_size_t)");

	elts[0][0].set_max_degree(d_red-1);
	elts[0][1].set_max_degree(d_red-1);
	elts[1][0].set_max_degree(d_red-1);
	elts[1][1].set_max_degree(d_red-1);

	(elts[0][0]).assign_one(); (elts[0][1]).assign_zero();
	(elts[1][0]).assign_zero(); (elts[1][1]).assign_one();

	lidia_size_t goal = U.degree() - d_red;

	if (V.degree() <= goal)
		return;

	Fp_polynomial Q, t;
	t.set_max_degree(d_red - 1);

	while (V.degree() > goal) {
		plain_div_rem(Q, U, U, V);
		// TPf: original code: plain_div_rem_vec(Q, U, U, V, TEMP);
		//      (speedups in plain_div_rem make plain_div_rem_vec obsolete)
		swap(U, V);

		LiDIA::multiply(t, Q, elts[1][0]);
		subtract(t, elts[0][0], t);
		elts[0][0].assign(elts[1][0]);
		elts[1][0].assign(t);

		LiDIA::multiply(t, Q, elts[1][1]);
		subtract(t, elts[0][1], t);
		elts[0][1].assign(elts[1][1]);
		elts[1][1].assign(t);
	}

}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

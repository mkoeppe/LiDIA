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
//                mainly based on V. Shoup's Fp_poly_modulus for Fp_polynomial
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf_polynomial.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//******************************************************************
//
//			  class gf_poly_modulus
//
//******************************************************************

// If you need to do a lot of arithmetic modulo a fixed f,
// build gf_poly_modulus F for f.  This pre-computes information about f
// that speeds up the computation a great deal.
// f should be monic, and deg(f) > 0.
// For mathematical background, see V. Shoup's paper


gf_poly_modulus::gf_poly_modulus(const gf_poly_modulus &p) :
	F_plain(p.F_plain), F(p.F), F_recip(p.F_recip)
{ }



gf_poly_modulus::gf_poly_modulus(const gf_polynomial &f)
{
	build(f);
}



void gf_poly_modulus::build(const gf_polynomial &f)
{
	debug_handler("gf_poly_modulus", "build(gf_polynomial&)");

	if (f.degree() <= 0 || !f.is_monic())
		lidia_error_handler("gf_poly_modulus", "build(gf_polynomial&)::"
				    "argument must be monic and of degree > 0");

	F_plain.assign(f);

	int fd = f.get_field().degree();
	if (f.degree() == 1 || (f.get_field().characteristic() == 2 &&
				(fd > 8 || (fd > 1 && f.degree() < 100*fd)))) {
		use_Fp = false;
		return;
	}

	use_Fp = true;

	gf_polynomial P1, P2;
	copy_reverse(P1, f, 0, f.degree());
	invert(P2, P1, f.degree());
	copy_reverse(P1, P2, 0, f.degree()-2);

	to_Kronecker(F_recip, P1, 0, P1.degree());
	to_Kronecker(F, f, 0, f.degree());
}



// x = a % f
// deg(a) <= 2(n-1), where n = f.degree()
void gf_poly_modulus::rem21(gf_polynomial& x, const gf_polynomial& a) const
{
	debug_handler("gf_poly_modulus", "rem21(gf_polynomial&, gf_polynomial&)");
	lidia_size_t da, n;
	da = a.degree();
	n = F_plain.degree();

	if (da > 2*n-2)
		lidia_error_handler("gf_poly_modulus", "rem21(gf_polynomial&, gf_polynomial&)::bad_args");

	if (da < n) {
		x.assign(a);
		return;
	}
	if (!use_Fp) {
		remainder(x, a, F_plain);
		return;
	}

	Fp_polynomial r1, r2;
	gf_polynomial P1;
	P1.ffield = a.ffield;
	gf_polynomial::build_frame(a.ffield);

	to_Kronecker(r1, a, n, 2*(n-1));
	multiply(r1, F_recip, r1);
	from_Kronecker(P1, r1, n-2, 2*n-4);

	to_Kronecker(r1, P1, 0, P1.degree());
	multiply(r1, F, r1);
	from_Kronecker(P1, r1, 0, n-1);

	lidia_size_t i, ds = P1.degree();
	x.ffield = a.ffield;
	x.set_degree(n-1);
	for (i = 0; i < n; i++) {
		if (i <= ds) subtract(x[i], a[i], P1[i]);
		else         x[i].assign(a[i]);
	}
	x.remove_leading_zeros();
	gf_polynomial::delete_frame();
}



// x = a % f, no restrictions on deg(a);  makes repeated calls to rem21
void remainder(gf_polynomial& x, const gf_polynomial& a, const gf_poly_modulus& F)
{
	debug_handler("gf_poly_modulus", "remainder(gf_polynomial&, gf_polynomial&, gf_poly_modulus&)");

	lidia_size_t da = a.degree();
	lidia_size_t n = F.modulus().degree();

	if (da <= 2*n-2) {
		F.rem21(x, a);
		return;
	}

	if (!F.use_Fp) {
		remainder(x, a, F.F_plain);
		return;
	}

	gf_polynomial buf;
	buf.ffield = (gf_polynomial::common_field(a.ffield, F.F_plain.ffield));
	gf_polynomial::build_frame(buf.ffield);

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
	gf_polynomial::delete_frame();
	x.assign(buf);
}



// x = (a * b) % f
// deg(a), deg(b) < deg_f
void multiply(gf_polynomial& x, const gf_polynomial& a, const gf_polynomial& b, const gf_poly_modulus& F)
{
	debug_handler("gf_poly_modulus", "multiply(gf_polynomial&, gf_polynomial&, gf_polynomial&, gf_poly_modulus&)");

	lidia_size_t  da, db, n;

	da = a.degree();
	db = b.degree();
	n = F.modulus().degree();

	if (da >= n || db >= n)
		lidia_error_handler("gf_poly_modulus", "multiply(gf_polynomial&, gf_polynomial&, gf_polynomial&, gf_poly_modulus&)::degree of gf_polynomials must be < degree of gf_poly_modulus");

	if (!F.use_Fp) {
		remainder(x, a*b, F.F_plain);
		return;
	}

	x.ffield = (gf_polynomial::common_field(
		gf_polynomial::common_field(a.ffield, b.ffield),
		F.F_plain.ffield));
	gf_polynomial::build_frame(x.ffield);

	Fp_polynomial r1, r2, prod;
	gf_polynomial P1;
	P1.ffield = x.ffield;

#if 0
	//HERE IS THE ORIGINAL AND READABLE CODE ...
	to_Kronecker(r1, a, 0, a.degree());
	to_Kronecker(r2, b, 0, b.degree());
	multiply(prod, r1, r2);
	from_Kronecker(P1, prod, n, da + db);

	to_Kronecker(r1, P1, 0, P1.degree());
	multiply(r2, F.F_recip, r1);
	from_Kronecker(P1, r2, n-2, 2*n-4);

	to_Kronecker(r1, P1, 0, P1.degree());
	multiply(r2, r1, F.F);
	from_Kronecker(x, r2, 0, n-1);

	from_Kronecker(P1, prod, 0, n-1);
	subtract(x, P1, x);
	gf_polynomial::delete_frame();
#endif
	//NOW, THE OPTIMIZED VERSION (EQUIVALENT TO THE UPPER ONE)

	lidia_size_t N = a.get_field().degree();
	//cmp. Kronecker substitution:
	//one coeff. of gf_polynomial  equals  (2*N-1) coeffs. of Fp_polynomial

	to_Kronecker(r1, a, 0, a.degree());
	to_Kronecker(r2, b, 0, b.degree());
	multiply(prod, r1, r2);
	from_Kronecker(P1, prod, n, da + db);
	trunc(prod, prod, (2*N-1)*n);
	//deletes the coefficients from prod which we have already converted

	to_Kronecker(r1, P1, 0, P1.degree());
	multiply(r2, F.F_recip, r1);
	from_Kronecker(P1, r2, n-2, 2*n-4);

	to_Kronecker(r1, P1, 0, P1.degree());
	multiply(r2, r1, F.F);

	//from_Kronecker(x, r2, 0, n-1);
	//from_Kronecker(P1, prod, 0, n-1);
	//subtract(x, P1, x);

	trunc(r2, r2, n*(2*N-1));
	subtract(prod, prod, r2);
	from_Kronecker(x, prod, 0, n-1);

	gf_polynomial::delete_frame();
}



// x = a^2 % f			a.degree() < f.degree()
void square(gf_polynomial& x, const gf_polynomial& a, const gf_poly_modulus& F)
{
	debug_handler("gf_poly_modulus", "square(gf_polynomial&, gf_polynomial&, gf_poly_modulus&)");


	lidia_size_t  da, n;

	da = a.degree();
	n = F.F_plain.degree();

	if (da >= n)
		debug_handler("gf_poly_modulus", "square(gf_polynomial&, gf_polynomial&, gf_poly_modulus&)::degree of gf_polynomial must be < degree of gf_poly_modulus");

	if (!F.use_Fp) {
		square(x, a);
		remainder(x, x, F.F_plain);
		return;
	}

	//compare :
	//multiply(gf_polynomial&, gf_polynomial&, gf_polynomial&, gf_poly_modulus&)

	Fp_polynomial r1, r2, prod;
	gf_polynomial P1;
	x.ffield = (gf_polynomial::common_field(a.ffield, F.F_plain.ffield));
	P1.ffield = x.ffield;
	gf_polynomial::build_frame(x.ffield);

	lidia_size_t N = a.get_field().degree();
	//cmp. Kronecker substitution:
	//one coeff. of gf_polynomial  equals  (2*N-1) coeffs. of Fp_polynomial

	to_Kronecker(r1, a, 0, a.degree());
	square(prod, r1);
	from_Kronecker(P1, prod, n, 2*a.degree());
	trunc(prod, prod, (2*N-1)*n);
	//deletes the coefficients from prod which we have already converted

	to_Kronecker(r1, P1, 0, P1.degree());
	multiply(r2, F.F_recip, r1);
	from_Kronecker(P1, r2, n-2, 2*n-4);

	to_Kronecker(r1, P1, 0, P1.degree());
	multiply(r2, r1, F.F);

	//from_Kronecker(x, r2, 0, n-1);
	//from_Kronecker(P1, prod, 0, n-1);
	//subtract(x, P1, x);

	trunc(r2, r2, n*(2*N-1));
	subtract(prod, prod, r2);
	from_Kronecker(x, prod, 0, n-1);

	gf_polynomial::delete_frame();

}



// x = a^e % f
void power(gf_polynomial& h, const gf_polynomial& g, const bigint& e, const gf_poly_modulus& F)
{
	debug_handler("gf_poly_modulus", "power(gf_polynomial&, gf_polynomial&, bigint&, gf_poly_modulus&)");

	gf_polynomial lg(g);

	if (e.is_negative())
		lidia_error_handler("gf_poly_modulus", "power(gf_polynomial&, gf_polynomial&, bigint&, gf_poly_modulus&)::exponent must be positive");

	int i, n = e.bit_length();

//###    gf_poly_multiplier G(lg, F);

	h.assign_one(gf_polynomial::common_field(g.ffield, F.F_plain.ffield));

	for (i = n - 1; i >= 0; i--) {
		square(h, h, F);
		if (e.bit(i))
//###	multiply(h, h, G, F);	//gf_poly_multiplier
			multiply(h, h, lg, F);
	}

}



// x = X^e % f
void power_x(gf_polynomial& h, const bigint& e, const gf_poly_modulus& F)
{
	debug_handler("gf_poly_modulus", "power_x(gf_polynomial&, bigint&, gf_poly_modulus&)");

	if (e.is_negative())
		lidia_error_handler("gf_poly_modulus", "power_x(gf_polynomial&, bigint&, gf_poly_modulus&)::exponent must be positive");

	int i, n = e.bit_length();

	h.assign_one(F.F_plain.get_field());

	for (i = n - 1; i >= 0; i--) {
		square(h, h, F);
		if (e.bit(i))
			multiply_by_x_mod(h, h, F.modulus());
	}

}



// x = (X + a)^e % f
void power_x_plus_a(gf_polynomial& h, const gf_element& a, const bigint& e, const gf_poly_modulus& F)
{
	debug_handler("gf_poly_modulus", "power_x_plus_a(gf_polynomial&, gf_element&, bigint&, gf_poly_modulus&)");

	if (e.is_negative())
		lidia_error_handler("gf_poly_modulus", "power_x_plus_a(gf_polynomial&, gf_element&, bigint&, gf_poly_modulus&)::exponent must be positive");


	gf_polynomial t1, t2;
	gf_element la(a); //allows input to alias output

	int i, n = e.bit_length();

	h.assign_one(gf_polynomial::common_field(a.get_field(),
						 F.F_plain.get_field()));

	for (i = n - 1; i >= 0; i--) {
		square(h, h, F);
		if (e.bit(i)) {
			multiply_by_x_mod(t1, h, F.modulus());
			multiply(t2, h, la);
			add(h, t1, t2);
		}
	}

}



#if 0
//
// The idea of gf_poly_multiplier (the analogon to 'poly_multiplier' for
// Fp_polynomial) does not speed up the multiplication significantly.
// Time is spent (with or without this class) in three multiplications of
// Fp_polynomial's. Conversions (see Kronecker substitution) do hardly cost
// anything - a this is the only thing we could save.
// Bad luck.
//

//******************************************************************
//
//			  class gf_poly_multiplier
//
//******************************************************************

// If you need to compute a * b % f for a fixed b, but for many a's
// (for example, computing powers of b modulo f), it is
// more efficient to first build a gf_poly_multiplier B for b,
// and then use the routine below.
// For mathematical background, see V. Shoup's paper

gf_poly_multiplier::gf_poly_multiplier() : mod(0)
{ }



gf_poly_multiplier::gf_poly_multiplier(const gf_poly_multiplier &B)
{
	if (mod != 0)
		build(B.modulus(), *mod);
}



gf_poly_multiplier::gf_poly_multiplier(const gf_polynomial& b,
				       const gf_poly_modulus& F)
{
	build(b, F);
}



gf_poly_multiplier::~gf_poly_multiplier()
{ }



void
gf_poly_multiplier::build(const gf_polynomial& b, const gf_poly_modulus& F)
{
	mod = &F;
	remainder(B_plain, b, F);
	to_Kronecker(B, b, 0, b.degree()); // store b  

	lidia_size_t n = F.modulus().degree();
	Fp_polynomial g;
	gf_polynomial h;
	h.ffield = b.ffield;
	gf_polynomial::build_frame(b.ffield);
	multiply(g, B, F.F_recip);
	from_Kronecker(h, g, n-1, 2*n-3);
	h.delete_frame();
	to_Kronecker(B_div_F, h, 0, h.degree()); // store (b*F_recip)/X^n  
}



void multiply(gf_polynomial &x, const gf_polynomial& a,
	      gf_poly_multiplier& B, gf_poly_modulus & F)
{
	if (B.mod != &F || B.mod == 0)
		lidia_error_handler("gf_poly_multiplier", "multiply(...)::"
				    "wrong Fp_poly_modulus, or poly_multiplier not initialized");

	lidia_size_t n = F.modulus().degree();
	lidia_size_t da = a.degree();
	if (da >= n)
		lidia_error_handler("gf_poly_multiplier", "multiply(...)::"
				    "degree of Fp_polynomial must be < degree of Fp_poly_modulus");

	Fp_polynomial r1, r2, r3, r4, r5, r6;
	gf_polynomial P1;
	x.ffield = (a.ffield);
	P1.ffield = (a.ffield);
	gf_polynomial::build_frame(a.ffield);

	to_Kronecker(r1, a, 0, da);
	multiply(r2, r1, B.B_div_F);
	from_Kronecker(P1, r2, n-1, 2*n-3);

	multiply(r3, r1, B.B);
	to_Kronecker(r4, P1, 0, P1.degree());
	multiply(r5, r4, F.F);

	subtract(r6, r3, r5);
	from_Kronecker(x, r6, 0, n-1);
	P1.delete_frame();
	x.delete_frame();
}
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

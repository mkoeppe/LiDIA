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
//	Author	:Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/udigit.h"
#include	"LiDIA/base/interface_lib.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// define USE_SINGLE_PREC, if you want to use single prec. subroutines for
// single prec. moduli (see negate, add, subtract)
#define USE_SINGLE_PREC
//#undef USE_SINGLE_PREC



void Fp_polynomial::make_monic()
{
	debug_handler_l("Fp_polynomial", "make_monic (void)", 5);

	if (is_monic())
		return;

	lidia_size_t i, d = degree();
	const bigint &p = modulus();
	if (p.is_zero() || d < 0) {
		lidia_error_handler("Fp_polynomial", "make_monic (void) :: "
				    "modulus = 0  or  degree < 0");
		return;
	}

	bigint lc_inv;
	InvMod(lc_inv, coeff[d], p);
	for (i = d-1; i >= 0; i--)
		MulMod(coeff[i], coeff[i], lc_inv, p);
	coeff[d].assign_one();
}



//***************************************************************
//
//			operators
//
//***************************************************************

Fp_polynomial operator - (const Fp_polynomial &a)
{
	Fp_polynomial x;
	negate(x, a);
	return x;
}



Fp_polynomial operator + (const Fp_polynomial &a, const Fp_polynomial &b)
{
	Fp_polynomial x;
	add(x, b, a);
	return x;
}



Fp_polynomial operator + (const Fp_polynomial &a, const bigint &b)
{
	Fp_polynomial x;
	add(x, a, b);
	return x;
}



Fp_polynomial operator + (const bigint &a, const Fp_polynomial &b)
{
	Fp_polynomial x;
	add(x, a, b);
	return x;
}



Fp_polynomial operator - (const Fp_polynomial &a, const Fp_polynomial &b)
{
	Fp_polynomial x;
	subtract(x, a, b);
	return x;
}



Fp_polynomial operator - (const bigint &a, const Fp_polynomial &b)
{
	Fp_polynomial x;
	negate(x, b);
	add(x, x, a);
	return x;
}



Fp_polynomial operator - (const Fp_polynomial &a, const bigint &b)
{
	Fp_polynomial x;
	subtract(x, a, b);
	return x;
}



Fp_polynomial operator * (const Fp_polynomial &a, const Fp_polynomial &b)
{
	Fp_polynomial x;
	multiply(x, a, b);
	return x;
}



Fp_polynomial operator * (const bigint &a, const Fp_polynomial &b)
{
	Fp_polynomial x;
	multiply_by_scalar(x, b, a);
	return x;
}



Fp_polynomial operator * (const Fp_polynomial &a, const bigint &b)
{
	Fp_polynomial x;
	multiply_by_scalar(x, a, b);
	return x;
}



Fp_polynomial operator / (const Fp_polynomial &a, const Fp_polynomial &b)
{
	Fp_polynomial x;
	divide(x, a, b);
	return x;
}



Fp_polynomial operator / (const bigint &a, const Fp_polynomial &b)
{
	Fp_polynomial x;
	divide(x, a, b);
	return x;
}



Fp_polynomial operator / (const Fp_polynomial &a, const bigint &b)
{
	Fp_polynomial x;
	divide(x, a, b);
	return x;
}



Fp_polynomial operator % (const Fp_polynomial &a, const Fp_polynomial &b)
{
	Fp_polynomial x;
	remainder(x, a, b);
	return x;
}



Fp_polynomial & Fp_polynomial::operator += (const Fp_polynomial &a)
{
	add(*this, *this, a);
	return *this;
}



Fp_polynomial & Fp_polynomial::operator += (const bigint &a)
{
	add(*this, *this, a);
	return *this;
}



Fp_polynomial & Fp_polynomial::operator -= (const Fp_polynomial &a)
{
	subtract(*this, *this, a);
	return *this;
}



Fp_polynomial & Fp_polynomial::operator -= (const bigint &a)
{
	subtract(*this, *this, a);
	return *this;
}



Fp_polynomial & Fp_polynomial::operator *= (const Fp_polynomial &a)
{
	multiply(*this, *this, a);
	return *this;
}



Fp_polynomial & Fp_polynomial::operator *= (const bigint &a)
{
	multiply(*this, *this, a);
	return *this;
}



Fp_polynomial & Fp_polynomial::operator /= (const Fp_polynomial &a)
{
	divide(*this, *this, a);
	return *this;
}



Fp_polynomial & Fp_polynomial::operator /= (const bigint &a)
{
	divide(*this, *this, a);
	return *this;
}



Fp_polynomial & Fp_polynomial::operator %= (const Fp_polynomial &a)
{
	remainder(*this, *this, a);
	return *this;
}



bigint Fp_polynomial::operator() (const bigint & value) const
{
	debug_handler_l ("Fp_polynomial", "operator () (bigint)", 5);

	if (Mp == 0) {
		lidia_error_handler("Fp_polynomial", "operator (bigint &)::modulus = 0");
		return bigint();
	}

	bigint result = 0;
	lidia_size_t deg = degree();

	if (value.is_zero())
		result.assign(const_term());
	else
		if (deg >= 0) {
			const bigint &p = modulus();
			result.assign(coeff[deg]);
			lidia_size_t i;
			for (i = deg - 1; i >= 0; i--) {
				multiply(result, result, value); //result *= value;
				add(result, result, coeff[i]); //result += coeff[i];
				remainder(result, result, p);
			}
			if (result.is_negative())
				add(result, result, p);
		}
	return result;
}



//***************************************************************
//
//	      	addition, subtraction
//
//***************************************************************

void Fp_polynomial::negate()
{
	debug_handler_l ("Fp_polynomial", "negate (void)", 5);
	LiDIA::negate(*this, *this);
}



void negate(Fp_polynomial &a, const Fp_polynomial &b)
{
	debug_handler_l ("Fp_polynomial", "negate (Fp_polynomial&, Fp_polynomial&)", 5);

	if (b.Mp == 0) {
		lidia_error_handler("Fp_polynomial", "negate(Fp_polynomial&, Fp_polynomial&)::modulus = 0");
		return;
	}

	if (&a != &b) {
		a.set_modulus(b);
		a.set_degree(b.degree());
	}

	lidia_size_t i;
	const bigint &p = b.modulus();
#ifdef USE_SINGLE_PREC
	register udigit pp = p.least_significant_digit();

	if (p.compare(pp) == 0 && pp < max_udigit_modulus()) {
		for (i = 0; i < b.c_length; i++)
			a.coeff[i].assign(subtract_mod(static_cast<udigit>(0),
							      b.coeff[i].least_significant_digit(), pp));
	}
	else
#endif
	{
		for (i = 0; i < b.c_length; i++)
			NegateMod(a.coeff[i], b.coeff[i], p);
	}
}



void add(Fp_polynomial &c,
	 const Fp_polynomial &a, const Fp_polynomial &b)
{
	debug_handler_l ("Fp_polynomial",
			 "add (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 5);

	a.comp_modulus(b, "add");

	const bigint *ap, *bp;
	bigint *cp;
	lidia_size_t deg_a = a.degree(), deg_b = b.degree();
	lidia_size_t i, min_deg_ab, max_deg_ab;

	if (deg_a > deg_b) {
		max_deg_ab = deg_a;
		min_deg_ab = deg_b;
	}
	else {
		max_deg_ab = deg_b;
		min_deg_ab = deg_a;
	}

	if (&c != &a && &c != &b)
		c.set_modulus(a);
	c.set_degree(max_deg_ab);

	ap = a.coeff;
	bp = b.coeff;
	cp = c.coeff;
	const bigint & p = a.modulus();
#ifdef USE_SINGLE_PREC
	register udigit pp = p.least_significant_digit();

	if (p.compare(pp) == 0 && pp < max_udigit_modulus()) {
		for (i = min_deg_ab + 1; i; i--, ap++, bp++, cp++)
			cp->assign(add_mod(ap->least_significant_digit(),
						  bp->least_significant_digit(), pp));
	}
	else
#endif
	{
		for (i = min_deg_ab + 1; i; i--, ap++, bp++, cp++)
			AddMod((*cp), (*ap), (*bp), p);
	}

	if (deg_a > min_deg_ab && &c != &a)
		for (i = deg_a - min_deg_ab; i; i--, cp++, ap++)
			*cp = *ap;
	else if (deg_b > min_deg_ab && &c != &b)
		for (i = deg_b - min_deg_ab; i; i--, cp++, bp++)
			*cp = *bp;
	else
		c.remove_leading_zeros();
}



void add(Fp_polynomial &c, const Fp_polynomial &a, const bigint &b)
{
	debug_handler_l ("Fp_polynomial",
			 "add (Fp_polynomial&, Fp_polynomial&, bigint&)", 5);

	c.assign(a);

	bigint tmp;
	add(tmp, c[0], b);

	if (c.degree() < 0)
		c.set_degree(0);
	Remainder(c.coeff[0], tmp, c.modulus());

	if (c.degree() == 0)
		c.remove_leading_zeros();
}



void subtract(Fp_polynomial &c,
	      const Fp_polynomial &a, const Fp_polynomial &b)
{
	debug_handler_l ("Fp_polynomial",
			 "subtract (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 5);

	a.comp_modulus(b, "subtract");

	const bigint *ap, *bp;
	bigint *cp;
	lidia_size_t deg_a = a.degree(), deg_b = b.degree();
	lidia_size_t i, min_deg_ab, max_deg_ab;

	if (deg_a > deg_b) {
		max_deg_ab = deg_a;
		min_deg_ab = deg_b;
	}
	else {
		max_deg_ab = deg_b;
		min_deg_ab = deg_a;
	}

	if (&c != &a && &c != &b)
		c.set_modulus(a);
	c.set_degree(max_deg_ab);

	ap = a.coeff;
	bp = b.coeff;
	cp = c.coeff;
	const bigint & p = a.modulus();
#ifdef USE_SINGLE_PREC
	register udigit pp = p.least_significant_digit();
	bool sgl_prec = (p.compare(pp) == 0 && pp < max_udigit_modulus());

	if (sgl_prec) {
		for (i = min_deg_ab + 1; i; i--, ap++, bp++, cp++)
			cp->assign(subtract_mod(ap->least_significant_digit(),
						       bp->least_significant_digit(), pp));
	}
	else
#endif
	{
		for (i = min_deg_ab + 1; i; i--, ap++, bp++, cp++)
			SubMod((*cp), (*ap), (*bp), p);
	}

	if (deg_a > min_deg_ab && &c != &a)
		for (i = deg_a - min_deg_ab; i; i--, cp++, ap++)
			*cp = *ap;
	else if (deg_b > min_deg_ab) {
#ifdef USE_SINGLE_PREC
		if (sgl_prec) {
			for (i = deg_b - min_deg_ab; i; i--, cp++, bp++)
				cp->assign(subtract_mod(static_cast<udigit>(0),
							       bp->least_significant_digit(), pp));
		}
		else
#endif
		{
			for (i = deg_b - min_deg_ab; i; i--, cp++, bp++)
				NegateMod((*cp), (*bp), p);
		}
	}
	else
		c.remove_leading_zeros();
}



void subtract(Fp_polynomial & x, const Fp_polynomial& a, const bigint& b)
{
	debug_handler_l ("Fp_polynomial", "subtract(Fp_polynomial &, "
			 "const Fp_polynomial&, const bigint&)", 5);
	add(x, a, -b);
}



//***************************************************************
//
//			multiplication routines
//
//***************************************************************

void multiply_by_scalar(Fp_polynomial &c,
			const Fp_polynomial &a, const bigint &b)
	//b may even be any of a's or c's coefficients
{
	debug_handler_l ("Fp_polynomial",
			 "multiply_by_scalar (Fp_polynomial&, Fp_polynomial&, bigint&)", 5);

	if (a.Mp == 0) {
		lidia_error_handler("Fp_polynomial", "multiply_by_scalar (Fp_polynomial&, Fp_polynomial&, bigint&)::modulus = 0");
		return;
	}

	const bigint& p = a.modulus();
	bigint tmp_b;
	Remainder(tmp_b, b, p);

	lidia_size_t i, deg_a = a.degree();

	c.set_modulus(a); //assigns zero

	if (!tmp_b.is_zero() && deg_a >= 0) {
		c.set_degree(deg_a);

		const bigint *ap;
		bigint *cp;
		for (i = 0, ap = a.coeff, cp = c.coeff; i <= deg_a; i++, ap++, cp++)
			MulMod((*cp), (*ap), tmp_b, p);
	}
}



void multiply(Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler_l("Fp_polynomial", "multiply (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 5);

	a.comp_modulus(b, "multiply");

	lidia_size_t crov =
		Fp_polynomial::crossovers.fftmul_crossover(a.modulus());
	if (a.degree() > crov && b.degree() > crov)
		fft_mul(x, a, b);
	else
		plain_mul(x, a, b);
}



void square(Fp_polynomial& x, const Fp_polynomial& a)
{
	debug_handler_l("Fp_polynomial", "square (Fp_polynomial&, Fp_polynomial&)", 5);

	lidia_size_t crov =
		Fp_polynomial::crossovers.fftmul_crossover(a.modulus());
	if (a.degree() > crov)  fft_sqr(x, a);
	else                    plain_sqr(x, a);
}



//***************************************************************
//
//			division routines
//
//***************************************************************

void div_rem(Fp_polynomial& q, Fp_polynomial& r, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler_l("Fp_polynomial", "div_rem(Fp_polynomial&, "
			"Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 5);

	a.comp_modulus(b, "div_rem");

	lidia_size_t crov =
		Fp_polynomial::crossovers.fftdiv_crossover(a.modulus());
	if (b.degree() > crov && a.degree() - b.degree() > crov)
		fft_div_rem(q, r, a, b);
	else
		plain_div_rem(q, r, a, b);
}



void divide(Fp_polynomial& q, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler_l("Fp_polynomial", "divide (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 5);

	a.comp_modulus(b, "divide");

	lidia_size_t crov =
		Fp_polynomial::crossovers.fftdiv_crossover(a.modulus());
	if (b.degree() > crov && a.degree() - b.degree() > crov)
		fft_div(q, a, b);
	else
		plain_div(q, a, b);
}



void divide(Fp_polynomial& q, const Fp_polynomial& a, const bigint& b)
{
	debug_handler_l("Fp_polynomial", "divide (Fp_polynomial&, Fp_polynomial&, bigint&)", 5);
	bigint tmp;
	InvMod(tmp, b, a.modulus()); //raises error if tmp == 0
	multiply_by_scalar(q, a, tmp);
}



void divide(Fp_polynomial& q, const bigint& a, const Fp_polynomial& b)
{
	Fp_polynomial tmp;
	tmp.set_modulus(b);
	tmp = a;
	divide(q, tmp, b);
}



void remainder(Fp_polynomial& r, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler_l("Fp_polynomial", "remainder (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 5);

	a.comp_modulus(b, "remainder");

	lidia_size_t crov =
		Fp_polynomial::crossovers.fftrem_crossover(a.modulus());
	if (b.degree() > crov && a.degree() - b.degree() > crov)
		fft_rem(r, a, b);
	else
		plain_rem(r, a, b);
}



void invert(Fp_polynomial& x, const Fp_polynomial& a, lidia_size_t m)
{
	debug_handler_l("Fp_polynomial", "inv(Fp_polynomial&, Fp_polynomial&, lidia_size_t)", 5);

	lidia_size_t crov = comparator< lidia_size_t >::max(
		Fp_polynomial::crossovers.inv_crossover(a.modulus()),
		Fp_polynomial::crossovers.log2_newton_crossover(a.modulus()));
	if (&x == &a) {
		Fp_polynomial la(a);
		if (m > crov)  newton_inv(x, la, m);
		else           plain_inv(x, la, m);
	}
	else {
		if (m > crov)  newton_inv(x, a, m);
		else           plain_inv(x, a, m);
	}
}



//***************************************************************
//
//		Modular Arithmetic without pre-conditioning
//
//***************************************************************

// arithmetic mod f.
// all inputs and outputs are polynomials of degree less than deg(f).
// ASSUMPTION: f is assumed monic, and deg(f) > 0.


void multiply_mod(Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& b, const Fp_polynomial& f)
{
	debug_handler_l("Fp_polynomial", "multiply_mod(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 6);

	a.comp_modulus(b, "multiply_mod");
	a.comp_modulus(f, "multiply_mod");

	Fp_polynomial t;
	multiply(t, a, b);
	remainder(x, t, f);
}



void square_mod(Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& f)
{
	debug_handler_l("Fp_polynomial", "square_mod(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 6);

	a.comp_modulus(f, "square_mod");

	Fp_polynomial t;
	square(t, a);
	remainder(x, t, f);
}



void multiply_by_x_mod(Fp_polynomial& h, const Fp_polynomial& a, const Fp_polynomial& f)
{
	debug_handler_l("Fp_polynomial", "multiply_by_x_mod(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 6);

	a.comp_modulus(f, "multiply_by_x_mod");

	lidia_size_t i, n, m;

	n = f.degree();
	m = a.degree();

	if (n > 0 && m == n-1 && f.is_monic() && &f != &h) {
		h.set_modulus(f);
		h.set_degree(n-1);
		const bigint& p = f.modulus();
		bigint z, t, *hh;
		const bigint *aa, *ff;
		hh = h.coeff;
		aa = a.coeff;
		ff = f.coeff;
		negate(z, aa[n-1]);
		for (i = n-1; i >= 1; i--) {
			MulMod(t, z, ff[i], p);
			AddMod(hh[i], aa[i-1], t, p);
		}
		MulMod(hh[0], z, ff[0], p);
		h.remove_leading_zeros();
	}
	else {
		if (m < n-1)
			shift_left(h, a, 1);
		else {
			Fp_polynomial tmp;
			shift_left(tmp, a, 1);
			remainder(h, tmp, f);
		}
	}
}



void invert_mod(Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& f)
{
	debug_handler_l("Fp_polynomial", "invert_mod(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 6);

	a.comp_modulus(f, "invert_mod");

	Fp_polynomial d, t;
	xgcd(d, x, t, a, f);
	if (!d.is_one())
		lidia_error_handler("Fp_polynomial", "invert_mod(Fp_polynomial&, "
				    "Fp_polynomial&, Fp_polynomial&)::can't compute multiplicative "
				    "inverse");
}



bool invert_mod_status(Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& f)
{
	debug_handler_l("Fp_polynomial", "invert_mod_status(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)", 6);

	a.comp_modulus(f, "invert_mod_status");

	Fp_polynomial d, t;
	xgcd(d, x, t, a, f);
	if (!d.is_one())
		x.assign(d);
	return (d.is_one());
}



// WARNING: obsolete.  Use power with Fp_poly_modulus.
void power_mod(Fp_polynomial& h, const Fp_polynomial& g, const bigint& e, const Fp_polynomial& f)
{
	debug_handler_l("Fp_polynomial", "power_mod(Fp_polynomial&, Fp_polynomial&, bigint&, Fp_polynomial&)", 6);

	g.comp_modulus(f, "power_mod");

	Fp_polynomial lg(g);
	lidia_size_t i, n = e.bit_length();

	Fp_polynomial tmp;
	tmp.set_modulus(f);
	tmp.set_max_degree(f.degree());
	tmp.assign_one();

	for (i = n - 1; i >= 0; i--) {
		square_mod(tmp, tmp, f);
		if (e.bit(i))
			multiply_mod(tmp, tmp, lg, f);
	}
	h.assign(tmp);
}



// WARNING: obsolete.  Use power_x with Fp_poly_modulus.
void power_x_mod(Fp_polynomial& h, const bigint& e, const Fp_polynomial& f)
{
	debug_handler_l("Fp_polynomial", "power_x_mod(Fp_polynomial&, bigint&, Fp_polynomial&)", 6);

	lidia_size_t i, n = e.bit_length();

	Fp_polynomial tmp;
	tmp.set_modulus(f);
	tmp.set_max_degree(f.degree());
	tmp.assign_one();

	for (i = n - 1; i >= 0; i--) {
		square_mod(tmp, tmp, f);
		if (e.bit(i))
			multiply_by_x_mod(tmp, tmp, f);
	}
	h.assign(tmp);
}



// WARNING: obsolete.  Use power_x_plus_a with Fp_poly_modulus.
void power_x_plus_a_mod(Fp_polynomial& h, const bigint& a, const bigint& e, const Fp_polynomial& f)
{
	debug_handler_l("Fp_polynomial", "power_x_plus_a_mod(Fp_polynomial&, bigint&, bigint&, Fp_polynomial&)", 6);

	lidia_size_t j = f.degree()-1;
	if (j == -1) {
		h.set_modulus(f); //assigns zero
		return;
	}
	Fp_polynomial t1, t2;
	t1.set_max_degree(j);
	t2.set_max_degree(j);

	Fp_polynomial tmp;
	tmp.set_modulus(f);
	tmp.set_max_degree(j);
	tmp.assign_one();

	lidia_size_t i, n = e.bit_length();

	for (i = n - 1; i >= 0; i--) {
		square_mod(tmp, tmp, f);
		if (e.bit(i)) {
			multiply_by_x_mod(t1, tmp, f);
			multiply_by_scalar(t2, tmp, a);
			add(tmp, t1, t2);
		}
	}
	h.assign(tmp);
}



//***************************************************************
//
//				   gcd's
//
//***************************************************************
//see  "gcd.cc"




//***************************************************************
//
//				other functions
//
//***************************************************************


// x = a % X^m
void trunc(Fp_polynomial &x, const Fp_polynomial &a, lidia_size_t m)
{
	debug_handler_l ("Fp_polynomial",
			 "trunc (Fp_polynomial&, Fp_polynomial&, lidia_size_t)", 5);
	if (&x == &a) {
		if (x.c_length > m)
			x.c_length = m;
	}
	else {
		x.set_modulus(a);
		lidia_size_t i, n = comparator< lidia_size_t >::min(a.c_length, m);
		x.set_degree(n-1);

		bigint* xp = x.coeff;
		const bigint* ap = a.coeff;

		for (i = 0; i < n; i++)
			xp[i].assign(ap[i]);
	}
	x.remove_leading_zeros();
}



// x = a/X^n
void shift_right(Fp_polynomial &x, const Fp_polynomial &a, lidia_size_t n)
{
	debug_handler_l ("Fp_polynomial",
			 "shift_right (Fp_polynomial&, Fp_polynomial&, lidia_size_t)", 5);

	if (&x != &a)
		x.set_modulus(a);
	lidia_size_t i, da = a.degree();

	if (da < n) {
		x.assign_zero();
		return;
	}

	if (&x != &a)
		x.set_max_degree(da-n);

	for (i = 0; i <= da-n; i++)
		x.coeff[i].assign(a.coeff[i+n]);

	x.c_length = da - n + 1;
}



// x = a*X^n
void shift_left(Fp_polynomial &x, const Fp_polynomial &a, lidia_size_t n)
{
	debug_handler_l ("Fp_polynomial",
			 "shift_left (Fp_polynomial&, Fp_polynomial&, lidia_size_t)", 5);

	if (&x != &a)
		x.set_modulus(a);

	if (a.is_zero()) {
		x.assign_zero();
		return;
	}

	lidia_size_t i, m = a.c_length;
	x.set_degree(m+n-1);

	for (i = m-1; i >= 0; i--)
		x.coeff[i+n].assign(a.coeff[i]);

	for (i = 0; i < n; i++)
		(x.coeff[i]).assign_zero();
}



// x = derivative of a
void derivative(Fp_polynomial &x, const Fp_polynomial &a)
{
	debug_handler_l ("Fp_polynomial",
			 "derivative (Fp_polynomial&, Fp_polynomial&)", 5);

	if (a.Mp == 0) {
		lidia_error_handler("Fp_polynomial", "derivative (Fp_polynomial&, Fp_polynomial&)::modulus = 0");
		return;
	}

	if (&x != &a)
		x.set_modulus(a);

	lidia_size_t d = a.degree();
	if (d <= 0) {
		x.assign_zero();
		return;
	}

	const bigint& p = a.modulus();
	x.set_degree(d - 1);
	lidia_size_t i;
	for (i = 0; i < d; i++) {
		multiply(x.coeff[i], a.coeff[i+1], i+1);
		Remainder(x.coeff[i], x.coeff[i], p);
	}
	x.remove_leading_zeros();
}



// computes x = a mod X^m-1
void cyclic_reduce(Fp_polynomial& x, const Fp_polynomial& a, lidia_size_t m)
{
	debug_handler_l("Fp_polynomial", "cyclic_reduce(Fp_polynomial&, Fp_polynomial&, lidia_size_t)", 6);
	lidia_size_t n = a.degree();
	lidia_size_t i, j;

	if (a.Mp == 0) {
		lidia_error_handler("Fp_polynomial", "cyclic_reduce(Fp_polynomial&, Fp_polynomal&, lidia_size_t)::modulus = 0");
		return;
	}

	if (n < m) {
		x.assign(a);
		return;
	}

	if (&x != &a) {
		x.set_modulus(a);
		x.set_degree(m-1);
	}

	const bigint &p = a.modulus();
	bigint accum;
	for (i = 0; i < m; i++) {
		accum.assign(a.coeff[i]);
		for (j = i + m; j <= n; j += m)
			add(accum, accum, a.coeff[j]);
		Remainder(x.coeff[i], accum, p);
	}

	if (&x == &a)
		x.set_degree(m-1);

	x.remove_leading_zeros();
}



void power(Fp_polynomial &x, const Fp_polynomial &a, lidia_size_t e)
{
	if (e < 0) {
		lidia_error_handler("Fp_polynomial", "power(Fp_polynomial&, "
				    "Fp_polynomial&, lidia_size_t)::exponent must not be negative");
		return;
	}

	if (e == 1) {
		x.assign(a);
		return;
	}

	Fp_polynomial m(a);
	int i, n = integer_log(e);
	x.set_modulus(a);
	x.assign_one();
	x.set_max_degree(a.degree() * e);

	for (i = n - 1; i >= 0; i--) {
		square(x, x);
		if (e & (1 << i))
			multiply(x, x, m);
	}
}



void
add_multiple(Fp_polynomial &f, const Fp_polynomial &g, const bigint &s,
	     lidia_size_t n, const Fp_polynomial &h)
	//f = g + s*x^n*h
{
	g.comp_modulus(h, "add_multiple");
	if (n < 0) {
		lidia_error_handler("Fp_polynomial",
				    "add_multiple(...)::bad argument");
		return;
	}

	if (h.is_zero() || s.is_zero()) {
		f.assign(g);
		return;
	}
	if (g.is_zero()) {
		shift_left(f, h, n);
		return;
	}
//    if (s.is_one() && n == 0)
//    {
//	add(f, g, h);
//	return;
//    }


	const bigint *gp, *hp;
	bigint *fp;

	lidia_size_t deg_g = g.degree(), deg_h = h.degree();
	lidia_size_t i, max_deg, min_deg;

	if (deg_g < deg_h+n) {
		min_deg = deg_g;
		max_deg = deg_h + n;
	}
	else {
		min_deg = deg_h + n;
		max_deg = deg_g;
	}

	Fp_polynomial hh;
	if (&f == &h) {
		hh.assign(h);
		hp = hh.coeff;
	}
	else
		hp = h.coeff;

	if (&f != &g && &f != &h)
		f.set_modulus(g);
	f.set_degree(max_deg);
	fp = f.coeff;
	gp = g.coeff;

	bigint tmp;

	if (&f != &g) {
		fp = f.coeff;
		gp = g.coeff;
		for (i = n; i; i--)
			*(fp++) = *(gp++);
	}
	else
		gp = fp = f.coeff + n;

	const bigint & p = g.modulus();
	for (i = min_deg-n+1; i; i--) {
		multiply(tmp, s, *(hp++));
		add(tmp, tmp, *(gp++));
		Remainder(*(fp++), tmp, p);
	}

	if (deg_h + n > min_deg)
		for (i = max_deg - min_deg; i; i--)
			MulMod(*(fp++), *(hp++), s, p);
	else if (deg_g > min_deg && &f != &g)
		for (i = max_deg - min_deg; i; i--)
			*(fp++) = *(gp++);
	else
		f.remove_leading_zeros();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

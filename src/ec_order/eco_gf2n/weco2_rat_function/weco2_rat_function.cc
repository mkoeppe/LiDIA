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
//	Author	: Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/weco2_rat_function.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// ---------- STATIC MEMBERS ----------

weco2_rat_function::ff_element *(weco2_rat_function::a6) = NULL;
weco2_rat_function::ff_pol_m weco2_rat_function::modulus;
weco2_rat_function::ff_pol weco2_rat_function::curve;



// ---------- CONSTRUCTORS / DESTRUCTORS ----------

weco2_rat_function::weco2_rat_function()
{
	zero = true;
}



weco2_rat_function::
weco2_rat_function(const ff_rat & xx, const ff_rat & yy)
{
	zero = false;
	x.assign(xx); y0.assign(yy); y1.assign_zero();
#ifdef DEBUG
	assert(on_curve());
#endif
}



weco2_rat_function::
weco2_rat_function(const ff_rat & xx, const ff_rat & yy0, const ff_rat & yy1)
{
	zero = false;
	x.assign(xx); y0.assign(yy0); y1.assign(yy1);
#ifdef DEBUG
	assert(on_curve());
#endif
}



weco2_rat_function::
weco2_rat_function(const weco2_rat_function & P)
{
	zero = P.zero;
	x.assign(P.x); y0.assign(P.y0); y1.assign(P.y1);
#ifdef DEBUG
	assert(on_curve());
#endif
}



weco2_rat_function::~weco2_rat_function()
{
}



// ---------- COMPARISONS ----------

bool operator == (const weco2_rat_function & P, const weco2_rat_function & Q)
{
	if (P.zero != Q.zero)
		return false;
	if (P.zero == true)
		return true;
	else
		return (equal(P.x, Q.x, weco2_rat_function::modulus) &&
			equal(P.y0, Q.y0, weco2_rat_function::modulus) &&
			equal(P.y1, Q.y1, weco2_rat_function::modulus));
}



bool weco2_rat_function::on_curve() const
{
	ff_rat h0, h1;

	if (is_zero())
		return true;

	if (!y1.is_zero ()) {              // test y1 * X + x(P) == 0
		shift(h1, y1, 1);
		h1.reduce(weco2_rat_function::modulus);
		add(h1, h1, x, weco2_rat_function::modulus);

		if (!h1.is_zero())
			return false;
	}

	add(h0, y0, x, weco2_rat_function::modulus);
	multiply(h0, h0, y0, weco2_rat_function::modulus);

	square(h1, y1, weco2_rat_function::modulus);
	multiply(h1, h1, weco2_rat_function::curve, weco2_rat_function::modulus);
	add(h0, h1, h0, weco2_rat_function::modulus);

	square(h1, x, weco2_rat_function::modulus);
	multiply(h1, h1, x, weco2_rat_function::modulus);
	add(h0, h0, h1, weco2_rat_function::modulus);

	add(h0, h0, gf2n_polynomial(*a6), weco2_rat_function::modulus);

	if (!h0.is_zero())
		return false;
	else
		return true;
}



// ---------- ARITHMETICAL OPERATIONS ----------

void negate (weco2_rat_function & R, const weco2_rat_function & P)
{
	if (P.is_zero())
		R.assign_zero();
	else {
		weco2_rat_function::ff_rat h1;
		add(h1, P.y0, P.x, weco2_rat_function::modulus);
		R.assign(P.x, h1, P.y1);
#ifdef DEBUG
		assert(R.on_curve());
#endif
	}
}



//***************************************************************
//  REMARK:  the following addition functions can probably improved
//           if the rational function operations are substituted
//           by formulas using Fp_polynomials (maybe one or two
//           multiplications are done twice at the moment)
//***************************************************************

// the non-trivial case of addition

void add_PQ (weco2_rat_function & R, const weco2_rat_function & P,
   	     const weco2_rat_function & Q)
{

#ifdef DEBUG
	assert(P.on_curve());
	assert(Q.on_curve());
#endif

	weco2_rat_function::ff_rat h, h1, h2, h3, lambda0;

	add(h1, P.x, Q.x, weco2_rat_function::modulus);
	add(lambda0, P.y0, Q.y0, weco2_rat_function::modulus);
	divide(lambda0, lambda0, h1, weco2_rat_function::modulus);

	// we know that lambda = lambda0 + Y/X

	add(h, lambda0, gf2n(1));
	multiply(h, lambda0, h, weco2_rat_function::modulus);
	shift(h, h, 2);
	add(h, h, weco2_rat_function::curve);
	shift(h3, h, -2);
	h3.reduce(weco2_rat_function::modulus);
	add(h3, h3, Q.x, weco2_rat_function::modulus); // h3 = x3 + x(P)
	add(h, h3, P.x, weco2_rat_function::modulus); // x-coordinate finished

	multiply(h2, h3, lambda0, weco2_rat_function::modulus);
	add(h2, h2, h, weco2_rat_function::modulus);
	add(h2, h2, P.y0, weco2_rat_function::modulus); // y0-part finished

	shift(lambda0, h3, -1);
	lambda0.reduce(weco2_rat_function::modulus);
	add(lambda0, lambda0, P.y1, weco2_rat_function::modulus); // y1-part finished

	R.assign(h, h2, lambda0);

#ifdef DEBUG
	assert(R.on_curve());
#endif
}



//****************************************************************

void add (weco2_rat_function & R, const weco2_rat_function & Q,
          const weco2_rat_function & P)
{
	if (P.is_zero()) {
		if (Q.is_zero())
			R.assign_zero();
		else
			R.assign(Q);
	}
	else {
		if (Q.is_zero())
			R.assign(P);
		else
			if (Q == P)
				multiply_by_2 (R, Q);
			else
				if (equal(P.x, Q.x, weco2_rat_function::modulus))
					R.assign_zero();
				else
					add_PQ(R, P, Q);
	}
}



void subtract (weco2_rat_function & R, const weco2_rat_function & Q,
               const weco2_rat_function & P)
{
	if (P.is_zero()) {
		if (Q.is_zero())
			R.assign_zero();
		else
			R.assign(Q);
	}
	else {
		if (Q.is_zero())
			negate(R, P);
		else if (Q == P)
			R.assign_zero();
		else
			if (equal(P.x, Q.x, weco2_rat_function::modulus))
				multiply_by_2(R, Q);
			else {
				weco2_rat_function tmp;
				negate(tmp, P);
				add_PQ(R, tmp, Q);
			}
	}
}



void multiply_by_2 (weco2_rat_function & Q, const weco2_rat_function & P)
{
	if (P.is_zero() || P.x.is_zero()) {
		Q.assign_zero();
		return;
	}

	weco2_rat_function::ff_rat h, h1, h2, h3, h4;

	square(h, P.x, weco2_rat_function::modulus); // h = x^2
	square(h1, h, weco2_rat_function::modulus);
	add(h1, h1, *(weco2_rat_function::a6)); // h1 = (x^4 + a6)
	divide(h2, h1, h, weco2_rat_function::modulus); // x-coordinate

	divide(h3, P.y1, P.x, weco2_rat_function::modulus);
	multiply(h3, h3, h2, weco2_rat_function::modulus); // y1-part

	divide(h4, P.y0, P.x, weco2_rat_function::modulus);
	add(h4, h4, P.x, weco2_rat_function::modulus);
	multiply(h4, h4, h1, weco2_rat_function::modulus);
	add(h4, h4, *(weco2_rat_function::a6));
	divide(h4, h4, h, weco2_rat_function::modulus); // y0-part

	Q.assign(h2, h4, h3);

#ifdef DEBUG
	assert(Q.on_curve());
#endif
}



//---------------------------------------------------------
//  The following function implements a left-to-right fast
//  exponentiation which bits 0, +1, -1. See the paper "Fast
//  multiplication on elliptic curves" by Volker Mueller


void multiply2(weco2_rat_function & Q, lidia_size_t n,
	       const weco2_rat_function & PP)
{
	weco2_rat_function P(PP);
	lidia_size_t mask;
	unsigned char state = 0;

	if (n < 0) {
		negate(P, P);
		n = -n;
	}

	Q.assign_zero();

	if (n == 0)
		return;

	if (n == 1) {
		Q.assign(P);
		return;
	}

	if (n < 65537)
		mask = 1 << 16;
	else
		mask = 1 << 31;

	while (!(n & mask))
		mask >>= 1;

	while (mask > 0) {
		switch (state) {
		case 0:
			if (mask & n)
				state = 1;
			else
				multiply_by_2(Q, Q);
			break;

		case 1:
			if (mask & n) {
				add(Q, Q, P);
				multiply_by_2(Q, Q);
				multiply_by_2(Q, Q);
				state = 2;
			}
			else {
				multiply_by_2(Q, Q);
				add(Q, Q, P);
				multiply_by_2(Q, Q);
				state = 0;
			}
			break;

		case 2:
			if (mask & n)
				multiply_by_2(Q, Q);
			else
				state = 3;
			break;

		case 3:
			if (mask & n) {
				multiply_by_2(Q, Q);
				subtract(Q, Q, P);
				multiply_by_2(Q, Q);
				state = 2;
			}
			else {
				subtract(Q, Q, P);
				multiply_by_2(Q, Q);
				multiply_by_2(Q, Q);
				state = 0;
			}
		}
		mask >>= 1;
	}

	switch (state) {
	case 0:
		break;
	case 1:
		multiply_by_2(Q, Q);
		add(Q, Q, P);
		break;
	case 2:
		subtract(Q, Q, P);
		break;
	case 3:
		subtract(Q, Q, P);
		multiply_by_2(Q, Q);
		break;
	}
}



void multiply (weco2_rat_function & R1, lidia_size_t n1,
               const weco2_rat_function & PP)
{
	weco2_rat_function P (PP);

	if (n1 < 0) {
		n1 = -n1; negate(P, PP);
	}

	R1.assign_zero();

	while (n1 > 0) {
		if (n1 % 4 == 3) {
			subtract(R1, R1, P);
			n1 ++;
		}
		else
			if (n1 & 1)
				add(R1, R1, P);

		n1 >>= 1;

		if (n1 > 0)
			multiply_by_2(P, P);
	}
}



// R1 = n1 * P, R2 = n2 * P, n1, n2 >= 0

void multiply (weco2_rat_function & R1, weco2_rat_function & R2,
               lidia_size_t n1, lidia_size_t n2,
               const weco2_rat_function & PP)
{

	weco2_rat_function P (PP);
	R1.assign_zero(); R2.assign_zero();

	if (n1 < 0 || n2 < 0)
		lidia_error_handler("weco2_rat_function", "multiply::either n1 or n2 is "
				    "negative");

	while (n1 > 0 || n2 > 0) {
		if (n1 % 4 == 3) {
			subtract(R1, R1, P);
			n1 ++;
		}
		else
			if (n1 & 1)
				add(R1, R1, P);
		n1 >>= 1;

		if (n2 % 4 == 3) {
			subtract(R2, R2, P);
			n2 ++;
		}
		else
			if (n2 & 1)
				add(R2, R2, P);
		n2 >>= 1;

		if (n1 > 0 || n2 > 0)
			multiply_by_2(P, P);
	}
}



// ---------- ASSIGNMENTS ----------

weco2_rat_function &
weco2_rat_function::operator = (const weco2_rat_function & wep)
{
	if (this != & wep) {
		if (wep.zero)
			this->zero = true;
		else {
			this->assign (wep.x, wep.y0, wep.y1);
			this->zero=false;
		}
	}
	return *this;
}



void weco2_rat_function::assign_xy ()
{
	zero = false;
	y1.assign_one();
	y0.assign_zero();
	x.assign_x();
}



void weco2_rat_function::
assign (const ff_rat & xx, const ff_rat & yy0)
{
	x.assign(xx);
	y0.assign(yy0);
	y1.assign_zero();
	zero = false;
}



void weco2_rat_function::
assign (const ff_rat & xx, const ff_rat & yy0, const ff_rat & yy1)
{
	x.assign(xx);
	y0.assign(yy0);
	y1.assign(yy1);
	zero = false;
}



void weco2_rat_function::
assign (const weco2_rat_function & w)
{
	x.assign(w.x);
	y0.assign(w.y0);
	y1.assign(w.y1);
	zero = w.zero;
}



void weco2_rat_function::
initialize(const weco2_rat_function::ff_element & aa,
           const weco2_rat_function::ff_pol_m & pm)
{
	if (aa.is_zero())
		lidia_error_handler("weco2_rat_function", "initialize::curve has "
				    "j-invariant zero");
	if (a6 != NULL)
		delete a6;
	weco2_rat_function::a6 = new weco2_rat_function::ff_element(aa);

	weco2_rat_function::modulus = pm;
	weco2_rat_function::curve.assign(ff_pol(3));
	weco2_rat_function::curve.set_coefficient(aa, 0);
	remainder(weco2_rat_function::curve, weco2_rat_function::curve,
		  pm.modulus());
}



std::ostream & operator << (std::ostream & o, const weco2_rat_function & P)
{
	if (P.is_zero())
		o << "O";
	else
		o << "(" << P.x << ", " << P.y0 << " + y * " << P.y1 << ")";
	return o;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

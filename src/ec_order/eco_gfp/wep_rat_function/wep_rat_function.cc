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
//	Author	:Markus Maurer (MM), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/wep_rat_function.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// ---------- STATIC MEMBERS ----------

wep_rat_function::ff_element wep_rat_function::a;
wep_rat_function::ff_element wep_rat_function::b;
wep_rat_function::ff_pol_m wep_rat_function::modulus;


// ---------- CONSTRUCTORS / DESTRUCTORS ----------

wep_rat_function::wep_rat_function()
{
	x.set_modulus(wep_rat_function::a.modulus());
	x.assign_one();
	y.set_modulus(wep_rat_function::a.modulus());
	y.assign_one();
	z.set_modulus(wep_rat_function::a.modulus());
	z.assign_zero();
}



wep_rat_function::
wep_rat_function(const ff_pol & xx, const ff_pol & yy)
{
	x.assign(xx);
	y.assign(yy);
	z.assign_one();
	if (!on_curve())
		lidia_error_handler("wep_rat_function", "ct::input not on curve");
}



wep_rat_function::~wep_rat_function()
{
}




// ---------- COMPARISONS ----------

bool operator == (const wep_rat_function & P, const wep_rat_function & Q)
{
	wep_rat_function::ff_pol h, h2, h1, h3;

	h.set_modulus(wep_rat_function::a.modulus());
	h1.set_modulus(wep_rat_function::a.modulus());
	h3.set_modulus(wep_rat_function::a.modulus());
	h2.set_modulus(wep_rat_function::a.modulus());

	square(h1, Q.z, wep_rat_function::modulus);
	multiply(h, P.x, h1, wep_rat_function::modulus);
	square(h2, P.z, wep_rat_function::modulus);
	multiply(h3, h2, Q.x, wep_rat_function::modulus);
	subtract(h, h, h3);

	if (! h.is_zero())
		return false;

	multiply(h1, h1, Q.z, wep_rat_function::modulus);
	multiply(h2, P.z, h2, wep_rat_function::modulus);
	multiply(h, P.y, h1, wep_rat_function::modulus);
	multiply(h3, h2, Q.y, wep_rat_function::modulus);
	subtract(h, h, h3);

	if (h.is_zero())
		return true;
	else
		return false;
}



bool wep_rat_function::on_curve() const
{
	ff_pol h, h2, h1, h3;

	if (z.is_zero())
		return true;

	h.set_modulus(a.modulus());
	h1.set_modulus(a.modulus());
	h2.set_modulus(a.modulus());
	h3.set_modulus(a.modulus());

	square(h, y, wep_rat_function::modulus);
	mult_by_curve(h, wep_rat_function::modulus.modulus());

	square(h1, x, wep_rat_function::modulus);
	multiply(h1, x, h1, wep_rat_function::modulus);
	square(h2, z, wep_rat_function::modulus);
	square(h3, h2, wep_rat_function::modulus);
	add(h2, wep_rat_function::a.mantissa() * x,
	    wep_rat_function::b.mantissa() * h2);
	multiply(h2, h2, h3, wep_rat_function::modulus);
	add(h1, h1, h2);
	subtract(h, h, h1);

	return (h.is_zero());
}



// ---------- ARITHMETICAL OPERATIONS ----------
// internal function to multiply with x^3 + a *x + b, but with
// shifts (should be faster than a multiply)

void negate (wep_rat_function & R, const wep_rat_function & P)
{
	if (P.is_zero())
		R.assign_zero();
	else
		R.assign(P.x, -P.y, P.z);
}


//**************************************************************
// Formulae as in IEEE P1363, Annex A
// ==> 12 multiplications, 4 squarings
//**************************************************************

void add_PQ (wep_rat_function & R, const wep_rat_function & P,
   	     const wep_rat_function & Q)
{
	wep_rat_function::ff_pol t1, t2, t3, t4, t5;
	wep_rat_function::ff_pol x3, y3, z3;

	t1.set_modulus(wep_rat_function::a.modulus());
	t2.set_modulus(wep_rat_function::a.modulus());
	t3.set_modulus(wep_rat_function::a.modulus());
	t4.set_modulus(wep_rat_function::a.modulus());
	t5.set_modulus(wep_rat_function::a.modulus());
	x3.set_modulus(wep_rat_function::a.modulus());
	y3.set_modulus(wep_rat_function::a.modulus());
	z3.set_modulus(wep_rat_function::a.modulus());

	square(t2, Q.z, wep_rat_function::modulus);
	multiply(t1, P.x, t2, wep_rat_function::modulus);
	multiply(t2, t2, Q.z, wep_rat_function::modulus);
	multiply(t2, P.y, t2, wep_rat_function::modulus);

	square(t4, P.z, wep_rat_function::modulus);
	multiply(t3, t4, Q.x, wep_rat_function::modulus);
	multiply(t4, t4, P.z, wep_rat_function::modulus);
	multiply(t4, Q.y, t4, wep_rat_function::modulus);

	subtract(t5, t1, t3);
	add(t1, t1, t3);
	subtract(t3, t2, t4);
	add(t2, t2, t4);

	if (t5.is_zero()) {
		if (t3.is_zero())
			multiply_by_2(R, P);
		else
			R.assign_zero();
		return;
	}

	multiply(z3, t5, Q.z, wep_rat_function::modulus);
	multiply(z3, z3, P.z, wep_rat_function::modulus);

	square(y3, t5, wep_rat_function::modulus);
	multiply(t4, t1, y3, wep_rat_function::modulus);
	square(x3, t3, wep_rat_function::modulus);
	mult_by_curve(x3, wep_rat_function::modulus.modulus());
	subtract(x3, x3, t4);

	subtract(t4, t4, 2 * x3);
	multiply(y3, y3, t5, wep_rat_function::modulus);
	multiply(y3, y3, t2, wep_rat_function::modulus);
	multiply(t1, t4, t3, wep_rat_function::modulus);
	subtract(y3, t1, y3);
	multiply_by_scalar(y3, y3, (inverse(bigmod(2))).mantissa());

	R.assign(x3, y3, z3);
}

// computes R = P + (X,Y) cheaper than general point addition

void add_P_XY (wep_rat_function & R, const wep_rat_function & P) 
{
  wep_rat_function::ff_pol t1, t2, t3, t4, t5;
  wep_rat_function::ff_pol x3, y3, z3;

  t1.set_modulus(wep_rat_function::a.modulus());
  t2.set_modulus(wep_rat_function::a.modulus());
  t3.set_modulus(wep_rat_function::a.modulus());
  t4.set_modulus(wep_rat_function::a.modulus());  
  t5.set_modulus(wep_rat_function::a.modulus());  
  x3.set_modulus(wep_rat_function::a.modulus()); 
  y3.set_modulus(wep_rat_function::a.modulus());  
  z3.set_modulus(wep_rat_function::a.modulus());

  square(t4, P.z, wep_rat_function::modulus);
  multiply_by_x_mod(t3, t4, wep_rat_function::modulus.modulus());
  multiply(t4, t4, P.z, wep_rat_function::modulus);
  
  subtract(t5, P.x, t3);  
  add(t1, P.x, t3);
  subtract(t3, P.y, t4);
  add(t2, P.y, t4);

  if (t5.is_zero())
    {
      if (t3.is_zero())
        multiply_by_2(R, P);
      else
        R.assign_zero();
      return;
    }

  multiply(z3, t5, P.z, wep_rat_function::modulus);

  square(y3, t5, wep_rat_function::modulus);
  multiply(t4, t1, y3, wep_rat_function::modulus);
  square(x3, t3, wep_rat_function::modulus);
  mult_by_curve(x3, wep_rat_function::modulus.modulus());
  subtract(x3, x3, t4);

  subtract(t4, t4, 2 * x3);
  multiply(y3, y3, t5, wep_rat_function::modulus);
  multiply(y3, y3, t2, wep_rat_function::modulus);
  multiply(t1, t4, t3, wep_rat_function::modulus);
  subtract(y3, t1, y3);
  multiply_by_scalar(y3, y3, (inverse(bigmod(2))).mantissa());
  
  R.assign(x3, y3, z3);
}


void add (wep_rat_function & R, const wep_rat_function & Q,
          const wep_rat_function & P)
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
			add_PQ(R, P, Q);
	}
}



void subtract (wep_rat_function & R, const wep_rat_function & Q,
               const wep_rat_function & P)
{
	wep_rat_function T(P);
	negate(T, T);
	add(R, Q, T);
}



//--------------------------------------------------------------------
// IEEE P1363: need 3 mults and 6 squarings

void multiply_by_2 (wep_rat_function & Q, const wep_rat_function & P)
{
	if (P.is_zero() || P.y.is_zero()) {
		Q.assign_zero();
		return;
	}

	wep_rat_function::ff_pol m, x3, y3, z3, s, t;

	x3.set_modulus(wep_rat_function::a.modulus());
	y3.set_modulus(wep_rat_function::a.modulus());
	z3.set_modulus(wep_rat_function::a.modulus());
	m.set_modulus(wep_rat_function::a.modulus());
	s.set_modulus(wep_rat_function::a.modulus());
	t.set_modulus(wep_rat_function::a.modulus());

	square(m, P.z, wep_rat_function::modulus);
	square(m, m, wep_rat_function::modulus);
	multiply_by_scalar(m, m, wep_rat_function::a.mantissa());
	square(s, P.x, wep_rat_function::modulus);
	multiply_by_scalar(s, s, 3);
	add(m, m, s);

	multiply(z3, P.y, P.z, wep_rat_function::modulus);
	multiply_by_scalar(z3, z3, 2);

	square(t, P.y, wep_rat_function::modulus);
	mult_by_curve(t, wep_rat_function::modulus.modulus());
	multiply(s, t, P.x, wep_rat_function::modulus);
	multiply_by_scalar(s, s, 4);

	square(x3, m, wep_rat_function::modulus);
	subtract(x3, x3, 2 * s);

	square(t, t, wep_rat_function::modulus);
	multiply_by_scalar(t, t, 8);

	subtract(y3, s, x3);
	multiply(y3, m, y3, wep_rat_function::modulus);
	subtract(y3, y3, t);

	mult_by_curve(x3, wep_rat_function::modulus.modulus());
	mult_by_curve(y3, wep_rat_function::modulus.modulus());
	mult_by_curve(z3, wep_rat_function::modulus.modulus());

	Q.assign(x3, y3, z3);
}



//---------------------------------------------------------
//  The following function implements a left-to-right fast
//  exponentiation which bits 0, +1, -1. See the paper "Fast
//  multiplication on elliptic curves" by Volker Mueller

void multiply(wep_rat_function & Q, lidia_size_t n,
	      const wep_rat_function & PP)
{
	wep_rat_function P(PP);
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



	

// R1 = n1 * P, R2 = n2 * P, n1, n2 >= 0

void multiply (wep_rat_function & R1, wep_rat_function & R2,
               lidia_size_t n1, lidia_size_t n2,
               const wep_rat_function & PP)
{

	wep_rat_function P (PP);
	R1.assign_zero(); R2.assign_zero();

	if (n1 < 0 || n2 < 0)
		lidia_error_handler("wep_rat_function", "multiply::either n1 or n2 "
				    "is negative");

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

wep_rat_function &
wep_rat_function::operator = (const wep_rat_function & wep)
{
	if (this != & wep)
		this->assign (wep.x, wep.y, wep.z);
	return *this;
}



void wep_rat_function::assign_xy ()
{
	z.assign_one();
	y.assign_one();
	x.assign_x();
}



							
void wep_rat_function::
assign (const ff_pol & xx, const ff_pol & yy, const ff_pol & zz)
{
	x.assign(xx);
	y.assign(yy);
	z.assign(zz);
}



void wep_rat_function::
assign (const ff_pol & xx, const ff_pol & yy)
{
	x.assign(xx);
	y.assign(yy);
	z.assign_one();
}



void wep_rat_function::
assign (const wep_rat_function & w)
{
	x.assign(w.x);
	y.assign(w.y);
	z.assign(w.z);
}



void wep_rat_function::
initialize(const wep_rat_function::ff_element & aa,
           const wep_rat_function::ff_element & bb,
           const wep_rat_function::ff_pol_m & pm)
{
	a.assign(aa); b.assign(bb);
	modulus = pm;
}



std::ostream & operator << (std::ostream & o, const wep_rat_function & P)
{
	o << "(" << P.x << " : " << P.y << " : " << P.z << ")";
	return o;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

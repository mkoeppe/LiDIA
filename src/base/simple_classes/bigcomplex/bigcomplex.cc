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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigcomplex.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void
bigcomplex::invert ()
{
	if (re.is_zero() && im.is_zero()) {
		lidia_error_handler("bigcomplex", "invert::division by zero.");
		return;
	}

	bigfloat inv_den, tmp1;

	square(tmp1, re);
	square(inv_den, im);
	add(tmp1, tmp1, inv_den);
	divide(inv_den, 1, tmp1);

	multiply(re, re, inv_den);
	im.negate();
	multiply(im, im, inv_den);
}



//
// procedural arithmetic
//

// (cr, ci) = (xr, xi) * (yr, yi)
//          = (xr*yr - xi*yi, xr*yi + xi*yr)
void
multiply (bigcomplex &c, const bigcomplex &x, const bigcomplex &y)
{
	bigfloat tmp, tmpr, tmpi;

	multiply(tmpr, x.re, y.re);
	multiply(tmpi, x.im, y.im);
	subtract(tmp, tmpr, tmpi);

	multiply(tmpr, x.re, y.im);
	multiply(tmpi, x.im, y.re);
	add(c.im, tmpr, tmpi);
	c.re.assign(tmp);
}



// from romine@xagsun.epm.ornl.gov
// (cr, ci) = (xr, xi) / (yr, yi)
void
divide (bigcomplex &c, const bigcomplex &x, const bigcomplex &y)
{
	bigfloat den = abs(y.re) + abs(y.im);
	if (den.is_zero()) {
		lidia_error_handler("bigcomplex", "operator/::division by zero.");
		return;
	}

	bigfloat inv_den, xrden, xiden, yrden, yiden;;
	divide(inv_den, 1, den);

	multiply(xrden, x.re, inv_den);
	multiply(xiden, x.im, inv_den);
	multiply(yrden, y.re, inv_den);
	multiply(yiden, y.im, inv_den);

	bigfloat inv_nrm, tmp1, tmp2;

	square(tmp1, yrden);
	square(tmp2, yiden);
	add(tmp1, tmp1, tmp2);
	divide(inv_nrm, 1, tmp1);

	multiply(tmp1, xrden, yrden);
	multiply(tmp2, xiden, yiden);
	add(tmp1, tmp1, tmp2);
	multiply(c.re, tmp1, inv_nrm);

	multiply(tmp1, xiden, yrden);
	multiply(tmp2, xrden, yiden);
	subtract(tmp1, tmp1, tmp2);
	multiply(c.im, tmp1, inv_nrm);

}



void
divide (bigcomplex &c, const bigfloat &x, const bigcomplex &y)
{
	if (y.re.is_zero() && y.im.is_zero()) {
		lidia_error_handler("bigcomplex", "operator/2::division by zero.");
		return;
	}

	bigfloat inv_den, tmp1, tmp2;

	square(tmp1, y.re);
	square(tmp2, y.im);
	add(tmp1, tmp1, tmp2);
	divide(inv_den, 1, tmp1);

	multiply(c.re, y.re, x);
	multiply(c.re, c.re, inv_den);

	multiply(c.im, y.im, x);
	multiply(c.im, c.im, inv_den);
	c.im.negate();
}



//
// functions
//

void
sqrt (bigcomplex & c, const bigcomplex & x)
{
	if (x.re.is_zero() && x.im.is_zero()) {
		c.re.assign_zero();
		c.im.assign_zero();
	}
	else {
		bigfloat s, t, d;

		s = hypot(x);
		t.assign(x.re);
		t.absolute_value();
		add(t, t, s);
		t.divide_by_2();
		sqrt(s, t);

		divide(d, x.im, s);
		d.divide_by_2();

		if (x.re.is_gt_zero()) {
			c.re.assign(s);
			c.im.assign(d);
		}
		else if (x.im.is_ge_zero()) {
			c.re.assign(d);
			c.im.assign(s);
		}
		else {
			d.negate();
			s.negate();
			c.re.assign(d);
			c.im.assign(s);
		}
	}
}



// Corrections based on reports from: thc@cs.brown.edu&saito@sdr.slb.com
void
power (bigcomplex & c, const bigcomplex & x, const bigcomplex & p)
{
	bigfloat h = hypot(x);
	if (h.is_le_zero()) {
		lidia_error_handler("bigcomplex", "power of number less equal zero");
		c.assign_zero();
		return;
	}
	bigfloat a, lr, li;
	atan2(a, x.im, x.re);
	power(lr, h, p.re);
	multiply(li, p.re, a);
	if (!p.im.is_zero()) {
		bigfloat tmp;
		multiply(tmp, p.im, a);
		exp(tmp, tmp);
		divide(lr, lr, tmp);
		log(tmp, h);
		multiply(tmp, tmp, p.im);
		add(li, li, tmp);
	}
	cos(c.re, li);
	multiply(c.re, c.re, lr);
	sin(c.im, li);
	multiply(c.im, c.im, lr);
}



void
power (bigcomplex & c, const bigcomplex & x, const bigfloat & p)
{
	if (p.is_zero()) {
		c.re.assign_one();
		c.im.assign_zero();
		return;
	}
	if (x.is_zero()) {
		c.re.assign_zero();
		c.im.assign_zero();
		return;
	}

	long pp = 0;
	p.longify(pp);
	if (p == pp) {
		// treat exponents which fit into a 'long' differently...
		bigcomplex b(x);
		c = bigcomplex(1.0, 0.0);
		if (pp < 0) {
			pp = -pp;
			b.invert(); // b = 1.0 / b
		}
		for (;;) {
			if (pp & 1)
				c *= b;
			if ((pp >>= 1) == 0)
				break;
			else
				b *= b;
		}
	}
	else {
		// exponents which do not fit into a 'long'
		bigfloat h = hypot(x);
		if (h.is_le_zero()) {
			lidia_error_handler("bigcomplex", "power of number less equal zero");
			c.assign_zero();
			return;
		}
		bigfloat lr, a;

		power(lr, h, p);
		atan2(a, x.im, x.re);

		bigfloat li;

		multiply(li, p, a);
		cos(c.re, li);
		multiply(c.re, c.re, lr);
		sin(c.im, li);
		multiply(c.im, c.im, lr);
	}
}



void
square (bigcomplex &c, const bigcomplex &y)
{
	bigfloat re2, im2, reim;

	square(re2, y.re);
	square(im2, y.im);
	multiply(reim, y.re, y.im);
	subtract(c.re, re2, im2);
	c.im.assign(reim);
	c.im.multiply_by_2();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

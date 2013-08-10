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
//	Author	:Volker Mueller (VM), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/ff2.h"
#include	"LiDIA/random_generator.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



ff1 ff2::tau;



ff2::ff2 ()
{
	debug_handler ("ff2", "ff2()");
}



ff2::ff2 (const ff2 & y)
{
	debug_handler ("ff2", "ff2(const ff2&)");
	a = y.a;
	b = y.b;
}



ff2::ff2 (const ff1 & y1, const ff1 & y2)
{
	debug_handler ("ff2", "ff2(const ff1&, const ff1&)");
	a = y1;
	b = y2;
}



ff2::ff2 (const ff1 & y1)
{
	debug_handler ("ff2", "ff2(const unsigned long)");
	a = y1;
	b.assign_zero();
}



ff2::~ff2 ()
{
	debug_handler ("ff2", "~ff2()");
}



void
ff2::init_field (const udigit & characteristic)
{
	debug_handler ("ff2", "init_field");

	ff1::set_characteristic (characteristic);
	(ff2::tau).assign_non_square();
}



udigit
ff2::characteristic ()
{
	debug_handler ("ff2", "characteristic()");
	return(ff1::characteristic());
}



ff2 & ff2::operator = (const ff2 & y)
{
	debug_handler ("ff2", "operator=(const ff2&)");
	a = y.a;
	b = y.b;
	return (*this);
}



void
ff2::assign (const ff2 & y)
{
	debug_handler ("ff2", "assign(const ff2&)");
	a = y.a;
	b = y.b;
}



void
ff2::assign(const ff1 & y1, const ff1 & y2)
{
	debug_handler ("ff2", "assign(const ff1&, const ff1&)");
	a = y1;
	b = y2;
}



void
ff2::assign_one ()
{
	debug_handler ("ff2", "assign_one()");
	a.assign_one  ();
	b.assign_zero ();
}



void
ff2::assign_zero()
{
	debug_handler ("ff2", "assign_zero()");
	a.assign_zero();
	b.assign_zero();
}



bool
ff2::is_one () const
{
	debug_handler ("ff2", "is_one()");
	return (a.is_one() && b.is_zero());
}



bool
ff2::is_zero () const
{
	debug_handler ("ff2", "is_zero()");
	return (a.is_zero() && b.is_zero());
}



std::istream &
operator >> (std::istream & in, ff2 & y)
{
	debug_handler ("ff2", "operator >> (istream&, ff2&)");
	in >> y.a;
	in >> y.b;
	return (in);
}



std::ostream &
operator << (std::ostream & out, const ff2 & y)
{
	debug_handler ("ff2", "operator >> (ostream&, const ff2&)");
	out << "[" << y.b;
	out << " * rho + ";
	out << y.a << ", rho^2 = " << y.tau << "]";
	return (out);
}



void
add (ff2 & c, const ff2 & a, const ff2 & b)
{
	debug_handler ("ff2", "add(ff2&, const ff2&, const ff2&)");
	add (c.a, a.a, b.a);
	add (c.b, a.b, b.b);
}



void
subtract (ff2 & c, const ff2 & a, const ff2 & b)
{
	debug_handler ("ff2", "subtract(ff2&, const ff2&, const ff2&)");
	subtract (c.a, a.a, b.a);
	subtract (c.b, a.b, b.b);
}



void
multiply (ff2 & c, const ff2 & a, const ff2 & b)
{
	debug_handler ("ff2", "multiply(ff2&, const ff2&, const ff2&)");

	ff1 h, r, h2;

	multiply (r, a.a, b.a);
	multiply (h, a.b, b.b);
	multiply (h, h, ff2::tau);
	add      (r, r, h);

	multiply (h, a.a, b.b);
	multiply (h2, a.b, b.a);
	add      (h, h2, h);
	c.assign(r, h);
}



void
divide (ff2 & c, const ff2 & a, const ff2 & b)
{
	debug_handler ("ff2", "divide(ff2&, const ff2&, const ff2&)");
	ff2 h;
	invert(h, b);
	multiply(c, a, h);
}



void
invert (ff2 & c, const ff2 & a)
{
	debug_handler ("ff2", "invert(ff2&, const ff2&)");

	if (a.is_zero()) {
		lidia_error_handler ("ff2::invert(ff2&, const ff2&)",
				     "Division by zero.");
	}
	else {
		ff1 h1, h2;

		square   (h1, a.a); // inverse is 1/Norm * (a - b*X)
		square   (h2, a.b);
		multiply (h2, h2, ff2::tau);
		subtract (h1, h1, h2);
		invert   (h2, h1);

		multiply (c.a, h2, a.a);
		negate   (h1, h2);
		multiply (c.b, h1, a.b);
	}
}



void
negate (ff2 & c, const ff2 & a)
{
	debug_handler ("ff2", "negate(ff2&, const ff2&)");
	negate (c.a, a.a);
	negate (c.b, a.b);
}



void
power (ff2 & c, const ff2 & a, long k)
{
	debug_handler ("ff2", "power(ff2&, const ff2&, long)");

	ff2 b;

	if (k == 0) {
		c.assign_one();
	}
	else {
		if (k < 0) {
			k = -k;
			negate(b, a);
		}
		else
			b = a;
		c.assign_one();

		while (k >= 1) {
			if (k & 1)
				multiply (c, c, b);
			square (b, b);
			k >>= 1;
		}
	}
}



//------------------------------------------------------------------
// determine all pairs (x, x^{-1}) of order d, set w to all such x
// i.e. store just one representative for (x, x^{-1})

void
nearly_all_of_order (udigit d, base_vector< ff2 > & w)
{
	debug_handler("ff2", "nearly_all_of_order(udigit, base_vector< ff2 >&)");

	ff2 generator;
	ff2 *elems;
	int i, num = 1, order;
	random_generator rg;
	udigit u1, u2, p;

	p = ff2::characteristic();

	do {
		rg >> u1; rg >> u2;
		generator.assign(ff1(u1 % p), ff1(u2 % p));
		if (generator.is_zero())
			order = 1;
		else order = generator.multiplicative_order();
	} while (order % d != 0);

	elems = new ff2[d];
	power (elems[0], generator, order/d);

	for (i = 2; i <= static_cast<int>(d/2); i++) {
		if (gcd(i, d) == 1)
			power(elems[num++], elems[0], i);
	}

	w.set_capacity (num);
	for (i = 0; i < num; i++)
		w[i] = elems[i];
	delete[] elems;
}



udigit
ff2::multiplicative_order () const
{
	debug_handler ("ff2", "multiplicative_order()");

	if (is_zero()) {
		lidia_error_handler ("ff2::multiplicative_order",
				     "Zero is not an element of (F_p^2)^*.");
		return 0;
	}

	udigit exp = 1;
	ff2    elem = *this;

	while (!elem.is_one()) {
		exp++;
		multiply (elem, elem, *this);
	}

	return exp;
}



void swap (ff2 & aa, ff2 & bb)
{
	swap (aa.a, bb.a);
	swap (aa.b, bb.b);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

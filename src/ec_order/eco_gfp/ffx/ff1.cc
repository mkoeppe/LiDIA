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
#include	"LiDIA/ff1.h"
#include	"LiDIA/sort_vector.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



udigit ff1::p;



// ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

ff1::ff1 ()
{
	debug_handler ("ff1", "ff1()");
}



ff1::ff1 (const ff1 & y)
{
	debug_handler ("ff1", "ff1(const ff1&)");
	x = y.x;
}



ff1::ff1 (const udigit & y)
{
	debug_handler ("ff1", "ff1(const ff1&)");
	x = y % p;
}



ff1::~ff1 ()
{
	debug_handler ("ff1", "~ff1()");
}



void
ff1::set_characteristic (udigit characteristic)
{
	debug_handler ("ff1", "set_characteristic(udigit)");

	if (characteristic <= 2)
		lidia_error_handler ("ff1::set_characteristic(udigit)",
				     "Characteristic is less than 3; not implemented.");

	p = characteristic;
}



ff1 &
ff1::operator = (const ff1 & y)
{
	debug_handler ("ff1", "operator=(const ff1&)");
	x = y.x;
	return (*this);
}



bool
ff1::operator == (const ff1 & y) const
{
	debug_handler ("ff1", "operator==(const ff1&) const");
	return (x == y.x);
}



void
ff1::assign_one ()
{
	debug_handler ("ff1", "assign_one()");
	x = 1;
}



void
ff1::assign_zero ()
{
	debug_handler ("ff1", "assign_zero()");
	x = 0;
}



void
ff1::assign_non_square ()
{
	debug_handler ("ff1", "assign_non_square()");
	x = 2;
	while (jacobi(x, p) == 1)
		x = next_prime(x+1);
}



void
ff1::assign(udigit a)
{
	debug_handler ("ff1", "assign(udigit)");
	if (a > p)
		x = divide(a, 0, a, p);
	else
		x = a;
}



bool
ff1::is_one () const
{
	debug_handler ("ff1", "is_one() const");
	return (x == 1);
}



bool
ff1::is_zero () const
{
	debug_handler ("ff1", "is_zero() const");
	return (x == 0);
}



std::istream &
operator >> (std::istream & in, ff1 & y)
{
	debug_handler ("ff1", "operator>>(istream&, ff1&)");
	in >> y.x;
	return (in);
}



std::ostream &
operator << (std::ostream & out, const ff1 & y)
{
	debug_handler ("ff1", "operator>>(ostream&, const ff1&)");
	out << y.x;
	return (out);
}



void
add (ff1 & c, const ff1 & a, const ff1 & b)
{
	debug_handler ("ff1", "add(ff1&, const ff1&, const ff1&)");
	c.x = add_mod (a.x, b.x, ff1::p);
}



void
subtract (ff1 & c, const ff1 & a, const ff1 & b)
{
	debug_handler ("ff1", "subtract(ff1&, const ff1&, const ff1&)");
	c.x = subtract_mod (a.x, b.x, ff1::p);
}



void
multiply (ff1 & c, const ff1 & a, const ff1 & b)
{
	debug_handler ("ff1", "multiply(ff1&, const ff1&, const ff1&)");
	c.x = multiply_mod(a.x, b.x, ff1::p);
}



void
divide (ff1 & c, const ff1 & a, const ff1 & b)
{
	debug_handler ("ff1", "divide(ff1&, const ff1&, const ff1&)");
	c.x = divide_mod (a.x, b.x, ff1::p);
}



void
invert (ff1 & c, const ff1 & a)
{
	debug_handler ("ff1", "invert(ff1&, const ff1&)");
	c.x = invert_mod (a.x, ff1::p);
}



void
negate (ff1 & c, const ff1 & a)
{
	debug_handler ("ff1", "negate(ff1&, const ff1&)");
	c.x = negate_mod (a.x, ff1::p);
}



void
power (ff1 & c, const ff1 & a, long e)
{
	debug_handler ("ff1", "power(ff1&, const ff1&, int)");
	if (e < 0) {
		invert(c, a);
		e = -e;
	}
	else
		c = a;
	c.x = power_mod (c.x, static_cast<udigit>(e), ff1::p);
}



//------------------------------------------------------------------
// determine all pairs (x, x^{-1}) of oder d, set w to all such x
// i.e. store just one representative for (x, x^{-1})

void
nearly_all_of_order (udigit d, base_vector< ff1 > & w)
{
	debug_handler ("ff1", "nearly_all_of_order(udigit, base_vector< ff1 >&)");
	ff1 generator;
	ff1 *elems;
	int i, num = 1, order;

	generator.x = 2;
	order = generator.multiplicative_order();

	while (order % d != 0) {
		generator.x = next_prime(generator.x+1);
		order = generator.multiplicative_order();
	}

	elems = new ff1[d];
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
ff1::multiplicative_order () const   // this can be done better
{
	debug_handler ("ff1", "multiplicative_order() const");
	if (x == 0) {
		lidia_error_handler ("ff1::multiplicative_order",
				     "Zero is not an element of Fp^*.");
		return 0;
	}
	else {
		udigit exp = 1;
		ff1  elem = *this;

		while (! (elem.is_one())) {
			exp++; multiply(elem, elem, *this);
		}
		return exp;
	}
}



bool
ff1::is_square () const
{
	debug_handler ("ff1", "is_square() const");

	if (x == 0)
		return true;
	else
		if (jacobi(x, ff1::p) != -1)
			return true;
		else
			return false;
}



void
sqrt (ff1 & R, const ff1 & a)    // It is assumed that a is a square.
{
	debug_handler ("ff1", "sqrt(ff1, const ff1&)");

	if (a.x == 0)
		R.x = 0;
	else
		R.x = sqrt_mod(a.x, ff1::p);
}



void swap (ff1 & a, ff1 & b)
{
	ff1 h;
	h = a; a = b; b = h;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

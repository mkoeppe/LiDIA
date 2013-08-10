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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf2n.h"
#include	"LiDIA/nmbrthry_functions.h"
#include	"LiDIA/random_generator.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// STATIC VARIABLES
// ----------------------------------------------------------------------

unsigned int gf2n::anzBI;
unsigned int gf2n::degree = 0;
rational_factorization gf2n::ord;
bool gf2n::ord_is_fact = false;
unsigned int *gf2n::exponents;
unsigned int gf2n::anz_exponents;
unsigned int gf2n::mulsel;
unsigned int gf2n::invsel;
gf2n_word *gf2n::B, *gf2n::C, *gf2n::F, *gf2n::G;



//
// CONSTRUCTORS and DESTRUCTOR
// ----------------------------------------------------------------------

gf2n::gf2n ()
{
	register unsigned int i;

	element = new gf2n_word[gf2n::anzBI];
	if (element == NULL) {
		lidia_error_handler ("gf2n", "gf2n()::not enough memory \n");
		return;
	}

	for (i = 0; i < gf2n::anzBI; i++)
		element[i] = static_cast<gf2n_word>(0);
}



gf2n::gf2n (const gf2n & a)
{
	register unsigned int i;
	element = new gf2n_word[gf2n::anzBI];
	if (element == NULL) {
		lidia_error_handler ("gf2n", "gf2n(const gf2n &)::not enough memory \n");
		return;
	}

	for (i = 0; i < gf2n::anzBI; i++)
		element[i] = a.element[i];
}



gf2n::gf2n (unsigned long ui)
{
	register unsigned int i;

	if (degree < (8*sizeof(gf2n_word)))
		if ((static_cast<gf2n_word>(1) << degree) <= static_cast<gf2n_word>(ui)) {
			lidia_error_handler ("gf2n", "gf2n(unsigned long)::input out of range ");
			return;
		}

	element = new gf2n_word[gf2n::anzBI];

	if (element == NULL) {
		lidia_error_handler ("gf2n", "gf2n(unsigned long)::not enough memory \n");
		return;
	}

	for (i = 1; i < gf2n::anzBI; i++)
		element[i] = static_cast<gf2n_word>(0);

	element[0] = static_cast<gf2n_word>(ui);
}



gf2n::gf2n (const bigint & bi)
{
	register unsigned int i;
	bigint h;

	shift_left(h, bigint(1), degree);
	if (h <= bi || bi.is_negative()) {
		lidia_error_handler ("gf2n", "gf2n(const bigint &)::input out of range ");
		return;
	}

	element = new gf2n_word[gf2n::anzBI];

	if (element == NULL) {
		lidia_error_handler ("gf2n", "gf2n(const bigint &)::not enough memory \n");
		return;
	}

	h.assign(bi);
	i = 0;

	while (!h.is_zero()) {
		element[i++] = static_cast<gf2n_word>(h.least_significant_digit());
		shift_right(h, h, 8*sizeof(gf2n_word));
	}

	while (i < gf2n::anzBI)
		element[i++] = static_cast<gf2n_word>(0);
}



gf2n::~gf2n ()
{
	delete[] element;
}



//
// NEW INIT AFTER CHANGE OF FIELD
// ----------------------------------------------------------------------

void
gf2n::re_initialize()
{
	register unsigned int i;

	delete[] element;
	element = new gf2n_word[gf2n::anzBI];
	if (element == NULL) {
		lidia_error_handler ("gf2n", "re_initialize::not enough memory \n");
		return;
	}

	for (i = 0; i < gf2n::anzBI; i++)
		element[i] = static_cast<gf2n_word>(0);

	if (gf2n::table_solve_quadratic_init == true)  // delete the tables
		gf2n::delete_table_for_solve_quadratic();
}



//
// ASSIGNMENTS
// ----------------------------------------------------------------------

gf2n & gf2n::operator = (const gf2n & a)
{
	for (register unsigned int i = 0; i < gf2n::anzBI; i++)
		element[i] = a.element[i];
	return (*this);
}



gf2n & gf2n::operator = (unsigned long ui)
{
	if (degree < 8*sizeof(unsigned long))
		if ((unsigned long) (1 << degree) <= ui) {
			lidia_error_handler ("gf2n", "operator = (unsigned long)::input out of range");
			return *this;
		}


	element[0] = static_cast<gf2n_word>(ui);
	for (register unsigned int i = 1; i < gf2n::anzBI; i++)
		element[i] = static_cast<gf2n_word>(0);

	return (*this);
}



gf2n & gf2n::operator = (const bigint & bi)
{
	register unsigned int i = 0;
	bigint h;

	shift_left(h, bigint(1), degree);
	if (h <= bi || bi.is_negative()) {
		lidia_error_handler ("gf2n", "operator = (const bigint &)::input out of range ");
		return *this;
	}

	h.assign(bi);

	while (!h.is_zero()) {
		element[i++] = static_cast<gf2n_word>(h.least_significant_digit());
		shift_right(h, h, 8*sizeof(gf2n_word));
	}

	while (i < gf2n::anzBI)
		element[i++] = static_cast<gf2n_word>(0);

	return (*this);
}



void gf2n::assign_zero()
{
	for (register unsigned int i = 0; i < gf2n::anzBI; i++)
		gf2n::element[i] = static_cast<gf2n_word>(0);
}



void gf2n::assign_one()
{
	for (register unsigned int i = 1; i < gf2n::anzBI; i++)
		gf2n::element[i] = static_cast<gf2n_word>(0);

	gf2n::element[0] = static_cast<gf2n_word>(1);
}



void gf2n::assign(const gf2n &a)
{
	for (register unsigned int i = 0; i < gf2n::anzBI; i++)
		gf2n::element[i] = a.element[i];
}



void swap(gf2n & a, gf2n & b)
{
	gf2n c(a);

	a.assign(b);
	b.assign(c);
}



//
// COMPARISONS
// ----------------------------------------------------------------------

bool operator == (const gf2n & a, const gf2n & b)
{
	register int i = gf2n::anzBI - 1;

	if (a.is_reduced() == false)
		partial_reduce2[gf2n::invsel](a.element);

	if (b.is_reduced() == false)
		partial_reduce2[gf2n::invsel](b.element);

	while (i >= 0 && a.element[i] == b.element[i])
		i--;

	if (i < 0)
		return true;
	else
		return false;
}



bool operator != (const gf2n & a, const gf2n & b)
{
	return (!(a == b));
}



bool gf2n::is_zero () const
{
	int i = gf2n::anzBI - 1;

	if ((*this).is_reduced() == false)
		partial_reduce2[gf2n::invsel](element);

	while (i >= 0 && element[i] == static_cast<gf2n_word>(0))
		i--;

	if (i < 0)
		return true;
	else
		return false;
}



bool gf2n::is_one ()  const
{
	int i = gf2n::anzBI - 1;

	if ((*this).is_reduced() == false)
		partial_reduce2[gf2n::invsel](element);

	while (i >= 0 && element[i] == static_cast<gf2n_word>(0))
		i--;

	if ((i > 0) || (i < 0))
		return false;
	else
		if (element[i] == static_cast<gf2n_word>(1))
			return true;
		else
			return false;
}



// ======================================================================
// ARITHMETIC                                        
// ======================================================================

gf2n operator * (const gf2n & a, const gf2n & b)
{
	gf2n c;
	multiply (c, a, b);
	return (c);
}



gf2n operator + (const gf2n & a, const gf2n & b)
{
	gf2n c;
	add (c, a, b);
	return (c);
}



gf2n operator - (const gf2n & a, const gf2n & b)
{
	gf2n c;
	subtract (c, a, b);
	return (c);
}



gf2n operator / (const gf2n & a, const gf2n & b)
{
	gf2n d(b);
	d.invert();
	multiply (d, a, d);
	return (d);
}



gf2n & gf2n::operator += (const gf2n & y)
{
	add (*this, *this, y);
	return *this;
}



gf2n & gf2n::operator -= (const gf2n & y)
{
	add (*this, *this, y);
	return *this;
}



gf2n & gf2n::operator *= (const gf2n & y)
{
	multiply (*this, *this, y);
	return *this;
}



gf2n & gf2n::operator /= (const gf2n & y)
{
	divide (*this, *this, y);
	return *this;
}



extern void (*uinvert[])(gf2n_word*, gf2n_word*);



//
// INVERSE
// ----------------------------------------------------------------------

gf2n inverse (const gf2n & x)
{
#ifdef __EMX__
	gf2n erg(x);
	erg.invert();
#else
	gf2n erg;

	if (x.is_one () == true)
		erg.assign_one ();
	else if (x.is_zero () == true) {
		lidia_error_handler ("gf2n", "inverse(const gf2n &)::can't invert 0 \n");
		return (gf2n(0));
	}
	else
		uinvert[gf2n::invsel] (erg.element, x.element);
#endif

	return erg;
}



//
// INVERT II
// ----------------------------------------------------------------------

void gf2n::invert ()
{
	if ((*this).is_one () == true)
		return;

	if ((*this).is_zero() == true) {
		lidia_error_handler ("gf2n", "invert()::can't invert 0 \n");
		return;
	}

	gf2n erg(*this);
	uinvert[gf2n::invsel](element, erg.element);
}



//
// SQRT
// ----------------------------------------------------------------------

void sqrt (gf2n & c, const gf2n & a)
{
	c.assign(a);

	for (register unsigned int i = 1; i < gf2n::degree; i++)
		square (c, c);
}



gf2n sqrt (const gf2n & a)
{
	gf2n c;
	sqrt(c, a);
	return c;
}



//
// RANDOM
// ----------------------------------------------------------------------

void gf2n::randomize(unsigned int d)
{
	static random_generator rg;
	long r;
	register unsigned int i;

	for (i = 0; i < anzBI; i++) {
		rg >> r;
		element[i] = static_cast<gf2n_word>(r);
	}
	partial_reduce2[gf2n::invsel](gf2n::element);

	if (d == degree)
		return;

	if (degree % d != 0) {
		lidia_error_handler("gf2n",
				    "randomize(unsigned int d)::d is no divisor of "
				    "field degree");
		return;
	}

	bigint e1, e2;

	shift_left(e1, bigint(1), gf2n::degree);
	dec(e1);
	shift_left(e2, bigint(1), d);
	dec(e2);
	divide(e1, e1, e2); // (*this)^e1 is an element of GF(2^d)
	power(*this, *this, e1);
}



gf2n randomize (const gf2n &x, unsigned int d)
{
	gf2n a(x);
	a.randomize(d);
	return a;
}



//**************************************************************************

bigint compute_order (const gf2n & a)
{
	bigint q, h, e_ord;
	unsigned int  j, i;
	gf2n res;

	if (a.is_zero()) {
		lidia_error_handler("gf2n", "compute_order(const gf2n&)::zero is no element of multiplicative group");
		return (0);
	}

	if (a.is_one())
		return bigint(1);

	if (gf2n::ord_is_fact == false) {
		// does this need a 64 bit adaption ?
		if (gf2n::degree > 30 && !(gf2n::degree & 1)) {
			i = gf2n::degree;

			while (!(i & 1) && i > 20)
				i >>= 1;

			shift_left(q, bigint(1), i);
			multiply(gf2n::ord, factor_completely(q-bigint(1)),
				 factor_completely(q+bigint(1)));
			i <<= 1;
			while (i < gf2n::degree) {
				square(q, q);
				multiply(gf2n::ord, gf2n::ord, factor_completely(q+bigint(1)));
				i <<= 1;
			}
		}
		else {
			shift_left(q, bigint(1), gf2n::degree);
			dec(q);
			gf2n::ord = factor_completely (q);
		}

		if (gf2n::ord.is_prime_factorization() == false) {
			lidia_error_handler("gf2n", "compute_order::can't factor group order");
			return (0);
		}

		gf2n::ord_is_fact = true;
	}

	// now compute order by dividing (q-1) by all prime factors of group order

	shift_left(q, bigint(1), gf2n::degree); // q = group order
	dec(q);
	e_ord.assign_one();

	for (i = 0; i < static_cast<unsigned int>(gf2n::ord.no_of_comp()); i++) {
		power(h, gf2n::ord.base(i), gf2n::ord.exponent(i));
		divide(h, q, h);

		power(res, a, h);
		j = 0;

		while (j < static_cast<unsigned int>(gf2n::ord.exponent(i)) && !res.is_one()) {
			j ++;
			multiply(e_ord, e_ord, gf2n::ord.base(i));
			power(res, res, gf2n::ord.base(i));
		}
		if (!res.is_one()) {
			lidia_error_handler("gf2n", "compute_order()::element order does not divide (q-1) !! \n");
			return 0;
		}
	}

	return (e_ord);
}



//
// gf2n - HASH
//

gf2n_word hash (const gf2n & a)
{
	if (!a.is_reduced())
		partial_reduce2[gf2n::invsel](a.element);

	return static_cast<gf2n_word>(a.element[0]);
}



unsigned int gf2n::relative_degree() const
{
	gf2n h;
	unsigned int i;

	if (is_one() || is_zero())
		return 1;

	square(h, *this);

	for (i = 1; i <= gf2n::degree; i++) {
		if (gf2n::degree % i == 0 && h == *this)   // *this in GF(2^i)
			return i;
		square(h, h);
	}
	lidia_error_handler("gf2n", "relative_degree::found no subfield");
	return 0;
}



//
// GENERATOR
// ----------------------------------------------------------------------

gf2n get_generator (unsigned int d)
{
	gf2n a, apow, ahelp;
	bigint *v, h1, gord, h;
	bool found = false;
	unsigned int i, length;

	if (gf2n::degree % d != 0) {
		lidia_error_handler("gf2n", "get_generator::d does not divide degree");
		return 1;
	}

	shift_left(gord, bigint(1), gf2n::degree); // order of mult. group
	dec(gord);


	if (gf2n::ord_is_fact == false)              // factor order, if necessary
	{
		gf2n::ord = factor_completely(gord);
		gf2n::ord_is_fact = true;
	}

	if (is_prime(gord, 8) && d == gf2n::degree) {
		do {
			a.randomize();
		} while (a.is_one() || a.is_zero());
		return a;
	}

	length = gf2n::ord.no_of_comp();
	v = new bigint[length];

	// compute in v the difference between ord/p_i for prime factors p_i
	// of ord

	divide(v[0], gord, gf2n::ord.base(length-1));
	h.assign(v[0]);
	i = 1;

	while (i < length) {
		divide(h1, gord, gf2n::ord.base(length-1-i));
		subtract(v[i++], h1, h);
		h.assign(h1);
	}


	while (found == false) {
		// check for generator 
		do {
			a.randomize();
		} while (a.is_zero() || a.is_one());

		power(apow, a, v[0]);
		if (apow.is_one())
			continue;

		i = 1;
		while (i < length) {
			power(ahelp, a, v[i++]);
			multiply(apow, apow, ahelp);
			if (apow.is_one())
				break;
		}
		if ((i == length) && (apow.is_one() == false))
			found = true;
	}
	delete[] v;

	if (d < gf2n::degree) {
		shift_left(h, bigint(1), d);
		dec(h);
		divide(h, gord, h);
		power(a, a, h);
	}

	return a;
}



int gf2n::trace() const
{
	gf2n a(*this), h(*this);
	int i, d;

	d = gf2n::get_absolute_degree();
	for (i = 1; i < d; i++) {
		square(a, a);
		add(h, h, a);
	}

	if (h.is_one())
		return 1;
	if (h.is_zero())
		return 0;
	else
		lidia_error_handler("gf2n", "trace::result not in GF(2)");
	return -1;
}



void power(gf2n & c, const gf2n & a, const long x)
{
	if (x < 0)
		lidia_power(c, inverse(a), static_cast<unsigned long>(-x));
	else
		lidia_power(c, a, static_cast<unsigned long>(x));
}



void power(gf2n & c, const gf2n & a, const bigint & x)
{
	if (x < 0)
		lidia_power(c, inverse(a), -x);
	else
		lidia_power(c, a, x);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

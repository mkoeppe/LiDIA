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
//	Author	: Nigel Smart (NiSm), Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigfloat_int.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void
bigfloat_int::normalize (bigfloat& Lo1, bigfloat &Lo2, bigfloat& Lo3, bigfloat &Lo4)
{
	Low = Lo1;
	Upp = Lo1;

	if (Lo1.compare(Lo2) == 1) {
		Low = Lo2;
	}
	else {
		Upp = Lo2;
	}
	if (Lo3.compare(Low) == -1)
		Low = Lo3;
	if (Lo3.compare(Upp) == 1)
		Upp = Lo3;
	if (Lo4.compare(Low) == -1)
		Low = Lo4;
	if (Lo4.compare(Upp) == 1)
		Upp = Lo4;
}



int
bigfloat_int::sign () const
{
	int s = 0;
	if (Upp.is_lt_zero()) {
		s = -1;
	}
	if (Low.is_gt_zero()) {
		s = 1;
	}
	return s;
}



void
bigfloat_int::absolute_value ()
{
	bigfloat u = Upp;
	bigfloat l = Low;
	if (Upp.is_lt_zero()) {
		u = -Low;
		l = -Upp;
	}
	else {
		if (Low.is_lt_zero()) {
			l = 0;
			u = -Low;
			if (Upp.compare(u) == 1) {
				u = Upp;
			}
		}
        }
	Upp = u;
	Low = l;
}



void
bigfloat_int::negate ()
{
	bigfloat t = Upp;
	bigfloat::set_mode(MP_RND_UP);
	Upp = -Low;
	bigfloat::set_mode(MP_RND_DOWN);
	Low = -t;
}



void
bigfloat_int::invert ()
{
	bigfloat_int z;
	bigfloat a1, a2, a3, a4;

	bigfloat::set_mode(MP_RND_UP);
	a1 = inverse(Low);
	a2 = inverse(Upp);
	bigfloat::set_mode(MP_RND_DOWN);
	a3 = inverse(Low);
	a4 = inverse(Upp);
	z.normalize(a1, a2, a3, a4);
	Low = z.Low;
	Upp = z.Upp;
}



void
bigfloat_int::dble (double & d)
{
	bigfloat c(Low+Upp);

	c.divide_by_2();

	if (c.doublify(d) == 1) {
		lidia_error_handler("bigfloat_int", "dble:Cannot assign to double");
	}
}



void
bigfloat_int::Approx ()
{
	bigfloat c(Low+Upp);

	c.divide_by_2();
	Low = c;
	Upp = c;
}



void
bigfloat_int::approx (bigfloat & c)
{
	c = Low+Upp;
	c.divide_by_2();
}



// Length of interval
void
bigfloat_int::int_length (bigfloat& a)
{
	a = abs(Low-Upp);
}



//
// arithmetic procedures
//

void
add (bigfloat_int & sum, const bigfloat_int & x, const bigfloat_int & y)
{
	bigfloat::set_mode(MP_RND_DOWN);
	add(sum.Low, x.Low, y.Low);
	bigfloat::set_mode(MP_RND_UP);
	add(sum.Upp, x.Upp, y.Upp);
}



void
add (bigfloat_int & sum, const bigfloat_int & x, long i)
{
	bigfloat::set_mode(MP_RND_DOWN);
	add(sum.Low, x.Low, i);
	bigfloat::set_mode(MP_RND_UP);
	add(sum.Upp, x.Upp, i);
}



void
add (bigfloat_int & sum, long i, const bigfloat_int & y)
{
	bigfloat::set_mode(MP_RND_DOWN);
	add(sum.Low, y.Low, i);
	bigfloat::set_mode(MP_RND_UP);
	add(sum.Upp, y.Upp, i);
}



void
subtract (bigfloat_int & dif, const bigfloat_int & x, const bigfloat_int & y)
{
	bigfloat::set_mode(MP_RND_DOWN);
	subtract(dif.Low, x.Low, y.Upp);
	bigfloat::set_mode(MP_RND_UP);
	subtract(dif.Upp, x.Upp, y.Low);
}



void
subtract (bigfloat_int & dif, const bigfloat_int & x, long i)
{
	bigfloat::set_mode(MP_RND_DOWN);
	subtract(dif.Low, x.Low, i);
	bigfloat::set_mode(MP_RND_UP);
	subtract(dif.Upp, x.Upp, i);
}



void
subtract (bigfloat_int & dif, long i, const bigfloat_int & y)
{
	bigfloat::set_mode(MP_RND_DOWN);
	subtract(dif.Low, i, y.Upp);
	bigfloat::set_mode(MP_RND_UP);
	subtract(dif.Upp, i, y.Low);
}



void
multiply (bigfloat_int & prod, const bigfloat_int & x, const bigfloat_int & y)
{
	bigfloat::set_mode(MP_RND_DOWN);
	multiply(prod.Low, x.Low, y.Low);
	bigfloat::set_mode(MP_RND_UP);
	multiply(prod.Upp, x.Upp, y.Upp);
}



void
multiply (bigfloat_int & prod, const bigfloat_int & x, long i)
{
	bigfloat::set_mode(MP_RND_DOWN);
	multiply(prod.Low, x.Low, i);
	bigfloat::set_mode(MP_RND_UP);
	multiply(prod.Upp, x.Upp, i);
}



void
multiply (bigfloat_int & prod, long i, const bigfloat_int & x)
{
	bigfloat::set_mode(MP_RND_DOWN);
	multiply(prod.Low, x.Low, i);
	bigfloat::set_mode(MP_RND_UP);
	multiply(prod.Upp, x.Upp, i);
}



void
divide (bigfloat_int & quot, const bigfloat_int & x, const bigfloat_int & y)
{
	bigfloat::set_mode(MP_RND_DOWN);
	divide(quot.Low, x.Low, y.Upp);
	bigfloat::set_mode(MP_RND_UP);
	divide(quot.Upp, x.Upp, y.Low);
}



void
divide (bigfloat_int & quot, const bigfloat_int & x, long i)
{
	bigfloat::set_mode(MP_RND_DOWN);
	divide(quot.Low, x.Low, i);
	bigfloat::set_mode(MP_RND_UP);
	divide(quot.Upp, x.Upp, i);
}



void
divide (bigfloat_int & quot, long i, const bigfloat_int & x)
{
	bigfloat::set_mode(MP_RND_DOWN);
	divide(quot.Low, i, x.Upp);
	bigfloat::set_mode(MP_RND_UP);
	divide(quot.Upp, i, x.Low);

}



void
truncate (bigint & y, const bigfloat_int & x)
{
	long s = x.sign();

	floor(y, abs(x));
	multiply(y, y, s);
}



void
round (bigint & y, const bigfloat_int & x)
{
	bigint t;

	round(x.lower()).bigintify(t);
	round(x.upper()).bigintify(y);
	if (t != y)
		lidia_error_handler("bigfloat_int", "round()::Not Enough Precision In Round");
}



void
sqrt (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat::set_mode(MP_RND_DOWN);
	sqrt(y.Low, x.Low);
	bigfloat::set_mode(MP_RND_UP);
	sqrt(y.Upp, x.Upp);
}



void
power (bigfloat_int & z, const bigfloat_int & x, const bigfloat_int & y)
{
	bigfloat l1, l2, l3, l4, l5, l6;
	bigfloat_int s;

	bigfloat::set_mode(MP_RND_DOWN);
	power(l1, x.Low, y.Low);
	power(l2, x.Low, y.Upp);
	power(l3, x.Upp, y.Low);
	power(l4, x.Upp, y.Upp);
	z.normalize(l1, l2, l3, l4);
	l5 = z.Low; l6 = z.Upp;
	bigfloat::set_mode(MP_RND_UP);
	power(l1, x.Low, y.Low);
	power(l2, x.Low, y.Upp);
	power(l3, x.Upp, y.Low);
	power(l4, x.Upp, y.Upp);
	s.normalize(l1, l2, l3, l4);
	z.normalize(s.Upp, s.Low, l5, l6);
}



void
power (bigfloat_int & x, const bigfloat_int & a, long n)
{
	bigfloat::set_mode(MP_RND_DOWN);
	power(x.Low, a.Low, n);
	bigfloat::set_mode(MP_RND_UP);
	power(x.Upp, a.Upp, n);
}



//
// constants
//

void
constant_E (bigfloat_int & x)
{
	bigfloat r, l;

	bigfloat::set_mode(MP_RND_DOWN);
	constant_E(l);
	bigfloat::set_mode(MP_RND_UP);
	constant_E(r);
	x.assign(l, r);
}



void
constant_Pi (bigfloat_int & x)
{
	bigfloat r, l;

	bigfloat::set_mode(MP_RND_DOWN);
	constant_Pi(l);
	bigfloat::set_mode(MP_RND_UP);
	constant_Pi(r);
	x.assign(l, r);
}



void
constant_Euler (bigfloat_int & x)
{
	bigfloat r, l;

	bigfloat::set_mode(MP_RND_DOWN);
	constant_Euler(l);
	bigfloat::set_mode(MP_RND_UP);
	constant_Euler(r);
	x.assign(l, r);
}



void
constant_Catalan (bigfloat_int & x)
{
	bigfloat r, l;

	bigfloat::set_mode(MP_RND_DOWN);
	constant_Catalan(l);
	bigfloat::set_mode(MP_RND_UP);
	constant_Catalan(r);
	x.assign(l, r);
}



//
// I/O
//

std::ostream &
operator << (std::ostream & out, const bigfloat_int & a)
{
	out << "(" << a.Low << ", " << a.Upp << ")";
	return out;
}



std::istream &
operator >> (std::istream & in, bigfloat_int & a)
{
	char ch;
	in >> ch;
	if (ch != '(')
		lidia_error_handler("bigfloat_int", " >>:: Expected (");
	in >> a.Low;
	in >> ch;
	if (ch != ',')
		lidia_error_handler("bigfloat_int", " >>:: Expected , ");
	in >> a.Upp;
	in >> ch;
	if (ch != ')')
		lidia_error_handler("bigfloat_int", " >>:: Expected)");
	return in;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

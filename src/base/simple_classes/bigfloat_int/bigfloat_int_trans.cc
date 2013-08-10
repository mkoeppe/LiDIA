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



#include	"LiDIA/bigfloat_int.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void
exp (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat::set_mode(MP_RND_DOWN);
	exp(y.Low, x.Low);
	bigfloat::set_mode(MP_RND_UP);
	exp(y.Upp, x.Upp);
}



void
log (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat::set_mode(MP_RND_DOWN);
	log(y.Low, x.Low);
	bigfloat::set_mode(MP_RND_UP);
	log(y.Upp, x.Upp);
}



void
sin (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat u1, l1, u2, l2;

	bigfloat::set_mode(MP_RND_DOWN);
	sin(l1, x.Low);
	sin(u1, x.Upp);
	bigfloat::set_mode(MP_RND_UP);
	sin(l2, x.Low);
	sin(u2, x.Upp);
	y.normalize(l1, l2, u1, u2);
	l1 = x.Low; u1 = x.Upp;
	u2 = 2*Pi();
	while (l1 > u2) {
		l1 = l1-u2;
		u1 = u1-u2;
	}
	while (l1.is_lt_zero()) {
		l1 = l1+u2;
		u1 = u1+u2;
	}
	l2 = Pi()/2;
	if (l1< l2 && u1 > l2) {
		y.Upp = 1.0;
	}
	l2 = u2-l2;
	if (l1< l2 && u1 > l2) {
		y.Low = -1.0;
	}
}



void
cos (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat u1, l1, u2, l2;

	bigfloat::set_mode(MP_RND_DOWN);
	cos(l1, x.Low);
	cos(u1, x.Upp);
	bigfloat::set_mode(MP_RND_UP);
	cos(l2, x.Low);
	cos(u2, x.Upp);
	y.normalize(l1, l2, u1, u2);
	l1 = x.Low; u1 = x.Upp;
	u2 = 2*Pi();
	while (l1 > u2) {
		l1 = l1-u2;
		u1 = u1-u2;
	}
	while (l1.is_lt_zero()) {
		l1 = l1+u2;
		u1 = u1+u2;
	}
	l2 = Pi();
	if (l1< l2 && u1 > l2) {
		y.Low = -1.0;
	}
	if (l1< u2 && u1 > u2) {
		y.Upp = 1.0;
	}
}



void
tan (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat u1, l1, u2, l2;

	bigfloat::set_mode(MP_RND_DOWN);
	tan(l1, x.Low);
	tan(u1, x.Upp);
	bigfloat::set_mode(MP_RND_UP);
	tan(l2, x.Low);
	tan(u2, x.Upp);
	y.normalize(l1, l2, u1, u2);
	l1 = x.Low; u1 = x.Upp;
	u2 = Pi()/2;
	while (l1 > u2) {
		l1 = l1-u2;
		u1 = u1-u2;
	}
	while (l1.is_lt_zero()) {
		l1 = l1+u2;
		u1 = u1+u2;
	}
	if (l1< u2 && u1 > u2) {
		lidia_error_handler("bigfloat_int", "tan::Interval contains a multiple of pi/2");
	}
}



void
cot (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat u1, l1, u2, l2;

	bigfloat::set_mode(MP_RND_DOWN);
	cot(l1, x.Low);
	cot(u1, x.Upp);
	bigfloat::set_mode(MP_RND_UP);
	cot(l2, x.Low);
	cot(u2, x.Upp);
	y.normalize(l1, l2, u1, u2);
	l1 = x.Low; u1 = x.Upp;
	u2 = Pi();
	while (l1 > u2) {
		l1 = l1-u2;
		u1 = u1-u2;
	}
	while (l1.is_lt_zero()) {
		l1 = l1+u2;
		u1 = u1+u2;
	}
	if (l1< u2 && u1 > u2) {
		lidia_error_handler("bigfloat_int", "tan::Interval contains a multiple of pi");
	}
}



void
asin (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat u1, l1, u2, l2;

	bigfloat::set_mode(MP_RND_DOWN);
	asin(l1, x.Low);
	asin(u1, x.Upp);
	bigfloat::set_mode(MP_RND_UP);
	asin(l2, x.Low);
	asin(u2, x.Upp);
	y.normalize(l1, l2, u1, u2);
}



void
acos (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat u1, l1, u2, l2;

	bigfloat::set_mode(MP_RND_DOWN);
	acos(l1, x.Low);
	acos(u1, x.Upp);
	bigfloat::set_mode(MP_RND_UP);
	acos(l2, x.Low);
	acos(u2, x.Upp);
	y.normalize(l1, l2, u1, u2);
}



void
atan (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat u1, l1, u2, l2;

	bigfloat::set_mode(MP_RND_DOWN);
	atan(l1, x.Low);
	atan(u1, x.Upp);
	bigfloat::set_mode(MP_RND_UP);
	atan(l2, x.Low);
	atan(u2, x.Upp);
	y.normalize(l1, l2, u1, u2);
}



void
atan2 (bigfloat_int & z, const bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat u1, l1, u2, l2, u3, l3;
	bigfloat_int s;

	bigfloat::set_mode(MP_RND_DOWN);
	atan2(l1, y.Low, x.Low);
	atan2(u1, y.Low, x.Upp);
	atan2(l2, y.Upp, x.Low);
	atan2(u2, y.Upp, x.Upp);
	z.normalize(l1, l2, u1, u2);
	u3 = z.Upp;
	l3 = z.Low;
	bigfloat::set_mode(MP_RND_UP);
	atan2(l1, y.Low, x.Low);
	atan2(u1, y.Low, x.Upp);
	atan2(l2, y.Upp, x.Low);
	atan2(u2, y.Upp, x.Upp);
	s.normalize(l1, u1, l2, u2);
	z.normalize(l3, u3, s.Upp, s.Low);
	if (x.Upp.is_le_zero() && y.Upp.is_approx_zero()
	    && y.Low.is_approx_zero()) {
		z.Low = z.Upp;
	}
}



void
acot (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int p1, p2;

	atan(p1, x);
	constant_Pi(p2);
	divide(p2, p2, 2);
	subtract(y, p2, p1);
}



void
sinh (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int tmp(x);

	exp(y, tmp);
	divide(tmp, 1, y);
	subtract(y, y, tmp);
	divide(y, y, 2);
}



void
cosh (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int tmp(x);

	exp(y, tmp);
	divide(tmp, 1, y);
	add(y, y, tmp);
	divide(y, y, 2);
}



void
tanh (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int tmp1, tmp2;

	exp(y, x);
	divide(tmp1, 1, y);
	add(tmp2, y, tmp1);
	subtract(y, y, tmp1);
	divide(y, y, tmp2);
}



void
coth (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int tmp1, tmp2;

	exp(y, x);
	divide(tmp1, 1, y);
	subtract(tmp2, y, tmp1);
	add(y, y, tmp1);
	divide(y, y, tmp2);
}



void
asinh (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int p1;

	square(p1, x);
	inc(p1);
	sqrt(p1, p1);
	add(y, x, p1);
	log(y, y);
}



void
acosh (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int p1;

	square(p1, x);
	dec(p1);
	if (p1.is_lt_zero())
		lidia_error_handler("bigfloat_int", "acosh::cannot handle arguments <= 1");
	sqrt(p1, p1);
	add(y, x, p1);
	log(y, y);
}



void
atanh (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int p1(x);

	p1.negate();
	inc(p1);
	divide(y, 2, p1);
	dec(y);
	if (p1.is_lt_zero())
		lidia_error_handler("bigfloat_int", "atanh::cannot handle arguments <= -1 and >= 1");
	log(y, y);
	divide(y, y, 2);
}



void
acoth (bigfloat_int & y, const bigfloat_int & x)
{
	bigfloat_int p1(x);

	dec(p1);
	divide(y, 2, p1);
	inc(p1);
	if (p1.is_lt_zero())
		lidia_error_handler("bigfloat_int", "acoth::cannot handle arguments >= -1 and <= 1");
	log(y, y);
	divide(y, y, 2);
}



void
besselj (bigfloat_int & y, int n, const bigfloat_int & x)
{
	bigfloat u1, l1, u2, l2;

	bigfloat::set_mode(MP_RND_DOWN);
	besselj(l1, n, x.Low);
	besselj(u1, n, x.Upp);
	bigfloat::set_mode(MP_RND_UP);
	besselj(l2, n, x.Low);
	besselj(u2, n, x.Upp);
	y.normalize(l1, l2, u1, u2);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

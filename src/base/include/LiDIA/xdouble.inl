// -*- C++ -*-
//
// Addition
//
/*
Copyright (C) 1997 Keith Martin Briggs

Except where otherwise indicated,
this program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

// 97 Aug 04 KMB added ldexp
// 97 Jul 11 Moved this stuff out of quad.h, created inline.h so it can
//	be #included even if we're not inlining, by just "#define inline"
// 97 Jul 12 Added all combos of doubledouble/double/int +-*/.  Only a couple actually
//	written; most just call the others by swapping arguments.  Generates
//	equivalent code with a good inlining compiler (tested with g++ 2.7.2).
//	- e.g., all subtraction is now done by addition of unary minus.
//	- removed if's from doubledouble*int. Zero-branch code is 2.5 faster (tested)
//	- generally cleaned up and re-organized the order of the functions,
//	  now they're grouped nicely by function.
//	- tested Newton's division.  Works but is terribly slow, because
//	  it needs to do several doubledouble + and * operations, and doubledouble /
//	  without it is about the same speed at doubledouble * anyway, so no win.
//	- moved recip here as an inline.
//	- checked and tested doubledouble/double (BUG?), seems fine.

//
// Addition
//
/*      (C) W. Kahan 1989
*       NOTICE:
*       Copyrighted programs may not be translated, used, nor
*       reproduced without the author's permission.  Normally that
*       permission is granted freely for academic and scientific
*       purposes subject to the following three requirements:
*       1. This NOTICE and the copyright notices must remain attached
*          to the programs and their translations.
*       2. Users of such a program should regard themselves as voluntary
*          participants in the author's researches and so are obliged
*          to report their experience with the program back to the author.
*       3. Neither the author nor the institution that employs him will
*          be held responsible for the consequences of using a program
*          for which neither has received payment.
*       Would-be commercial users of these programs must correspond
*       with the author to arrange terms of payment and warranty.
*/



//
// Addition
//
XDBL_INLINE xdouble operator + (const xdouble& x, const xdouble& y)
{
	x86_FIX
		double H, h, T, t, S, s, e, f;
	S = x.hi+y.hi; T = x.lo+y.lo; e = S-x.hi; f = T-x.lo; s = S-e; t = T-f;
	s = (y.hi-e)+(x.hi-s); t = (y.lo-f)+(x.lo-t);
	e = s+T; H = S+e; h = e+(S-H); e = t+h;
	xdouble z; z.hi = H+e; z.lo = e+double(H-z.hi);
	END_x86_FIX
		return z;
}

XDBL_INLINE xdouble operator + (const double& x, const xdouble& y)
{
	x86_FIX
		double H, h, S, s, e;
	S = x+y.hi; e = S-x; s = S-e;
	s = (y.hi-e)+(x-s); H = S+(s+y.lo); h = (s+y.lo)+(S-H);
	xdouble z; z.hi = H+h; z.lo = h+double(H-z.hi);
	END_x86_FIX
		return z;
}

XDBL_INLINE xdouble operator + (const xdouble& x, const double& d)
{ return (static_cast<double>(d) + x); }
XDBL_INLINE xdouble operator + (const unsigned long& ul, const xdouble& x)
{ return (static_cast<double>(ul) + x); }
XDBL_INLINE xdouble operator + (const xdouble& x, const unsigned long& ul)
{ return (static_cast<double>(ul) + x); }
XDBL_INLINE xdouble operator + (const long& l, const xdouble& x)
{ return (static_cast<double>(l) + x); }
XDBL_INLINE xdouble operator + (const xdouble& x, const long& l)
{ return (static_cast<double>(l) + x); }
XDBL_INLINE xdouble operator + (const unsigned int& ui, const xdouble& x)
{ return (static_cast<double>(ui) + x); }
XDBL_INLINE xdouble operator + (const xdouble& x, const unsigned int& ui)
{ return (static_cast<double>(ui) + x); }
XDBL_INLINE xdouble operator + (const int& i, const xdouble& x)
{ return (static_cast<double>(i) + x); }
XDBL_INLINE xdouble operator + (const xdouble& x, const int& i)
{ return (static_cast<double>(i) + x); }
XDBL_INLINE xdouble operator + (const unsigned short& us, const xdouble& x)
{ return (static_cast<double>(us) + x); }
XDBL_INLINE xdouble operator + (const xdouble& x, const unsigned short& us)
{ return (static_cast<double>(us) + x); }
XDBL_INLINE xdouble operator + (const short& s, const xdouble& x)
{ return (static_cast<double>(s) + x); }
XDBL_INLINE xdouble operator + (const xdouble& x, const short& s)
{ return (static_cast<double>(s) + x); }
XDBL_INLINE xdouble operator + (const unsigned char& uc, const xdouble& x)
{ return (static_cast<double>(uc) + x); }
XDBL_INLINE xdouble operator + (const xdouble& x, const unsigned char& uc)
{ return (static_cast<double>(uc) + x); }
XDBL_INLINE xdouble operator + (const char& c, const xdouble& x)
{ return (static_cast<double>(c) + x); }
XDBL_INLINE xdouble operator + (const xdouble& x, const char& c)
{ return (static_cast<double>(c) + x); }

//
// Assignment & Addition
//
XDBL_INLINE xdouble& xdouble::operator += (const xdouble& y)
{
	x86_FIX
		double H, h, T, t, S, s, e, f;
	S = hi+y.hi; T = lo+y.lo; e = S-hi; f = T-lo; s = S-e; t = T-f;
	s = (y.hi-e)+(hi-s); t = (y.lo-f)+(lo-t); f = s+T; H = S+f; h = f+(S-H);
	hi = H+(t+h); lo = (t+h)+double(H-hi);
	END_x86_FIX
		return (*this);
}

XDBL_INLINE xdouble& xdouble::operator += (const double& d)
{
	x86_FIX
		double H, h, S, s, e, f;
	S = hi+d; e = S-hi; s = S-e;
	s = (d-e)+(hi-s); f = s+lo; H = S+f; h = f+(S-H);
	hi = H+h; lo = h+double(H-hi);
	END_x86_FIX
		return (*this);
}
XDBL_INLINE xdouble& xdouble::operator += (const unsigned long& ul)
{ *this += static_cast<double>(ul); return (*this); }
XDBL_INLINE xdouble& xdouble::operator += (const long& l)
{ *this += static_cast<double>(l); return (*this); }
XDBL_INLINE xdouble& xdouble::operator += (const unsigned int& ui)
{ *this += static_cast<double>(ui); return (*this); }
XDBL_INLINE xdouble& xdouble::operator += (const int& i)
{ *this += static_cast<double>(i); return (*this); }
XDBL_INLINE xdouble& xdouble::operator += (const unsigned short& us)
{ *this += static_cast<double>(us); return (*this); }
XDBL_INLINE xdouble& xdouble::operator += (const short& s)
{ *this += static_cast<double>(s); return (*this); }
XDBL_INLINE xdouble& xdouble::operator += (const unsigned char& uc)
{ *this += static_cast<double>(uc); return (*this); }
XDBL_INLINE xdouble& xdouble::operator += (const char& c)
{ *this += static_cast<double>(c); return (*this); }

//
// Unary minus
//
XDBL_INLINE xdouble operator - (const xdouble & x)
{
	return (xdouble(-x.hi, -x.lo));
}

//
// Absolute value
//
XDBL_INLINE xdouble fabs(const xdouble & x)
{
	if (x.h() >= 0.0)
		return (x);
	else
		return (-x);
}

//
// Subtraction
//
XDBL_INLINE xdouble operator - (const xdouble& x, const xdouble& y)
{
	return (x + xdouble(-y.hi, -y.lo));
}
XDBL_INLINE xdouble operator - (const double& x, const xdouble& y)
{
	return (x + xdouble(-y.hi, -y.lo));
}

//
// Subtraction
//
XDBL_INLINE xdouble operator - (const xdouble& x, const double& d)
{ return (-d + x); }
XDBL_INLINE xdouble operator - (const unsigned long& ul, const xdouble& x)
{ return (static_cast<double>(ul) - x); }
XDBL_INLINE xdouble operator - (const xdouble& x, const unsigned long& ul)
{ return (-(static_cast<double>(ul)) + x); }
XDBL_INLINE xdouble operator - (const long& l, const xdouble& x)
{ return (static_cast<double>(l) - x); }
XDBL_INLINE xdouble operator - (const xdouble& x, const long& l)
{ return (-(static_cast<double>(l)) + x); }
XDBL_INLINE xdouble operator - (const unsigned int& ui, const xdouble& x)
{ return (static_cast<double>(ui) - x); }
XDBL_INLINE xdouble operator - (const xdouble& x, const unsigned int& ui)
{ return (-(static_cast<double>(ui)) + x); }
XDBL_INLINE xdouble operator - (const int& i, const xdouble& x)
{ return (static_cast<double>(i) - x); }
XDBL_INLINE xdouble operator - (const xdouble& x, const int& i)
{ return (-(static_cast<double>(i)) + x); }
XDBL_INLINE xdouble operator - (const unsigned short& us, const xdouble& x)
{ return (static_cast<double>(us) - x); }
XDBL_INLINE xdouble operator - (const xdouble& x, const unsigned short& us)
{ return (-(static_cast<double>(us)) + x); }
XDBL_INLINE xdouble operator - (const short& s, const xdouble& x)
{ return (static_cast<double>(s) - x); }
XDBL_INLINE xdouble operator - (const xdouble& x, const short& s)
{ return (-(static_cast<double>(s)) + x); }
XDBL_INLINE xdouble operator - (const unsigned char& uc, const xdouble& x)
{ return (static_cast<double>(uc) - x); }
XDBL_INLINE xdouble operator - (const xdouble& x, const unsigned char& uc)
{ return (-(static_cast<double>(uc)) + x); }
XDBL_INLINE xdouble operator - (const char& c, const xdouble& x)
{ return (static_cast<double>(c) - x); }
XDBL_INLINE xdouble operator - (const xdouble& x, const char& c)
{ return (-(static_cast<double>(c)) + x); }

//
// Assignment & Subtraction
//
XDBL_INLINE xdouble& xdouble::operator -= (const xdouble& x)
{ *this = *this - x; return (*this); }
XDBL_INLINE xdouble& xdouble::operator -= (const double& d)
{ *this = *this - d; return (*this); }
XDBL_INLINE xdouble& xdouble::operator -= (const unsigned long& ul)
{ *this = *this - ul; return (*this); }
XDBL_INLINE xdouble& xdouble::operator -= (const long& l)
{ *this = *this - l; return (*this); }
XDBL_INLINE xdouble& xdouble::operator -= (const unsigned int& ui)
{ *this = *this - ui; return (*this); }
XDBL_INLINE xdouble& xdouble::operator -= (const int& i)
{ *this = *this - i; return (*this); }
XDBL_INLINE xdouble& xdouble::operator -= (const unsigned short& us)
{ *this = *this - us; return (*this); }
XDBL_INLINE xdouble& xdouble::operator -= (const short& s)
{ *this = *this - s; return (*this); }
XDBL_INLINE xdouble& xdouble::operator -= (const unsigned char& uc)
{ *this = *this - uc; return (*this); }
XDBL_INLINE xdouble& xdouble::operator -= (const char& c)
{ *this = *this - c; return (*this); }

//
// Multiplication
//
XDBL_INLINE xdouble operator *(const xdouble& x, const xdouble& y)
{
	x86_FIX
		double hx, tx, hy, ty, C, c;
	C = xdouble::Split*x.hi; hx = C-x.hi; c = xdouble::Split*y.hi;
	hx = C-hx; tx = x.hi-hx; hy = c-y.hi;
	C = x.hi*y.hi; hy = c-hy; ty = y.hi-hy;
	c = ((((hx*hy-C)+hx*ty)+tx*hy)+tx*ty)+(x.hi*y.lo+x.lo*y.hi);
	xdouble z; z.hi = C+c; hx = C-z.hi; z.lo = c+hx;
	END_x86_FIX
		return z;
}

// double*xdouble
XDBL_INLINE xdouble operator *(const double& x, const xdouble& y)
{
	x86_FIX
		double hx, tx, hy, ty, C, c;
	C = xdouble::Split*x; hx = C-x; c = xdouble::Split*y.hi; hx = C-hx;
	tx = x-hx; hy = c-y.hi; C = x*y.hi; hy = c-hy; ty = y.hi - hy;
	c = ((((hx*hy-C)+hx*ty)+tx*hy)+tx*ty)+x*y.lo;
	xdouble z; z.hi = C+c; z.lo = c+double(C-z.hi);
	END_x86_FIX
		return z;
}

//
// Multiplication
//
XDBL_INLINE xdouble operator * (const xdouble& x, const double& d)
{ return (d * x); }
XDBL_INLINE xdouble operator * (const unsigned long& ul, const xdouble& x)
{ return (static_cast<double>(ul) * x); }
XDBL_INLINE xdouble operator * (const xdouble& x, const unsigned long& ul)
{ return (static_cast<double>(ul) * x); }
XDBL_INLINE xdouble operator * (const long& l, const xdouble& x)
{ return (static_cast<double>(l) * x); }
XDBL_INLINE xdouble operator * (const xdouble& x, const long& l)
{ return (static_cast<double>(l) * x); }
XDBL_INLINE xdouble operator * (const unsigned int& ui, const xdouble& x)
{ return (static_cast<double>(ui) * x); }
XDBL_INLINE xdouble operator * (const xdouble& x, const unsigned int& ui)
{ return (static_cast<double>(ui) * x); }
XDBL_INLINE xdouble operator * (const int& i, const xdouble& x)
{ return (static_cast<double>(i) * x); }
XDBL_INLINE xdouble operator * (const xdouble& x, const int& i)
{ return (static_cast<double>(i) * x); }
XDBL_INLINE xdouble operator * (const unsigned short& us, const xdouble& x)
{ return (static_cast<double>(us) * x); }
XDBL_INLINE xdouble operator * (const xdouble& x, const unsigned short& us)
{ return (static_cast<double>(us) * x); }
XDBL_INLINE xdouble operator * (const short& s, const xdouble& x)
{ return (static_cast<double>(s) * x); }
XDBL_INLINE xdouble operator * (const xdouble& x, const short& s)
{ return (static_cast<double>(s) * x); }
XDBL_INLINE xdouble operator * (const unsigned char& uc, const xdouble& x)
{ return (static_cast<double>(uc) * x); }
XDBL_INLINE xdouble operator * (const xdouble& x, const unsigned char& uc)
{ return (static_cast<double>(uc) * x); }
XDBL_INLINE xdouble operator * (const char& c, const xdouble& x)
{ return (static_cast<double>(c) * x); }
XDBL_INLINE xdouble operator * (const xdouble& x, const char& c)
{ return (static_cast<double>(c) * x); }

//
// Assignment & Multiplication
//
XDBL_INLINE xdouble& xdouble::operator *= (const xdouble& x)
{
	x86_FIX
		double hx, tx, hy, ty, C, c;
	C = Split*hi; hx = C-hi; c = Split*x.hi;
	hx = C-hx; tx = hi-hx; hy = c-x.hi;
	C = hi*x.hi; hy = c-hy; ty = x.hi-hy;
	c = ((((hx*hy-C)+hx*ty)+tx*hy)+tx*ty)+(hi*x.lo+lo*x.hi);
	hx = C+c; hi = hx; lo = c+double(C-hx);
	END_x86_FIX
		return (*this);
}

XDBL_INLINE xdouble& xdouble::operator *= (const double& d)
{
	x86_FIX
		double hx, tx, hy, ty, C, c;
	C = Split*hi; hx = C-hi; c = Split*d;
	hx = C-hx; tx = hi-hx; hy = c-d;
	C = hi*d; hy = c-hy; ty = d-hy;
	c = ((((hx*hy-C)+hx*ty)+tx*hy)+tx*ty)+(lo*d);
	hx = C+c; hi = hx; lo = c+double(C-hx);
	END_x86_FIX
		return (*this);
}

XDBL_INLINE xdouble& xdouble::operator *= (const unsigned long& ul)
{ *this *= static_cast<double>(ul); return (*this); }
XDBL_INLINE xdouble& xdouble::operator *= (const long& l)
{ *this *= static_cast<double>(l); return (*this); }
XDBL_INLINE xdouble& xdouble::operator *= (const unsigned int& ui)
{ *this *= static_cast<double>(ui); return (*this); }
XDBL_INLINE xdouble& xdouble::operator *= (const int& i)
{ *this *= static_cast<double>(i); return (*this); }
XDBL_INLINE xdouble& xdouble::operator *= (const unsigned short& us)
{ *this *= static_cast<double>(us); return (*this); }
XDBL_INLINE xdouble& xdouble::operator *= (const short& s)
{ *this *= static_cast<double>(s); return (*this); }
XDBL_INLINE xdouble& xdouble::operator *= (const unsigned char& uc)
{ *this *= static_cast<double>(uc); return (*this); }
XDBL_INLINE xdouble& xdouble::operator *= (const char& c)
{ *this *= static_cast<double>(c); return (*this); }


//
// Division
//
// Reciprocal
XDBL_INLINE xdouble recip(const xdouble& y)
{
	x86_FIX
		double  hc, tc, hy, ty, C, c, U, u;
	C = 1.0/y.h();
	c = xdouble::Split*C;
	hc = c-C;
	u = xdouble::Split*y.h();
	hc = c-hc; tc = C-hc; hy = u-y.h(); U = C*y.h(); hy = u-hy; ty = y.h()-hy;
	u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
	c = ((((1.0-U)-u))-C*y.l())/y.h();
	xdouble z; z.hi = C+c; z.lo = double(C-z.hi)+c;
	END_x86_FIX
		return (z);
}

XDBL_INLINE xdouble operator / (const xdouble& x, const xdouble& y)
{
	x86_FIX
		double hc, tc, hy, ty, C, c, U, u;
	C = x.hi/y.hi; c = xdouble::Split*C; hc = c-C; u = xdouble::Split*y.hi; hc = c-hc;
	tc = C-hc; hy = u-y.hi; U = C * y.hi; hy = u-hy; ty = y.hi-hy;
	u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
	c = ((((x.hi-U)-u)+x.lo)-C*y.lo)/y.hi;
	xdouble z; z.hi = C+c; z.lo = double(C-z.hi)+c;
	END_x86_FIX
		return (z);
}

//
// double/xdouble:
//
XDBL_INLINE xdouble operator / (const double& x, const xdouble& y)
{
	x86_FIX
		double  hc, tc, hy, ty, C, c, U, u;
	C = x/y.hi; c = xdouble::Split*C; hc = c-C; u = xdouble::Split*y.hi; hc = c-hc;
	tc = C-hc; hy = u-y.hi; U = C*y.hi; hy = u-hy; ty = y.hi-hy;
	u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
	c = ((((x-U)-u))-C*y.lo)/y.hi;
	xdouble z; z.hi = C+c; z.lo = double(C-z.hi)+c;
	END_x86_FIX
		return (z);
}

//
// xdouble/double
//
XDBL_INLINE xdouble operator / (const xdouble& x, const double& y)
{
	x86_FIX
		double hc, tc, hy, ty, C, c, U, u;
	C = x.hi/y; c = xdouble::Split*C; hc = c-C; u = xdouble::Split*y; hc = c-hc;
	tc = C-hc; hy = u-y; U = C*y; hy = u-hy; ty = y-hy;
	u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
	c = (((x.hi-U)-u)+x.lo)/y;
	xdouble z; z.hi = C+c; z.lo = double(C-z.hi)+c;
	END_x86_FIX
		return (z);
}

//
// Division
//
XDBL_INLINE xdouble operator / (const xdouble& x, const unsigned long& ul)
{ return (x / static_cast<double>(ul)); }
XDBL_INLINE xdouble operator / (const unsigned long& ul, const xdouble& x)
{ return (static_cast<double>(ul) / x); }
XDBL_INLINE xdouble operator / (const xdouble& x, const long& l)
{ return (x / static_cast<double>(l)); }
XDBL_INLINE xdouble operator / (const long& l, const xdouble& x)
{ return (static_cast<double>(l) / x); }
XDBL_INLINE xdouble operator / (const unsigned int& ui, const xdouble& x)
{ return (static_cast<double>(ui) / x); }
XDBL_INLINE xdouble operator / (const xdouble& x, const unsigned int& ui)
{ return (x / static_cast<double>(ui)); }
XDBL_INLINE xdouble operator / (const int& i, const xdouble& x)
{ return (static_cast<double>(i) / x); }
XDBL_INLINE xdouble operator / (const xdouble& x, const int& i)
{ return (x / static_cast<double>(i)); }
XDBL_INLINE xdouble operator / (const unsigned short& us, const xdouble& x)
{ return (static_cast<double>(us) / x); }
XDBL_INLINE xdouble operator / (const xdouble& x, const unsigned short& us)
{ return (x / static_cast<double>(us)); }
XDBL_INLINE xdouble operator / (const short& s, const xdouble& x)
{ return (static_cast<double>(s) / x); }
XDBL_INLINE xdouble operator / (const xdouble& x, const short& s)
{ return (x / static_cast<double>(s)); }
XDBL_INLINE xdouble operator / (const unsigned char& uc, const xdouble& x)
{ return (static_cast<double>(uc) / x); }
XDBL_INLINE xdouble operator / (const xdouble& x, const unsigned char& uc)
{ return (x / static_cast<double>(uc)); }
XDBL_INLINE xdouble operator / (const char& c, const xdouble& x)
{ return (static_cast<double>(c) / x); }
XDBL_INLINE xdouble operator / (const xdouble& x, const char& c)
{ return (x / static_cast<double>(c)); }


//
// Assignment & Division
//
XDBL_INLINE xdouble& xdouble::operator /= (const xdouble& x)
{
	x86_FIX
		double hc, tc, hy, ty, C, c, U, u;
	C = hi/x.hi; c = Split*C; hc = c-C; u = Split*x.hi; hc = c-hc;
	tc = C-hc; hy = u-x.hi; U = C * x.hi; hy = u-hy; ty = x.hi-hy;
	u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
	c = ((((hi-U)-u)+lo)-C*x.lo)/x.hi;
	u = C+c; hi = u; lo = double(C-u)+c;
	END_x86_FIX
		return (*this);
}

XDBL_INLINE xdouble& xdouble::operator /= (const double& d)
{
	x86_FIX
		double hc, tc, hy, ty, C, c, U, u;
	C = hi/d; c = Split*C; hc = c-C; u = Split*d; hc = c-hc;
	tc = C-hc; hy = u-d; U = C * d; hy = u-hy; ty = d-hy;
	u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
	c = (((hi-U)-u)+lo)/d;
	u = C+c; hi = u; lo = double(C-u)+c;
	END_x86_FIX
		return (*this);
}

XDBL_INLINE xdouble& xdouble::operator /= (const unsigned long& ul)
{ *this /= static_cast<double>(ul); return (*this); }
XDBL_INLINE xdouble& xdouble::operator /= (const long& l)
{ *this /= static_cast<double>(l); return (*this); }
XDBL_INLINE xdouble& xdouble::operator /= (const unsigned int& ui)
{ *this /= static_cast<double>(ui); return (*this); }
XDBL_INLINE xdouble& xdouble::operator /= (const int& i)
{ *this /= static_cast<double>(i); return (*this); }
XDBL_INLINE xdouble& xdouble::operator /= (const unsigned short& us)
{ *this /= static_cast<double>(us); return (*this); }
XDBL_INLINE xdouble& xdouble::operator /= (const short& s)
{ *this /= static_cast<double>(s); return (*this); }
XDBL_INLINE xdouble& xdouble::operator /= (const unsigned char& uc)
{ *this /= static_cast<double>(uc); return (*this); }
XDBL_INLINE xdouble& xdouble::operator /= (const char& c)
{ *this /= static_cast<double>(c); return (*this); }

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
//	$Id: xdouble.cc,v 2.5 2004/06/15 10:19:28 lidiaadm Exp $
//
//	Author	: Keith Briggs (KB)
//	Changes	: See CVS log
//
//==============================================================================================

// C++ functions for xdouble (double+double) precision.
// Use with xdouble.h

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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/xdouble.h"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef LIDIA_NAMESPACE
using std::floor;	// required for the macro lidia_rint
#endif



//
// Useful constants
//
const xdouble xdouble::Log2 = "0.6931471805599453094172321214581765680755";
const xdouble xdouble::Log10 = "2.302585092994045684017991454684364207601";
const xdouble xdouble::Pi = "3.1415926535897932384626433832795028841972";
const xdouble xdouble::TwoPi = "6.2831853071795864769252867665590057683943";
const xdouble xdouble::Pion2 = "1.5707963267948966192313216916397514420985";
const xdouble xdouble::Pion4 = "0.7853981633974483096156608458198757210493";
const xdouble xdouble::_Pi = "0.3183098861837906715377675267450287240689";


const double xdouble::Split = 134217729.0L;


//
// xdbl_long_abs() replaces abs(), because with CC,
// abs() is only defined for int-parameter, but
// it is used with long-parameter in this code.
//

static inline long xdbl_long_abs (long l)
{
	return ((l > 0) ? l : (-l));
}



//
// member constants
//
#if 0
// TPf: not used...
void
base_and_prec ()
{
//
// Linnainmaa ACM TOMS 7, 272 Thm 3
//
	int p;
	x86_FIX
		std::cerr << "Base and precision determination by Linnainmaa's method:" << std::endl;
	{
		double U, R, u, r, beta;
		U = 4.0/3.0;
		U -= 1; U *= 3; U -= 1; U = std::fabs(U);
		R = U/2+1; R -= 1;
		if (R != 0.0) U = R;
		u = 2; u /= 3; u -= 0.5; u *= 3; u -= 0.5;
		u = std::fabs(u); r = u/2+0.5; r -= 0.5;
		if (r != 0.0) u = r;
		beta = U/u; p = static_cast<int>(-std::log(u)/std::log(beta)+0.5);
		std::cout << "Type double: base is " << beta << ", precision is " << p << std::endl;
	}
	{
		xdouble U, R, u, r, beta;
		U = 4; U /= 3;
		U -= 1; U *= 3; U -= 1; U = fabs(U);
		R = U/2+1; R -= 1;
		if (R.h() != 0.0) U = R;
		u = 2; u /= 3; u -= 0.5; u *= 3; u -= 0.5;
		u = fabs(u); r = u/2+0.5; r -= 0.5;
		if (r.h() != 0.0) u = r;
		beta = U/u;
		p = int_xdbl(-log(u)/log(beta)+0.5);
		std::cout << "Type xdouble:   base is " << beta << ", precision is " << p << std::endl;
	}
	END_x86_FIX
		return;
}
#endif // if 0



//
// number of decimal digits to which x and y agree
//
int
digits (const xdouble& x, const xdouble& y)
{
	xdouble diff = fabs(x-y);
	if (diff.h() == 0.0)
		return (32);
	int d = -int_xdbl(floor((0.4*log((diff/fabs(x))))).h());
	return (d < 32?d:32);
}



//
// define copysign if needed
//
#ifndef HAVE_COPYSIGN
inline double
copysign (double x, double y)
{
	if (y >= 0.0)
		return std::fabs(x);
	else
		return -std::fabs(x);
}
#endif



inline xdouble
copysign (const xdouble& x, const double y)
{
	if (y >= 0.0)
		return fabs(x);
	else
		return -fabs(x);
}



//
// String conversions
// xdouble -> string  (modified code of operatior <<)
// string -> xdouble  (previous atoq (compiler problems))
//
void
string_to_xdouble (const char *s, xdouble& xd)
{
	xdouble result = 0.0; int n, sign, ex = 0;
	/* eat whitespace */ while (*s == ' ' || *s == '\t' || *s == '\n') s++;
	switch (*s) { // get sign of mantissa
	case '-': { sign = -1; s++; break; }
	case '+': s++; // no break
	default: sign = 1;
	}
	// get digits before decimal point
	while (n = (*s++)-'0', n >= 0 && n < 10) result = 10.0*result+n;
	s--;
	if (*s == '.')   /* get digits after decimal point */ {
		s++;
		while (n = (*s++)-'0', n >= 0 && n < 10) { result = 10.0*result+n; --ex; }
		s--;
	}
	if (*s == 'e' || *s == 'E') /* get exponent */ ex += atoi(++s);
	// exponent adjustment
	// std::cerr<<"atodd: result="<<std::endl; result.dump(""); std::cerr<<"atodd: ex="<<ex<<std::endl;
	while (ex-- > 0) result *= 10;
	while (++ex < 0) result /= 10;
	xd = ((sign >= 0)?result:-result);
}



void
xdouble_to_string (const xdouble& xd, char *s)
{
	char *p = s;
	if (xd.h() == 0.0) {
		sprintf(s, "0.0");
		return;
	}
	long Digits = 34;
	xdouble ten = 10.0, y = fabs(xd);
	double q = std::log10(y.h());
	long m, n = static_cast<long>(std::floor(q));
	if (n < 0)
		n++;
	xdouble l = powint(ten, n);
	y = y / l;
	if (sign(xd) < 0)
		*p++ = '-';
	long d = Digits > 34 ? 34 : Digits;
	d = d < 3 ? 3 : d;
	for (long i = 1; i <= d; i++) {
		if (i == 2)
			*p++ = '.';
		m = static_cast<long>(std::floor(y.h()));
		sprintf(p, "%ld", m);
		p++;
		y = (y - xdouble(m)) * ten;
		if (y.h() < 0.0)
			break; // x must be an longeger

	}
	if (n != 0) {
		*p++ = 'e';
		sprintf(p, "%ld", n);
	}
	else
		*p = 0;
}



//
// Constructor from string
//
xdouble::xdouble (char *s)
{
	xdouble y;
	string_to_xdouble(s, y);
	hi = y.hi;
	lo = y.lo;
}



//
// xdouble = string
//
xdouble &
xdouble::operator = (char *s)
{
	xdouble y;
	string_to_xdouble(s, y);
	hi = y.hi;
	lo = y.lo;
	return *this;
}



//
// Debugging use
//
void
xdouble::dump (char *s = "")
{
	std::cerr << s << "xdouble(" << hi << ", " << lo << ")\n";
}



//
// Input
//
std::istream &
operator >> (std::istream& s, xdouble& x)
{
	xdouble result = 0.0; int n, sign = 1, ex = 0; char ch;
	s >> std::ws; // eat whitespace
	s >> ch;
	if (ch == '-') { sign = -1; s >> ch; } else if (ch == '+') { s >> ch; }
	// get digits before decimal point
	n = ch-'0';
	while (n >= 0 && n < 10) {
		result = 10*result+n;
		s.get(ch); // cannot use s >> ch; as that will eat whitespace
		if (ch == ' ' || ch == '\t' || ch == '\n') goto fixup;
		n = ch-'0';
	}
	if (ch == '.') {  // get digits after decimal point
		s >> ch; n = ch-'0';
		while (n >= 0 && n < 10) {
			result = 10*result+n;
			s.get(ch);
			n = ch-'0'; --ex;
			if (ch == ' ' || ch == '\t' || ch == '\n') goto fixup;
		}
	}
	n = 0;
	if (ch == 'e' || ch == 'E') { s >> n; ex += n; } // get exponent
	else s.putback(ch);
 fixup:
	if (sign < 0) result = -result;
	// exponent adjustment
	while (ex-- > 0) result *= 10;
	while (++ex < 0) result /= 10;
	x = result;
	return (s);
}



//
// Output
//
std::ostream &
operator << (std::ostream& s, const xdouble& x) {
	if (x.h() == 0.0) { s << "0.0 "; return s; }
	if (x.h() != x.h()) { s << "NaN "; return s; }
	int Digits = s.precision();
	xdouble ten = 10.0, y = fabs(x); double q = std::log10(y.h());
	int m, n = static_cast<int>(std::floor(q));
	if (n < 0) n++;
	xdouble l = powint(ten, n);
	y = y/l;
	if (sign(x) < 0) s << "-"; //else s << " ";
	int d = Digits > 34 ? 34 : Digits;
	d = d < 3 ? 3 : d;
	for (int i = 1; i <= d; i++) {
		if (2 == i) s << ".";
		m = int_xdbl(floor(y));
		if (m< 0 || m > 9)
		{ lidia_error_handler("xdouble", "operator << (s, x) "
				      "Internal error");
		}
		s << m;
		y = (y-xdouble(m))*ten;
		if (y.h() <= 0.0) break; // x must be an integer
		//if (y.h()==0.0) break; // ???????????
	}
	if (n != 0) s << "e" << n; else s << "";
	s << "";
	return s;
}



//
// ieee - Functions (Part I)
//
//
// rint (round to nearest int)
//
xdouble
rint (const xdouble& x)
{
	return floor(x + xdouble(0.5));
}



//
// Another floor. V. Shoup 97 Sep 23...
//
xdouble
floor_s (const xdouble& x) {
	double fhi = std::floor(x.h());
	if (fhi != x.h())
		return xdouble(fhi);
	else
		return xdouble(fhi, std::floor(x.l()));
}



// Floor.  See Graham, Knuth and Patashnik `Concrete Mathematics' page 70.
// Greatest integer <= x
// maybe floor_s is better?
xdouble
floor (const xdouble& x)
{
	double fh = std::floor(x.h()), fl = std::floor(x.l());
	double t1 = x.h()-fh;
	double t2 = x.l()-fl;
	double t3 = t1+t2;
	t3 = x.h()-fh+x.l()-fl;
	int t = static_cast<int>(std::floor(t3));
	switch (t) {
	case 0: return xdouble(fh)+xdouble(fl);
	case 1: return xdouble(fh)+xdouble(fl+1);
		// case 2 can only arise when fp({hi})<=1, fp({lo})<=1, but fp({hi}+{lo})=2
	case 2: return xdouble(fh)+xdouble(fl+1);
	default: lidia_error_handler("xdouble::floor(xdouble)",
				     "Invalid case.");
	// MM: MS Visual C++ 5.0 needs a return
	return xdouble(fh);
	}
	// never get here

}



xdouble
ceil (const xdouble& x)
{
	return -floor(-x);
}



xdouble
trunc (const xdouble& x)
{
	if (x.h() >= 0.0)
		return floor(x);
	else
		return -floor(-x);
}



xdouble
fmod (const xdouble& x, const int n)
{
	return x-n*floor(x/n);
}



// Same as the ANSI modf for doubles.
// BEWARE: contains ugly, magical code
//
xdouble
modf (const xdouble &D, xdouble *id) {
	int sign = 1;
	xdouble qFrac, d = D;
	double lowFrac, highFrac, lowInt, highInt;
	if (d < xdouble(0.0)) { sign = -1; d = -d; }
	if (d < xdouble(1.0)) { // it's all fraction, no integer
		*id = 0; return sign*d; }
	if (d + 1 == d) { // It's all integer, no fraction
		*id = sign*d; return 0.0; }
	highFrac = std::modf(d.h(), &highInt);
	lowFrac = std::modf(d.l(), &lowInt);
	/* special case: if d.l is opposite in sign to d.h, then adding them
	** together changes the integer part.
	*/
	if (highInt == d.h() && d.l() < 0) {
#if 1
		if (lowFrac != 0) {
			*id = xdouble(highInt) + xdouble(lowInt) - 1;
			qFrac = 1 + xdouble(lowFrac);
		} else {
			qFrac = lowFrac;
			*id = xdouble(highInt) + xdouble(lowInt);
		}
#else
		// loop method, slow but should work
		*id = 0;
		top = QUAD_BITS;
		while (true) {
			// Invariant: *id <= d
			if (!(*id <= d)) printf("Oops, !(*id <= d): %f %f\n\n", *id, d);
			if (d < *id + 1) break;
			for (i = top-1; i >= 0 && (*id + std::ldexp(1., i) > d); i--);
			*id += std::ldexp(1., i);
			top = i;
		}
		*id *= sign;
		qFrac = d - sign*(*id);
		printf("loop Int = %.16g %.16g, frac = %.16g %.16g\n",
		       id->h, id->l, qFrac.h, qFrac.l);
#endif
	}
	else {
		*id = xdouble(highInt) + xdouble(lowInt);
		qFrac = xdouble(highFrac) + xdouble(lowFrac);
	}
	if (sign == -1) { *id = -*id; return -qFrac; }
	else return qFrac;
}



//
// Signum
//
int
sign (const xdouble& x)
{
	if (x.h() > 0.0)
		return 1;
	else
		if (x.h() < 0.0)
			return -1;
		else
			return 0;
}



//
// Comparison
//

bool
operator > (const xdouble & x, const xdouble & y)
{
	return (x.h() > y.h()) || (x.h() == y.h() && x.l() > y.l());
}



bool
operator >= (const xdouble & x, const xdouble & y)
{
	return (x.h() >= y.h()) || (x.h() == y.h() && x.l() >= y.l());
}



bool
operator < (const xdouble & x, const xdouble & y)
{
	return (x.h() < y.h()) || (x.h() == y.h() && x.l() < y.l());
}



bool
operator <= (const xdouble & x, const xdouble & y)
{
	return (x.h() <= y.h()) || (x.h() == y.h() && x.l() <= y.l());
}



bool
operator == (const xdouble & x, const xdouble & y)
{
	return x.h() == y.h() && x.l() == y.l();
}



bool
operator != (const xdouble & x, const xdouble & y)
{
	return x.h() != y.h() || x.l() != y.l();
}


#ifdef XDBL_NO_INLINE
#include	"LiDIA/xdouble.inl"
#endif






//
// ieee - Functions (Part II)
//

//
// Square  (faster than x*x)
//
xdouble
sqr (const xdouble & x)
{
	x86_FIX
		double hx, tx, C, c;
	C = xdouble::Split*x.h(); hx = C-x.h(); hx = C-hx; tx = x.h()-hx;
	C = x.h()*x.h();
	c = ((((hx*hx-C)+2.0*hx*tx))+tx*tx)+2.0*x.h()*x.l();
	hx = C+c;
	xdouble r = xdouble(hx, c+(C-hx));
	END_x86_FIX
		return r;
}



//
// cube
//
xdouble
cub (const xdouble & x)
{
	xdouble z = x * sqr(x);
	return z;
}



//
// ldexp
//
xdouble
ldexp(const xdouble& x, const int exp)
{
	// x*2^exp
	return xdouble(std::ldexp(x.h(), exp), std::ldexp(x.l(), exp));
}



xdouble
exp (const xdouble& x)
{
/* Uses some ideas by Alan Miller
   Method:
   x    x.log2(e)    nint[x.log2(e)] + frac[x.log2(e)]
   e = 2 = 2

   iy    fy
   = 2   . 2
   Then
   fy    y.loge(2)
   2 = e

   Now y.loge(2) will be less than 0.3466 in absolute value.
   This is halved and a Pade' approximant is used to approximate e^x over
   the region (-0.1733, +0.1733).   This approximation is then squared.
*/
	if (x.h() < -744.4400719213812) return 0.0; // exp(x) < 1e-300
	int iy;
	xdouble y, temp, ysq, sum1, sum2;
	y = x/xdouble::Log2;
	temp = iy = int_xdbl(rint(y));
	y = (y-temp)*xdouble::Log2;
	y = ldexp(y, -1);
	ysq = sqr(y);
	sum1 = y*((((ysq+3960.)*ysq+2162160.)*ysq+302702400.)*ysq+8821612800.);
	sum2 = (((90.*ysq+110880.)*ysq+30270240.)*ysq+2075673600.)*ysq+17643225600.;
/*
  sum2 + sum1         2.sum1
  Now approximation = ----------- = 1 + ----------- = 1 + 2.temp
  sum2 - sum1       sum2 - sum1

  Then (1 + 2.temp)^2 = 4.temp.(1 + temp) + 1
*/
	temp = sum1/(sum2-sum1);
	y = temp*(temp+1);
	y = ldexp(y, 2);
	return ldexp(y+1, iy);
}



//
// See Higham "Accuracy and Stability of Numerical Algorithms", p 511
//
xdouble
hypot (const xdouble& a, const xdouble& b)
{
	xdouble p, q, r, s, aa, ab, four = 4.0;
	aa = fabs(a); ab = fabs(b);
	if (aa > ab) { p = aa; q = ab; } else { q = aa; p = ab; } // now p >= q
	if (0.0 == p.h()) return q;
	while (true) {
		r = sqr(q/p);
		if (four == (aa = r+four)) return p;
		s = r/aa; p += 2*s*p; q *= s;
	} // Never get here
}



//
// square root
//
xdouble
sqrt (const xdouble& y)
{
	x86_FIX
		double c, p, q, hx, tx, u, uu, cc, hi = y.h();
	if (hi < 0.0) {
		lidia_error_handler("xdouble", "sqrt(y) :: y < 0");
	}
	if (0.0 == hi) return y;
	c = std::sqrt(hi);
	p = xdouble::Split*c; hx = c-p; hx += p; tx = c-hx;
	p = hx*hx;
	q = 2.0*hx*tx;
	u = p+q;
	uu = (p-u)+q+tx*tx;
	cc = (((y.h()-u)-uu)+y.l())/(c+c);
	u = c+cc;
	xdouble r = xdouble(u, cc+(c-u));
	END_x86_FIX
		return r;
}



//
// Natural logarithm
//
xdouble
log (const xdouble& x)
{
	// Newton method
	if (x.h() < 0.0) {
		lidia_error_handler("xdouble", "log(x) :: x < 0 !");
	}
	if (x.h() == 0.0) {
		lidia_error_handler("xdouble", "log(x) :: x == 0 !");
	}
	xdouble s = std::log(x.h()), e = exp(s); // s = double approximation to result
	return s+(x-e)/e; // Newton correction, good enough
	//xdouble d=(x-e)/e; return s+d-0.5*sqr(d);  // or 2nd order correction
}



//
// logarithm base 10
//
xdouble
log10 (const xdouble& t) {
	static const xdouble one_on_log10 = "0.4342944819032518276511289189166050822944";
	return one_on_log10*log(t);
}



//
// xdouble^xdouble
//
xdouble
pow (const xdouble& a, const xdouble& b)
{
	return exp(b*log(a));
}



//
// xdouble^int
//
xdouble
powint (const xdouble& u, const long c)
{
	if (c < 0) return recip(powint(u, -c));
	switch (c) {
	case 0: return u.h() == 0.0?xdouble(pow(0.0, 0.0)):xdouble(1.0); // let math.h handle 0^0.
	case 1: return u;
	case 2: return sqr(u);
	case 3: return sqr(u)*u;
	default: { // binary method
		int n = c, m; xdouble y = 1.0, z = u; if (n < 0) n = -n;
		while (true) {
			m = n; n /= 2;
			if (n+n != m) { // if m odd
				y *= z; if (0 == n) return y;
			}
			z = sqr(z);
		}
	}
	}
}



//
// Like Miller's modr
// a=n*b+rem, |rem|<=b/2, exact result.
xdouble
modr (const xdouble& a, const xdouble& b, int& n, xdouble& rem)
{
	if (0.0 == b.h()) {
		lidia_error_handler("xdouble", "modr(a, b, n, rem) :: "
				    "xdouble: divisor is zero in modr!");
	}
	xdouble temp;
	temp = a/b;
	n = static_cast<int>(lidia_rint(temp.h()));
	temp = n*xdouble(b.h());
	rem = xdouble(a.h());
	temp = rem-temp;
	rem = xdouble(a.l());
	temp = rem+temp;
	rem = n*xdouble(b.l());
	rem = temp-rem;
	return rem;
}



xdouble
sin (const xdouble& x)
{
	static const xdouble tab[9] = { // tab[b] : = sin(b*Pi/16)...
		0.0,
		"0.1950903220161282678482848684770222409277",
		"0.3826834323650897717284599840303988667613",
		"0.5555702330196022247428308139485328743749",
		"0.7071067811865475244008443621048490392850",
		"0.8314696123025452370787883776179057567386",
		"0.9238795325112867561281831893967882868225",
		"0.9807852804032304491261822361342390369739",
		1.0
	};
	static const xdouble sinsTab[7] = { // Chebyshev coefficients
		"0.9999999999999999999999999999993767021096",
		"-0.1666666666666666666666666602899977158461",
		"8333333333333333333322459353395394180616.0e-42",
		"-1984126984126984056685882073709830240680.0e-43",
		"2755731922396443936999523827282063607870.0e-45",
		"-2505210805220830174499424295197047025509.0e-47",
		"1605649194713006247696761143618673476113.0e-49"
	};
	if (std::fabs(x.h()) < 1.0e-7) return x*(1.0-sqr(x)/3);
	int a, b; xdouble sins, coss, k1, k3, t2, s, s2, sinb, cosb;
	// reduce x: -Pi < x <= Pi
	k1 = x/xdouble::TwoPi; k3 = k1-rint(k1);
	// determine integers a and b to minimize |s|, where  s=x-a*Pi/2-b*Pi/16
	t2 = 4*k3;
	a = int_xdbl(rint(t2));
	b = int_xdbl(rint((8*(t2-a)))); // must have |a| <= 2 and |b| <= 7 now
	s = xdouble::Pi*(k3+k3-(8*a+b)/16.0); // s is now reduced argument. -Pi/32 < s < Pi/32
	s2 = sqr(s);
	// Chebyshev series on -Pi/32..Pi/32, max abs error 2^-98=3.16e-30, whereas
	// power series has error 6e-28 at Pi/32 with terms up to x^13 included.
	sins = s*(sinsTab[0]+(sinsTab[1]+(sinsTab[2]+(sinsTab[3]+(sinsTab[4]+
								  (sinsTab[5]+sinsTab[6]*s2)*s2)*s2)*s2)*s2)*s2);
	coss = sqrt(1.0-sqr(sins)); // ok as -Pi/32 < s < Pi/32
	// here sinb=sin(b*Pi/16) etc.
	sinb = (b >= 0) ? tab[b] : -tab[-b]; cosb = tab[8-xdbl_long_abs(b)];
	if (0 == a) {
		return  sins*cosb+coss*sinb;
	} else if (1 == a) {
		return -sins*sinb+coss*cosb;
	} else if (-1 == a) {
		return  sins*sinb-coss*cosb;
	} else  { // |a |= 2
		return -sins*cosb-coss*sinb;
	}
	// i.e. return sins*(cosa*cosb-sina*sinb)+coss*(sina*cosb+cosa*sinb);
}



//
// sin and cos.   Faster than separate calls of sin and cos.
//
void
sincos (const xdouble& x, xdouble& sinx, xdouble& cosx)
{
	static const xdouble tab[9] = { // tab[b] : = sin(b*Pi/16)...
		0.0,
		"0.1950903220161282678482848684770222409277",
		"0.3826834323650897717284599840303988667613",
		"0.5555702330196022247428308139485328743749",
		"0.7071067811865475244008443621048490392850",
		"0.8314696123025452370787883776179057567386",
		"0.9238795325112867561281831893967882868225",
		"0.9807852804032304491261822361342390369739",
		1.0
	};
	static const xdouble sinsTab[7] = { // Chebyshev coefficients
		"0.9999999999999999999999999999993767021096",
		"-0.1666666666666666666666666602899977158461",
		"8333333333333333333322459353395394180616.0e-42",
		"-1984126984126984056685882073709830240680.0e-43",
		"2755731922396443936999523827282063607870.0e-45",
		"-2505210805220830174499424295197047025509.0e-47",
		"1605649194713006247696761143618673476113.0e-49"
	};
	if (std::fabs(x.h()) < 1.0e-11) { sinx = x; cosx = 1.0-0.5*sqr(x); return; }
	int a, b; xdouble sins, coss, k1, k3, t2, s, s2, sinb, cosb;
	k1 = x/xdouble::TwoPi; k3 = k1-rint(k1);
	t2 = 4*k3;
	a = int_xdbl(rint(t2));
	b = int_xdbl(rint((8*(t2-a))));
	s = xdouble::Pi*(k3+k3-(8*a+b)/16.0);
	s2 = sqr(s);
	sins = s*(sinsTab[0]+(sinsTab[1]+(sinsTab[2]+(sinsTab[3]+(sinsTab[4]+
								  (sinsTab[5]+sinsTab[6]*s2)*s2)*s2)*s2)*s2)*s2);
	coss = sqrt(1.0-sqr(sins)); // ok, sins is small
	sinb = (b >= 0) ? tab[b] : -tab[-b]; cosb = tab[8-xdbl_long_abs(b)];
	// sin(x)=
	// sin(s)(cos(1/2 a Pi) cos(1/16 b Pi) - sin(1/2 a Pi) sin(1/16 b Pi))
	// cos(s)(sin(1/2 a Pi) cos(1/16 b Pi) + cos(1/2 a Pi) sin(1/16 b Pi))
	// cos(x)=
	// cos(s)(cos(1/2 a Pi) cos(1/16 b Pi) - sin(1/2 a Pi) sin(1/16 b Pi))
	//-sin(s)(sin(1/2 a Pi) cos(1/16 b Pi) + cos(1/2 a Pi) sin(1/16 b Pi))
	if (0 == a) {
		sinx = sins*cosb+coss*sinb;
		cosx = coss*cosb-sins*sinb;
	} else if (1 == a) {
		sinx = -sins*sinb+coss*cosb;
		cosx = -coss*sinb-sins*cosb;
	} else if (-1 == a) {
		sinx = sins*sinb-coss*cosb;
		cosx = coss*sinb+sins*cosb;
	} else  { // |a |= 2
		sinx = -sins*cosb-coss*sinb;
		cosx = -coss*cosb+sins*sinb;
	}
	return;
}



//
// cos
//
xdouble
cos (const xdouble& x)
{
	return sin(xdouble::Pion2-x);
}



//
// hyperbolic
//
xdouble
sinh (const xdouble& x)
{
	if (std::fabs(x.h()) < 1.0e-5) { // avoid cancellation in e^x-e^(-x), use Taylor series...
		xdouble q = sqr(x); return x*(1+q/6*(1+q/20*(1+q/42))); }
	xdouble t = exp(x);
	return 0.5*(t-recip(t));
}



xdouble
cosh(const xdouble& x)
{
	xdouble t = exp(x);
	return 0.5*(t+recip(t));
}



xdouble
tanh (const xdouble& z)
{
	xdouble e;
	if (z.h() > 0.0) {
		e = exp(-2.0*z);
		return (1.0-e)/(1.0+e);
	}
	else {
		e = exp(2.0*z);
		return (e-1.0)/(1.0+e);
	}
}



xdouble
atan (const xdouble& x)
{
	double xh = x.h();
	if (0.0 == xh) return xdouble(0.0);
	// Asymptotic formula for large |x|...
	if (std::fabs(xh) > 1.0e6)
	{ xdouble r = recip(x), r2 = sqr(r);
	return xdouble::Pion2-r+r2*r*(1.0-(6*r2)/10)/3; }
	xdouble s, c, a = std::atan(xh); // a = double approx to result
	sincos(a, s, c);
	return a+c*(c*x-s); // Newton step
}



xdouble
atan2(const xdouble& qy, const xdouble& qx)
{
	// Based on GNU libc atan2.c
	static const double one = 1.0, zero = 0.0; double x, y;
	x = qx.h(); y = qy.h();
	double signx, signy;
	if (x != x) return qx; // x = NaN
	if (y != y) return qy;
	signy = copysign(one, y);
	signx = copysign(one, x);
	if (y == zero) return signx == one ? qy : copysign(xdouble::Pi, signy);
	if (x == zero) return copysign(xdouble::Pion2, signy);
	if (__isinf(x)) {
		if (__isinf(y)) return copysign(signx == one ? xdouble::Pion4 : 3.0*xdouble::Pion4, signy);
		else            return copysign(signx == one ? xdouble(0.0) : xdouble::Pi, signy);
	}
	if (__isinf(y)) return copysign(xdouble::Pion2, signy);
	xdouble aqy = fabs(qy);
	if (x < 0.0) // X is negative.
		return copysign(xdouble::Pi-atan(aqy/(-qx)), signy);
	return copysign(atan(aqy/qx), signy);
}



//
// arcsin
//
xdouble
asin (const xdouble& x)
{
	if (fabs(x) > xdouble(1.0))
	{ lidia_error_handler("xdouble", "asin(x) :: |x| > 1 !"); }
	return atan2(x, sqrt(1.0-sqr(x)));
}



// KMB 96 Oct 28 version 0.3 97 Dec 22 OK now.
// erfc written 97 Dec 22, NOT CHECKED!!!!!!!
// Based on series and continued fraction for incomplete gamma function.
xdouble
erf (const xdouble& x)
{
	if (0.0 == x.h()) return 0.0;
	int i; xdouble y, r;
	// const below is evalf(1/sqrt(Pi),40)
	static const xdouble oneonrootpi = "0.564189583547756286948079451560772585844";
	const double cut = 1.5; // switch to continued fraction here
	y = fabs(x);
	if (y.h() > 26.0) { // erf is +or- 1 to 300 dp.
		r = 1;
	} else if (y.h() < cut) { // use power series
		xdouble ap = 0.5, s, t, x2;
		s = t = 2.0; x2 = sqr(x);
		for (i = 0; i < 200; i++) {
			ap += 1;
			t *= x2/ap;
			s += t;
			if (std::fabs(t.h()) < 1e-35*std::fabs(s.h())) {
				r = x*oneonrootpi*s/exp(x2);
				break;
			}
		}
		if (i > 199)
			lidia_error_handler("xdouble", "erf(const xdouble)::"
					    "no convergence in erf power series (internal error)");
	} else { // |x| >= cut, use continued fraction, Lentz method
		double an, small_value = 1e-300;
		xdouble b, c, d, h, del, x2;
		x2 = sqr(x);
		b = x2+0.5;
		c = 1.0e300;
		d = recip(b);
		h = d;
		for (i = 1; i < 300; i++) {
			an = i*(0.5-i);
			b += 2.0;
			d = an*d+b;
			if (std::fabs(d.h()) < small_value) d = small_value;
			c = b+an/c;
			if (std::fabs(c.h()) < small_value) c = small_value;
			d = recip(d);
			del = d*c;
			h *= del;
			if (del.h() == 1.0 && del.l() < 1.0e-30) break;
			if (299 == i)
				lidia_error_handler("xdouble", "erf(const xdouble)::"
						    "no convergence in erf continued fraction (internal error)");
		}
		r = 1.0-oneonrootpi*sqrt(x2)/exp(x2)*h;
	}
	if (x.h() > 0.0) return r; else return -r;
}



xdouble
erfc (const xdouble& x)
{
	if (0.0 == x.h()) return 1.0;
	if (x.h() < 0.0) return 1.0-erf(x);
	int i; xdouble y, r;
	// const below is evalf(1/sqrt(Pi),40)
	static const xdouble oneonrootpi = "0.564189583547756286948079451560772585844";
	const double cut = 1.5; // switch to continued fraction here
	y = fabs(x);
	if (y.h() < cut) { // use power series
		xdouble ap = 0.5, s, t, x2;
		s = t = 2.0; x2 = sqr(x);
		for (i = 0; i < 200; i++) {
			ap += 1;
			t *= x2/ap;
			s += t;
			if (std::fabs(t.h()) < 1e-35*std::fabs(s.h())) {
				r = 1.0-x*oneonrootpi*s/exp(x2);
				break;
			}
		}
		if (i > 199)
			lidia_error_handler("xdouble", "erfc(const xdouble)::"
					    "no convergence in erfc power series (internal error)");
	} else { // |x| >= cut, use continued fraction
		double an, small_value = 1e-300;
		xdouble b, c, d, h, del, x2;
		x2 = sqr(x);
		b = x2+0.5;
		c = 1e300;
		h = d = recip(b);
		for (i = 1; i < 300; i++) {
			an = i*(0.5-i);
			b += 2.0;
			d = an*d+b;
			if (std::fabs(d.h()) < small_value) d = small_value;
			c = b+an/c;
			if (std::fabs(c.h()) < small_value) c = small_value;
			d = recip(d);
			del = d*c;
			h *= del;
			if (del.h() == 1.0 && del.l() < 1.0e-30) break;
			if (299 == i)
				lidia_error_handler("xdouble", "erfc(const xdouble)::"
						    "no convergence in erfc continued fraction (internal error)");
		}
		r = oneonrootpi*sqrt(x2)/exp(x2)*h;
	}
	return r;
}



xdouble
gamma (const xdouble &x)
{
	const int n = 43; // don't really need so many!
	static const xdouble c[n] = { // Taylor coefficients for 1/gamma(1+x)-x...
		"+0.5772156649015328606065120900824024310421593359",
		"-0.6558780715202538810770195151453904812797663805",
		"-0.0420026350340952355290039348754298187113945004",
		"+0.1665386113822914895017007951021052357177815022",
		"-0.0421977345555443367482083012891873913016526841",
		"-0.0096219715278769735621149216723481989753629422",
		"+0.0072189432466630995423950103404465727099048009",
		"-0.0011651675918590651121139710840183886668093337",
		"-0.0002152416741149509728157299630536478064782419",
		"+0.0001280502823881161861531986263281643233948920",
		"-0.0000201348547807882386556893914210218183822948",
		"-0.0000012504934821426706573453594738330922423226",
		"+0.0000011330272319816958823741296203307449433240",
		"-0.0000002056338416977607103450154130020572836512",
		"+0.0000000061160951044814158178624986828553428672",
		"+0.0000000050020076444692229300556650480599913030",
		"-0.0000000011812745704870201445881265654365055777",
		"+1.0434267116911005104915403323122501914007098231E-10",
		"+7.7822634399050712540499373113607772260680861813E-12",
		"-3.6968056186422057081878158780857662365709634513E-12",
		"+5.1003702874544759790154813228632318027268860697E-13",
		"-2.0583260535665067832224295448552374197460910808E-14",
		"-5.3481225394230179823700173187279399489897154781E-15",
		"+1.2267786282382607901588938466224224281654557504E-15",
		"-1.1812593016974587695137645868422978312115572918E-16",
		"+1.1866922547516003325797772429286740710884940796E-18",
		"+1.4123806553180317815558039475667090370863507503E-18",
		"-2.2987456844353702065924785806336992602845059314E-19",
		"+1.7144063219273374333839633702672570668126560625E-20",
		"+1.3373517304936931148647813951222680228750594717E-22",
		"-2.0542335517666727893250253513557337966820379352E-22",
		"+2.7360300486079998448315099043309820148653116958E-23",
		"-1.7323564459105166390574284515647797990697491087E-24",
		"-2.3606190244992872873434507354275310079264135521E-26",
		"+1.8649829417172944307184131618786668989458684290E-26",
		"+2.2180956242071972043997169136268603797317795006E-27",
		"+1.2977819749479936688244144863305941656194998646E-28",
		"+1.1806974749665284062227454155099715185596846378E-30",
		"-1.1245843492770880902936546742614395121194117955E-30",
		"+1.2770851751408662039902066777511246477487720656E-31",
		"-7.3914511696151408234612893301085528237105689924E-33",
		"+1.1347502575542157609541652594693063930086121959E-35",
		"+4.6391346410587220299448049079522284630579686797E-35"
	};
	xdouble ss = x, f = 1.0, sum = 0.0, one = 1.0;
	while (ss > one) { ss -= 1; f *= ss; }
	while (ss < one) { f /= ss; ss += 1; }
	if (ss == one) return f;
	ss -= 1.0;
	for (int i = n-1; i >= 0; i--) sum = c[i]+ss*sum;
	return f/(ss*sum+1);
}

//
// end of xdouble.cc
//



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

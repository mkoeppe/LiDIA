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
//	Author	: Andreas M"uller (AM), Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



sdigit
gcd (sdigit a, sdigit b)
{
	sdigit r;

	if (b == 0)
		return a;
	if (a == 0)
		return b;

	while (true) {
		r = a % b;
		if (!r)
			return (b);
		a = b;
		b = r;
	}
}



#if 0
unsigned int
lidia_gcd(unsigned int a, unsigned int b)
{
	unsigned int r;

	if (b == 0)
		return a;
	if (a == 0)
		return b;

	while (true) {
		r = a % b;
		if (!r) return (b);
		a = b; b = r;
	}
}
#endif


#if 0
int
int_xgcd(int &xv, int &yv, int a, int b)
{

	char sa, sb, f;
	int s, t, u, v, xu, yu, hxv, hyv;

	if ((sa = (a < 0)))
		a = (-a);
	if (!b) {
		yv = 0;
		xv = (sa) ? -1 : 1;
		return (a);
	}
	if ((sb = (b < 0)))
		b = (-b);
	if (!a) {
		xv = 0;
		yv = (sb) ? -1 : 1;
		return (b);
	}
	if ((f = (b > a))) {
		t = a;
		a = b;
		b = t;
	}
	xu = hyv = 1;
	yu = hxv = 0;
	u = a;
	v = b;
	s = 0;
	for (;;) {
		if ((u & 1) || (v & 1))
			break;
		s++;
		u >>= 1;
		v >>= 1;
	}
	for (;;) {
		if (u & 1)
			break;
		u >>= 1;
		if ((xu & 1) || (yu & 1)) {
			xu += b;
			yu -= a;
		}
		xu >>= 1;
		yu >>= 1;
	}
	for (;;) {
		if (v & 1)
			break;
		v >>= 1;
		if ((hxv & 1) || (hyv & 1)) {
			hxv += b;
			hyv -= a;
		}
		hxv >>= 1;
		hyv >>= 1;
	}
	for (;;) {
		if (u == v)
			break;
		if (u > v) {
			u -= v;
			xu -= hxv;
			yu -= hyv;
			for (;;) {
				if (u & 1)
					break;
				u >>= 1;
				if ((xu & 1) || (yu & 1)) {
					xu += b;
					yu -= a;
				}
				xu >>= 1;
				yu >>= 1;
			}
		}
		else {
			v -= u;
			hxv -= xu;
			hyv -= yu;
			for (;;) {
				if (v & 1)
					break;
				v >>= 1;
				if ((hxv & 1) || (hyv & 1)) {
					hxv += b;
					hyv -= a;
				}
				hxv >>= 1;
				hyv >>= 1;
			}
		}
	}

	if (f) {
		t = hxv;
		hxv = hyv;
		hyv = t;
	}
	if (sa)
		hxv = -hxv;
	if (sb)
		hyv = -hyv;
	xv = hxv;
	yv = hyv;
	return (u << s);
}
#endif



sdigit
xgcd(sdigit & xv, sdigit &yv, sdigit a, sdigit b)
{

	char sa, sb, f;
	long s, t, u, v, xu, yu, hxv, hyv;

	if ((sa = (a < 0)))
		a = (-a);
	if (!b) {
		yv = 0;
		xv = (sa) ? -1 : 1;
		return (a);
	}
	if ((sb = (b < 0)))
		b = (-b);
	if (!a) {
		xv = 0;
		yv = (sb) ? -1 : 1;
		return (b);
	}
	if ((f = (b > a))) {
		t = a;
		a = b;
		b = t;
	}
	xu = hyv = 1;
	yu = hxv = 0;
	u = a;
	v = b;
	s = 0;
	for (;;) {
		if ((u & 1) || (v & 1))
			break;
		s++;
		u >>= 1;
		v >>= 1;
	}
	for (;;) {
		if (u & 1)
			break;
		u >>= 1;
		if ((xu & 1) || (yu & 1)) {
			xu += b;
			yu -= a;
		}
		xu >>= 1;
		yu >>= 1;
	}
	for (;;) {
		if (v & 1)
			break;
		v >>= 1;
		if ((hxv & 1) || (hyv & 1)) {
			hxv += b;
			hyv -= a;
		}
		hxv >>= 1;
		hyv >>= 1;
	}
	for (;;) {
		if (u == v)
			break;
		if (u > v) {
			u -= v;
			xu -= hxv;
			yu -= hyv;
			for (;;) {
				if (u & 1)
					break;
				u >>= 1;
				if ((xu & 1) || (yu & 1)) {
					xu += b;
					yu -= a;
				}
				xu >>= 1;
				yu >>= 1;
			}
		}
		else {
			v -= u;
			hxv -= xu;
			hyv -= yu;
			for (;;) {
				if (v & 1)
					break;
				v >>= 1;
				if ((hxv & 1) || (hyv & 1)) {
					hxv += b;
					hyv -= a;
				}
				hxv >>= 1;
				hyv >>= 1;
			}
		}
	}

	if (f) {
		t = hxv;
		hxv = hyv;
		hyv = t;
	}
	if (sa)
		hxv = -hxv;
	if (sb)
		hyv = -hyv;
	xv = hxv;
	yv = hyv;
	return (u << s);
}



unsigned int
power(unsigned int a, unsigned int b)
{
	unsigned int i, n = a;

	for (i = 1; i < b; i++)
		n *= a;

	return n;
}



#define MAXint       2147483647
#define MAXIH        1073741823
#define SBASIS            65536
#define MAXDBL       1000000000000000.0

int h_mul_mod(int a, int b, int m)
{
	register unsigned int n, n1, n2;
	double nd, d;

	if (m < SBASIS) {
		n = a * b;
		return(n % m);
	}

	d = 1.0 * a * b;

	if (d < MAXint) {
		n = static_cast<int>(d);
		n = n % m;
	}
	else
		if (d < MAXDBL) {
			n = static_cast<int>(d / m);
			nd = 1.0 * n;
			nd = nd * m;
			d = d - nd;	// (1.0 * n * m)
			n = static_cast<int>(d);
		}
		else {
			n1 = h_mul_mod(a % SBASIS, b, m);
			n2 = h_mul_mod(SBASIS , h_mul_mod(a / SBASIS, b, m) , m);
			if ((n1 < MAXIH) && (n2 < MAXIH))
				n = (n1 + n2) % m;
			else {
				d = 1.0 * n1 + 1.0 * n2;
				if (d < MAXint) {
					n = static_cast<int>(d);
					n = n % m;
				}
				else {
					n = static_cast<int>(d / m);
					d = d - (1.0 * n * m);
					n = static_cast<int>(d);
				}
			}
		}

	return(n);
}



#if defined(sparc) || defined(i386) || defined(mips) || defined(hppa)
int mul_mod(int a, int b, int m)
{
	int q, r;
	if (m <= 1)
		lidia_error_handler("mul_mod", "modul <= 1");

	q = static_cast<int>(static_cast<double>(a) * static_cast<double>(b) /m);
	r = a*b-q*m;

	if
		(r < 0) r += m;
	else
		if (r > m) r -= m;

	return r;
}



#else



int mul_mod(int a, int b, int m)
{
	register unsigned int n, n1, n2;
	double nd, d;

	if (m <= 1)
		lidia_error_handler("mul_mod", "modul <= 1");

	a = a % m;
	b = b % m;

	if (a < 0)
		a = a + m;
	if (b < 0)
		b = b + m;

	if (m < SBASIS) {
		n = a * b;
		return(n % m);
	}

	d = 1.0 * a * b;

	if (d < MAXint) {
		n = static_cast<int>(d);
		n = n % m;
	}
	else
		if (d < MAXDBL) {
			n = static_cast<int>(d / m);
			nd = 1.0 * n;
			nd = nd * m;
			d = d - nd;	// (1.0 * n * m)
			n = static_cast<int>(d);
		}
		else {
			n1 = h_mul_mod(a % SBASIS, b, m);
			n2 = h_mul_mod(SBASIS , h_mul_mod(a / SBASIS, b, m) , m);
			if ((n1 < MAXIH) && (n2 < MAXIH))
				n = (n1 + n2) % m;
			else {
				d = 1.0 * n1 + 1.0 * n2;
				if (d < MAXint) {
					n = static_cast<int>(d);
					n = n % m;
				}
				else {
					n = static_cast<int>(d / m);
					d = d - (1.0 * n * m);
					n = static_cast<int>(d);
				}
			}
		}

	return(n);
}
#endif



int power_mod(int bas, int exp, int modul)
{
	register int erg, ergz;

	if (exp < 0)
		lidia_error_handler("power_mod", "exponent < 0");

	if (modul <= 1)
		lidia_error_handler("power_mod", "modul <= 1");

	if (exp == 0)
		return(1);

	ergz = bas % modul;
	if (ergz == 0)
		return(0);
	if (ergz < 0)
		ergz += modul;

	erg = 1;

	while (true) {
		if (exp & 1) {
			erg = h_mul_mod(erg, ergz, modul);
			if (exp == 1)
				return(erg);
		}
		exp >>= 1;
		ergz = h_mul_mod(ergz, ergz, modul);
	}
}



int legendre(int n, int p)
{
	if (n == 0)
		return 0;
	else
		return (power_mod(n, (p-1) >> 1, p) == 1 ? 1 : -1);
}



int invert(int a, int b)
{
	register int r, q,
		xs = 0,
		xss = 1,
		lb = b;

	if (b == 0)
		return(-1);

	while (true) {
		r = a % lb;
		if (r == 0)
			break;
		q = a / lb;
		a = lb;
		lb = r;
		r = xss - q * xs;
		xss = xs;
		xs = r;
	}

	if (abs(lb) != 1)
		return (-1);

	r = (lb > 0) ? xs : -xs;

	if (r < 0) r += b;

	return r;
}

//--------------------------------------------------------------------

int jacobi(int a, int b)
{
  int k = 1, v, h;
  static int table[] = {0, 1, 0, -1, 0, -1, 0, 1};
  
  if (b == 0)   // trivial case are checked first
    {
      if (a < 0)
	a = -a;
      if (a == 1) 
	return (1);
      else 
	return (0);
    }

  if (!(a & 1) && !(b & 1))  // both even
    return (0);

  if (b < 0)
    {
      b = -b;
      if (a < 0)
	k = -1;
    }

  v = 0;
  while (!(b & 1))
    {
      v++;
      b >>= 1;
    }

  if (v & 1)
    k *= table[a & 7];

  if (a < 0)
    {
      if (b & 2)
	k = -k;
      a = -a;
    }

  while (a != 0)
    {
      v = 0;
      while (!(a & 1))
	{
	  v++;
	  a >>= 1;
	}
      if (v & 1)
	k *= table[b & 7];

      if (a < b)
	{
	  h = b;
	  b = a;
	  a = h;
	  if (a & b & 2)
	    k = -k;
	}
      a = a - b;
    }
  if (b == 1)
    return k;
  else 
    return 0;
}



int ressol(int a, int p)
{
	register int r, n, c, z, k, t, s;

	if (p == 2)    return(a);

	// (1) p = 3 mod 4
	if ((p & 0x3) == 3)    return(power_mod(a, (p+1) >> 2 , p));

	// (2) initialisation:
	// compute k, s : p = 2^s (2k+1) + 1
	k = p - 1;
	s = 0;

	while (!(k & 1)) {
		++s;
		k >>= 1;
	}

	k = (k - 1) >> 1;

	// (3) initial values
	r = power_mod(a, k, p);

	n = mul_mod(r, r, p);
	n = mul_mod(n, a, p);
	r = mul_mod(r, a, p);

	if (n == 1)
		return(r);

	// (4) no quadratic residue
	z = 2;
	while (jacobi(z, p) == 1)
		++z;

	c = power_mod(z, (k << 1) + 1, p);

	// (5) iteration
	while (n > 1) {
		k = n;
		t = s;
		s = 0;

		while (k != 1) {
			k = mul_mod(k, k, p);
			s++;
		}

		t = t - s;
		if(t == 0) {
		    lidia_error_handler("ressol",
					"a is not a quadratic residue mod p");
		}
		
		c = power_mod(c, 1 << (t-1), p);
		r = mul_mod(r, c, p);
		c = mul_mod(c, c, p);
		n = mul_mod(n, c, p);
	}

	return(r);
}



static unsigned char log2_table[] = {
	0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
	6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
};



int
integer_log(unsigned long x)
{
	int l = 0;
#if defined(__alpha)
	if (x & 0xffffffff00000000) l += 32, x >>= 32;
#endif
	if (x & 0xffff0000) l += 16, x >>= 16;
	if (x & 0xff00) l += 8, x >>= 8;
	return (l + log2_table[x]);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

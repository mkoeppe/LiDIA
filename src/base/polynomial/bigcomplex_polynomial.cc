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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigcomplex.h"
#include	"LiDIA/polynomial.h"
#include	"LiDIA/random_generator.h"

#ifndef LIDIA_INCLUDE_CC
# include	"LiDIA/base/poly_intern.cc"
# include	"LiDIA/field_polynomial.cc"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#define DV_CP LDBL_UNIPOL+25
#define DM_CP "polynomial< bigcomplex >"
#define LP_ERROR poly_error_msg


template <>
void base_polynomial< bigcomplex >::remove_leading_zeros()
{
	debug_handler_c(DM_CP, "in member - function remove_leading_zeros ()",
			DV_CP, std::cout << "original degree is " << deg << std::endl << std::flush);
	bigcomplex c, *np;
	lidia_size_t d, i;

	d = deg;
	np = coeff + d;
	while (d >= 0 && np->real().is_approx_zero()
	       && np->imag().is_approx_zero())
		d--, np--;

	if (d < 0) {
		deg = d;
		delete[] coeff;
		coeff = NULL;
	}
	else if (d != deg) {
		debug_handler_c(DM_CP, "in member - function remove_leading_zeros()",
				DV_CP + 1, std::cout << "new degree is " << d << std::endl << std::flush);
		deg = d;
		np = new bigcomplex[d + 1];
		memory_handler(np, DM_CP, "remove_leading_zeros() :: "
			       "Error in memory allocation (np)");
		for (i = 0; i <= d; i++)
			np[i] = coeff[i];
		delete[] coeff;
		coeff = np;
	}
}



// instantiate class base_polynomial
template class base_polynomial< bigcomplex >;

// instantiate friend functions -- arithmetical functions
template void negate(base_polynomial< bigcomplex > & c,
                     const base_polynomial< bigcomplex > &a);

template void add(base_polynomial< bigcomplex > & c,
                  const base_polynomial< bigcomplex > & a,
                  const base_polynomial< bigcomplex > & b);
template void add(base_polynomial< bigcomplex > & c,
                  const base_polynomial< bigcomplex > & a,
                  const bigcomplex & b);
template void add(base_polynomial< bigcomplex > & c,
                  const bigcomplex & b,
                  const base_polynomial< bigcomplex > & a);

template void subtract(base_polynomial< bigcomplex > & c,
                       const base_polynomial< bigcomplex > & a,
                       const base_polynomial< bigcomplex > & b);
template void subtract(base_polynomial< bigcomplex > & c,
                       const base_polynomial< bigcomplex > & a,
                       const bigcomplex & b);
template void subtract(base_polynomial< bigcomplex > & c,
                       const bigcomplex & b,
                       const base_polynomial< bigcomplex > & a);

template void multiply(base_polynomial< bigcomplex > & c,
                       const base_polynomial< bigcomplex > & a,
                       const base_polynomial< bigcomplex > & b);
template void multiply(base_polynomial< bigcomplex > & c,
                       const base_polynomial< bigcomplex > & a,
                       const bigcomplex & b);
template void multiply(base_polynomial< bigcomplex > & c,
                       const bigcomplex & b,
                       const base_polynomial< bigcomplex > & a);

template void power(base_polynomial< bigcomplex > & c,
                    const base_polynomial< bigcomplex > & a, const bigint & b);

template void derivative(base_polynomial< bigcomplex > &c,
                         const base_polynomial< bigcomplex > &a);
template base_polynomial< bigcomplex >
derivative(const base_polynomial< bigcomplex > &a);

// instantiate friend functions -- operators
template bool operator == (const base_polynomial< bigcomplex > &a,
			   const base_polynomial< bigcomplex > &b);
template bool operator != (const base_polynomial< bigcomplex > &a,
			   const base_polynomial< bigcomplex > &b);
template base_polynomial< bigcomplex >
operator -(const base_polynomial< bigcomplex > &a);

template base_polynomial< bigcomplex >
operator +(const base_polynomial< bigcomplex > &a,
	   const base_polynomial< bigcomplex > &b);

template base_polynomial< bigcomplex >
operator +(const base_polynomial< bigcomplex > &a,
	   const bigcomplex & b);

template base_polynomial< bigcomplex >
operator +(const bigcomplex & b,
	   const base_polynomial< bigcomplex > &a);

template base_polynomial< bigcomplex >
operator -(const base_polynomial< bigcomplex > &a,
	   const base_polynomial< bigcomplex > &b);

template base_polynomial< bigcomplex >
operator -(const base_polynomial< bigcomplex > &a,
	   const bigcomplex &b);

template base_polynomial< bigcomplex >
operator -(const bigcomplex &a,
	   const base_polynomial< bigcomplex > &b);

template base_polynomial< bigcomplex >
operator *(const base_polynomial< bigcomplex > &a,
	   const base_polynomial< bigcomplex > &b);

template base_polynomial< bigcomplex >
operator *(const base_polynomial< bigcomplex > &a,
	   const bigcomplex &b);

template base_polynomial< bigcomplex >
operator *(const bigcomplex &b,
	   const base_polynomial< bigcomplex > &a);

template std::istream & operator >> (std::istream &,
				     base_polynomial< bigcomplex > &);
template std::ostream & operator << (std::ostream &,
				     const base_polynomial< bigcomplex > &);

// instantiate friend functions -- access functions
template bigcomplex lead_coeff(const base_polynomial< bigcomplex > &a);
template bigcomplex const_term(const base_polynomial< bigcomplex > &a);

// instantiate class field_polynomial
template class field_polynomial< bigcomplex >;

// instantiate friend functions -- arithmetical functions
template void div_rem(field_polynomial< bigcomplex > &q,
                      field_polynomial< bigcomplex > &r,
                      const base_polynomial< bigcomplex > &a,
                      const base_polynomial< bigcomplex > &b);
template void divide(field_polynomial< bigcomplex > & c,
                     const base_polynomial< bigcomplex > & a,
                     const bigcomplex & b);
template void divide(field_polynomial< bigcomplex > & q,
		     const base_polynomial< bigcomplex > & a,
		     const base_polynomial< bigcomplex > & b);
template void remainder(field_polynomial< bigcomplex > & r,
			const base_polynomial< bigcomplex > & a,
			const base_polynomial< bigcomplex > & b);
template void power_mod(field_polynomial< bigcomplex > & c,
                        const base_polynomial< bigcomplex > & a,
                        const bigint & b,
                        const base_polynomial< bigcomplex > & f);

template void integral(field_polynomial< bigcomplex > &c,
                       const base_polynomial< bigcomplex > &a);

template void gcd(field_polynomial< bigcomplex > &d,
		  const base_polynomial< bigcomplex > &aa,
		  const base_polynomial< bigcomplex > &bb);

template void xgcd(field_polynomial< bigcomplex > &d,
		   field_polynomial< bigcomplex > &x,
		   field_polynomial< bigcomplex > &y,
		   const base_polynomial< bigcomplex > &aa,
		   const base_polynomial< bigcomplex > &bb);

bigcomplex root(const base_polynomial< bigcomplex > & p,
                const bigcomplex & start)
{
	debug_handler_l(DM_CP, "in member - function "
			"root (const base_polynomial< bigcomplex > &, "
			"const bigcomplex &)", DV_CP + 2);
	if (p.degree() <= 0) {
		lidia_error_handler("polynomial< bigcomplex >", "root(...): You can't "
				    "compute zeros of a constant polynomial");
		return bigcomplex();
	}
	int c, i = p.degree(), j;
	bigcomplex u, u1, dx, x(start), y, px(p[i]), p1x(0.0), p2x(0.0);
	bigfloat m, m1;
	bigfloat tmp;
	bigfloat lambda;

	for (j = i-1; j >= 0; j--) {
		multiply(p2x, x, p2x);
		add(p2x, p2x, p1x);
		multiply(p1x, x, p1x);
		add(p1x, p1x, px);
		multiply(px, x, px);
		add(px, px, p[j]);
	}

	u.assign(px);
	m.assign(norm(u));

	for (;;) {
		c = 0;
		tmp.assign(hypot(u));

		if (tmp.is_approx_zero())
			divide(dx, u, p1x);
		else {
			multiply(tmp, tmp, hypot(p2x));
			invert(lambda, tmp);
			multiply(lambda, lambda, norm(p1x) << 1);
			if (exponent(lambda) + lambda.bit_length() > 0)
				divide(dx, u, p1x);
			else {
				multiply(dx, lambda, u);
				divide(dx, dx, p1x);
			}
		}
		i = dx.real().is_approx_zero();
		j = dx.imag().is_approx_zero();
		if (i && j)
			break;
		do {
			subtract(y, x, dx);
			u1.assign(p(y));
			m1.assign(norm(u1));
			if (m1.is_approx_zero() || m1 < m)
				break;
			else {
				c++;
				divide(dx, dx, bigfloat(4.0));
				if (c > 20) {
					do {
						random_generator rg;
						unsigned long	r1, r2;

						rg >> r1 >> r2;
						add(x, start, bigcomplex(1 + 4*(static_cast<double>(r1)/(static_cast<double>(INT_MAX) + 1)),
									 static_cast<double>(r2)/(static_cast<double>(INT_MAX) + 1)));
						i = p.degree();
						px.assign(p[i]);
						p1x.assign_zero();
						p2x.assign_zero();
						for (j = i-1; j >= 0; j--) {
							multiply(p2x, x, p2x);
							add(p2x, p2x, p1x);
							multiply(p1x, x, p1x);
							add(p1x, p1x, px);
							multiply(px, x, px);
							add(px, px, p[j]);
						}
						u.assign(px);
						multiply(tmp, hypot(u), hypot(p2x));
						invert(lambda, tmp);
						multiply(lambda, lambda, (norm(p1x) << 1));
						if (exponent(lambda) + lambda.bit_length() >= 0)
							lambda = 1;
					} while (exponent(lambda) + lambda.bit_length() < -64);
					multiply(x, lambda, u);
					divide(x, x, p1x);
					m.assign(norm(u));
					c = 0;
				}
			}
		} while (m1 >= m);
		x.assign(y);
		u.assign(u1);
		m.assign(m1);
		i = p.degree();
		px.assign(p[i]);
		p1x.assign_zero();
		p2x.assign_zero();
		for (j = i-1; j >= 0; j--) {
			multiply(p2x, x, p2x);
			add(p2x, p2x, p1x);
			multiply(p1x, x, p1x);
			add(p1x, p1x, px);
			multiply(px, x, px);
			add(px, px, p[j]);
		}
	}
	return x;
}



//
// A root finding algorithm over C by H.Cohen
//

void
cohen(const base_polynomial< bigcomplex > & P, bigcomplex * wurz, int f, int &cnt)
{
	debug_handler_l(DM_CP, "in member - function "
			"cohen (const base_polynomial< bigcomplex > &, "
			"bigcomplex *, int, int &)", DV_CP + 3);

	polynomial< bigcomplex > R(bigcomplex(1.0)), R1;
	static bigfloat _two(2.0), _four(4.0);
	polynomial< bigcomplex > Q(P);
	bigcomplex x, x1(1.3, 0.31415926535897932384);
	int n = P.degree(), c, r;

	while (n > 2) {
		x = root(Q, x1);
		c = x.imag().is_approx_zero();
		if (c)
			x.assign(bigcomplex(x.real()));
		wurz[cnt].assign(x);
		cnt++;
		x.negate();
		r = R.deg;
		R1.set_degree(r + 1);
		R1.coeff[r+1].assign(R.coeff[r]);
		while (r > 0) {
			multiply(R1.coeff[r], x, R.coeff[r]);
			add(R1.coeff[r], R1.coeff[r], R.coeff[r-1]);
			r--;
		}
		multiply(R1.coeff[0], x, R.coeff[0]);
		R.assign(R1);
		n--;
		if (f && !c) {
			x.negate();
			conj(x, x);
			wurz[cnt].assign(x);
			cnt++;
			x.negate();
			r = R.deg;
			R1.set_degree(r + 1);
			R1.coeff[r+1].assign(R.coeff[r]);
			while (r > 0) {
				multiply(R1.coeff[r], x, R.coeff[r]);
				add(R1.coeff[r], R1.coeff[r], R.coeff[r-1]);
				r--;
			}
			multiply(R1.coeff[0], x, R.coeff[0]);
			R.assign(R1);
			n--;
		}
		divide(Q, P, R);
	}

	if (n == 2) {
		bigcomplex y;
		multiply(y, _four, Q.coeff[2]);
		multiply(y, y, Q.coeff[0]);
		square(x, Q.coeff[1]);
		subtract(x, x, y);
		sqrt(y, x);
		subtract(wurz[cnt], y, Q.coeff[1]);
		multiply(y, _two, Q.coeff[2]);
		divide(wurz[cnt], wurz[cnt], y);
		cnt++;
		sqrt(y, x);
		add(wurz[cnt], y, Q.coeff[1]);
		wurz[cnt].negate();
		multiply(y, _two, Q.coeff[2]);
		divide(wurz[cnt], wurz[cnt], y);
		cnt++;
		return;
	}
	if (n == 1) {
		negate(wurz[cnt], Q.coeff[0]);
		cnt++;
		return;
	}
}



//
// A root finding algorithm over C by H.Cohen
// MODIFICATION: The polynomial does not have to be squarefree.
//

void
roots(const base_polynomial< bigcomplex > & P, bigcomplex * wurz)
{
	debug_handler_l(DM_CP, "in member - function "
			"roots (const base_polynomial< bigcomplex > &, "
			"bigcomplex *)", DV_CP + 3);

	polynomial< bigcomplex > PP = P, SQRF, GCD;
	int n = P.degree(), i, f = 1, deg, count = 0;
	static bigfloat _two = 2.0, _four = 4.0;
	bigcomplex x;

	random_generator::seed(1);
	for (i = 0; i <= n; i++)
		if (PP.coeff[i].imag().length()) {
			f = 0;
			break;
		}
	if (PP.coeff[n] != bigcomplex(1.0)) {
		invert(x, PP.coeff[n]);
		for (i = 0; i <= n; i++)
			multiply(PP.coeff[i], PP.coeff[i], x);
	}
	if (n == 2) {
		bigcomplex y;
		multiply(y, _four, PP.coeff[2]);
		multiply(y, y, PP.coeff[0]);
		square(x, PP.coeff[1]);
		subtract(x, x, y);

		sqrt(y, x);
		subtract(wurz[count], y, PP.coeff[1]);
		multiply(y, _two, PP.coeff[2]);
		divide(wurz[count], wurz[count], y);
		count++;

		sqrt(y, x);
		add(wurz[count], y, PP.coeff[1]);
		wurz[count].negate();
		multiply(y, _two, PP.coeff[2]);
		divide(wurz[count], wurz[count], y);
	}
	else if (n == 1) {
		negate(wurz[count], PP.coeff[0]);
		divide(wurz[count], wurz[count], PP.coeff[1]);
	}
	else {
		SQRF = derivative(PP);
		GCD = gcd(PP, SQRF);
		divide(SQRF, PP, GCD);
		//********************
		if (SQRF.deg == 1) {
			negate(x, SQRF.coeff[0]);
			divide(x, x, SQRF.coeff[1]);
			for (i = 0; i < n; i++)
				wurz[i].assign(x);
			return;
		}
		//********************
		if (PP.deg == SQRF.deg)
			cohen(PP, wurz, f, count);
		else
			for (;;) {
				deg = SQRF.deg;
				n -= deg;
				if (deg == 2) {
					bigcomplex y;
					multiply(y, _four, SQRF.coeff[2]);
					multiply(y, y, SQRF.coeff[0]);
					square(x, SQRF.coeff[1]);
					subtract(x, x, y);

					sqrt(y, x);
					subtract(wurz[count], y, SQRF.coeff[1]);
					multiply(y, _two, SQRF.coeff[2]);
					divide(wurz[count], wurz[count], y);
					count++;

					sqrt(y, x);
					add(wurz[count], y, SQRF.coeff[1]);
					wurz[count].negate();
					multiply(y, _two, SQRF.coeff[2]);
					divide(wurz[count], wurz[count], y);
					count++;
				}
				else if (deg == 1) {
					negate(wurz[count], SQRF.coeff[0]);
					divide(wurz[count], wurz[count], SQRF.coeff[1]);
					count++;
				}
				else
					cohen(SQRF, wurz, f, count);
				if (n <= 0)
					break;
				divide(PP, PP, SQRF);
				SQRF = derivative(PP);
				GCD = gcd(PP, SQRF);
				divide(SQRF, PP, GCD);
			}
	}
}



bigcomplex
integral(const bigfloat & a, const bigfloat & b,
         const base_polynomial< bigcomplex > & P)
{
	debug_handler_l(DM_CP, "in friend - function "
			"integral (const bigfloat &, const bigfloat &, "
			"const base_polynomial< bigcomplex > &)", DV_CP + 1);

	int i = P.degree();

	if (i < 0) return bigcomplex(bigfloat(0.0));

	polynomial< bigcomplex > IntP;
	IntP.set_degree(i + 1);

	bigcomplex result;

	while (i > -1) {
		divide(IntP.coeff[i + 1], P[i], bigfloat(i + 1.0));
		i--;
	}

	IntP.coeff[0].assign(bigfloat(0.0));
	subtract(result, IntP(b), IntP(a));

	return result;
}



#undef DV_CP
#undef DM_CP
#undef LP_ERROR



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

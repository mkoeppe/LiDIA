// -*- C++ -*-
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

// Elliptic Curve Challenge Slave Program
//
// File: udigit128.h, gmp version
//
// $Id$


#ifndef LIDIA_HAVE_GMP_H
#define LIDIA_HAVE_GMP_H

#include	"LiDIA/kernel/gmp-includes.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigmod.h"

#include	<stream.h>
#include	<ctype.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// define this to the maximum number of bits you want to have

#ifndef C_BIT_SIZE
#define C_BIT_SIZE 128
#endif

//
// the number of mp_limbs required for representing a C_BIT_SIZE-bit number

const int C_LG = (C_BIT_SIZE / BITS_PER_MP_LIMB) +
((C_BIT_SIZE % BITS_PER_MP_LIMB) ? 1 : 0);

//
// the number of mp_limbs required for representing the product of
// 2 C_BIT_SIZE-bit numbers

const int C_LG2 = 2 * C_LG;

//
// rename mp_limb_t to digit

typedef mp_limb_t digit;

class udigit128;

extern udigit128 modulus; // modulus and normalised modulus
extern udigit128 n_modulus;
extern int n_mod_length; // length of normalised modulus
extern int toshift; // shiftting distance for normalization

//
// udigit128 represents fixed-size unsigned integers


class udigit128
{
private:

	digit m[C_LG];

	void clear(int s = 0, int e = C_LG)
	{
		for (int i = s; i < C_LG; i++)
			m[i] = 0;
	}

	void copy(const digit *p, int s = 0, int e = C_LG)
	{
		int i;

		for (i = s; i < e; i++)
			m[i] = p[i];
		if (e < C_LG)
			for (i = e; i < C_LG; i++)
				m[i] = 0;
	}

public:

	udigit128()
	{
		clear();
	}
	udigit128(digit x)
	{
		m[0] = x;
		clear(1);
	}
	udigit128(const char *s)
	{
		assign_string(s);
	}

	udigit128(const bigint & x)
	{
		char s[300];
		bigint_to_string(x, s); assign_string(s);
	}
	udigit128(const bigmod & x)
	{
		char s[300];
		bigint_to_string(x.mantissa(), s); assign_string(s);
	}

	udigit128(const udigit128 & x)
	{
		copy(x.m);
	}
	~udigit128() { }

	int length() const
	{
		for (int i = C_LG; i > 0; i--)
			if (m[i - 1])
				return i;
		return 0;
	}

	bool is_zero() const
	{
		return (length() == 0);
	}

	bool is_one() const
	{
		return (m[0] == 1) && (length() == 1);
	}

	unsigned long least_significant_digit()
	{
		return m[0];
	}

	void assign(const udigit128 & x)
	{
		copy(x.m);
	}



	udigit128 & operator = (const udigit128 & a)
	{
		copy(a.m);
		return *this;
	}

	udigit128 & operator = (const bigint & aa)
	{
		bigint a(aa);
		for (int i = 0; i < C_LG; i++) {
			if (a.is_zero())
				clear(i);
			else {
				m[i] = a.least_significant_digit();
				shift_right(a, a, 32);
			}
		}
		return *this;
	}

	udigit128 & operator = (unsigned long a)
	{
		m[0] = a;
		clear(1);
		return *this;
	}


	int compare(const udigit128 & y) const
	{
		return mpn_cmp(m, y.m, C_LG);
	}

	friend int cmp(const udigit128 & a, const udigit128 &b)
	{
		for (int i = C_LG-1; i >= 0; i--) {
			if (a.m[i] != b.m[i])
				if (a.m[i] > b.m[i])
					return 1;
				else
					return -1;
		}
		return 0;
	}

	friend bool operator == (const udigit128 &a, const udigit128 &b)
	{
		return (cmp(a, b) == 0);
	}

	friend bool operator != (const udigit128 &a, const udigit128 &b)
	{
		return (cmp(a, b) != 0);
	}

	friend bool operator > (const udigit128 &a, const udigit128 &b)
	{
		return (cmp(a, b) == 1);
	}

	void divide_by_2()
	{
		mpn_rshift(m, m, C_LG, 1);
	}


	//--------------------------------------------------------------
	// integer arithmetic
	//--------------------------------------------------------------

	//--------------------------------------------------------------
	friend void add(udigit128 &c, const udigit128 & x, const udigit128 & y)
	{
		if (&x == &y)
			mpn_lshift(c.m, x.m, C_LG, 1);
		else
			mpn_add_n(c.m, x.m, y.m, C_LG);
	}

	friend void add(udigit128 &r, const udigit128 & x, digit y)
	{
		mpn_add_1(r.m, x.m, C_LG, y);
	}

	//--------------------------------------------------------------
	friend void subtract(udigit128 & r, const udigit128 & x, const udigit128 & y)
	{
		mpn_sub_n(r.m, x.m, y.m, C_LG);
	}

	friend void subtract(udigit128 &r, const udigit128 & x, digit y)
	{
		mpn_sub_1(r.m, x.m, C_LG, y);
	}

	//--------------------------------------------------------------

	friend void multiply(udigit128 &r, const udigit128 & x,
			     const udigit128 & y)
	{
		digit accu[C_LG2]; mpn_mul(accu, x.m, C_LG, y.m, C_LG);
		r.copy(accu);
	}
	friend void multiply(udigit128 &r, const udigit128 & x, digit y)
	{
		mpn_mul_1(r.m, x.m, C_LG, y);
	}

	//--------------------------------------------------------------

	friend digit divide(udigit128 & q, const udigit128 & x, digit y)
	{
		return mpn_divmod_1(q.m, x.m, C_LG, y);
	}

	//---- special division function where 0 <= a/b < 2^32 ---------

	friend unsigned long divide(const udigit128 & a, const udigit128 & b)
	{
		register int n, alength = a.length(), blength = b.length();
		unsigned long h;

		if (cmp(a, b) <= 0)
			if (cmp(a, b) < 0)
				return 0;
			else
				return 1;

		count_leading_zeros(n, b.m[blength - 1]);
		digit paccu[alength], paccu2[alength+1];

		if (n != 0) {
			mpn_lshift(paccu, b.m, alength, n);
			mpn_lshift(paccu2, a.m, alength, n);
			h = paccu2[alength-1]/paccu[alength-1];
		}
		else {
			h = a.m[alength-1] / b.m[alength-1];
		}

		return h;
	}


	//--------------------------------------------------------------
	// modular arithmetic
	//--------------------------------------------------------------

	// multiplies (*this) by 2^n, so that the most significant bit
	// of (*this) is set. Returns the length of (*this) and the number
	// of bits n. The following stands:
	// new_this = old_this * 2^n && length(new_this) == length(old_this)
	//                           && length(new_this) <= C_LG
	int normalize(int & n)
	{
		register int l = length();
		count_leading_zeros(n, m[l - 1]);
		if (n != 0)
			mpn_lshift(m, m, l, n);
		return l;
	}

	//--------------------------------------------------------------
	friend void add_mod(udigit128 &r, const udigit128 & x, const udigit128 & y,
			    const udigit128 & p)
	{
		if (&x == &y)
			mpn_lshift(r.m, x.m, C_LG, 1);
		else
			mpn_add_n(r.m, x.m, y.m, C_LG);
		if (cmp(r, p) >= 0)
			mpn_sub_n(r.m, r.m, p.m, C_LG);
	}

	friend void add_mod(udigit128 &r, const udigit128 & x, digit y,
			    const udigit128 & p)
	{
		mpn_add_1(r.m, x.m, C_LG, y);
		if (cmp(r, p) >= 0)
			mpn_sub_n(r.m, r.m, p.m, C_LG);
	}

	//--------------------------------------------------------------
	friend void subtract_mod(udigit128 &r, const udigit128 & x,
				 const udigit128 & y, const udigit128 & p)
	{
		if (cmp(x, y) >= 0)
			mpn_sub_n(r.m, x.m, y.m, C_LG);
		else {
			if (&r != &x) {
				// look out, you overwrite x  !!
				mpn_sub_n(r.m, p.m, y.m, C_LG);
				mpn_add_n(r.m, r.m, x.m, C_LG);
			}
			else {
				mpn_add_n(r.m, r.m, p.m, C_LG);
				mpn_sub_n(r.m, r.m, y.m, C_LG);
			}
		}
		if (cmp(r, p) >= 0)
			mpn_sub_n(r.m, r.m, p.m, C_LG);
	}

	friend void subtract_mod(udigit128 &r, const udigit128 & x, digit & y,
				 const udigit128 & p)
	{
		mpn_sub_1(r.m, p.m, C_LG, y);
		mpn_add_n(r.m, r.m, x.m, C_LG);
		if (cmp(r, p) >= 0)
			mpn_sub_n(r.m, r.m, p.m, C_LG);
	}

	friend negate_mod(udigit128 &r, const udigit128 & x, const udigit128 & p)
	{
		if (x.is_zero())
			r.assign(x);
		else
			subtract(r, p, x);
	}

	//--------------------------------------------------------------

	friend void multiply_mod(udigit128 &r, const udigit128 & x,
				 const udigit128 & y, const udigit128 & p)
	{
		register int n, plength = p.length();
		count_leading_zeros(n, p.m[plength - 1]);
		digit paccu[C_LG+1];
		if (n != 0)
			mpn_lshift(paccu, p.m, plength, n);
		else {
			int i;
			for (i = 0; i < plength; i++) {
				paccu[i] = p.m[i];
			}
		}


		// compute x * y in the C_LG2-digit-buffer accu
		digit accu[C_LG2+2];
		for (int i = 0; i < C_LG2; i++) accu[i] = 0;
		mpn_mul_n(accu, x.m, y.m, C_LG);
		if (n != 0)
			mpn_lshift(accu, accu, C_LG2, n);

		// compute the length of x * y
		int alength = C_LG2;
		while (alength && accu[alength - 1] == 0) alength--;

		// if length(accu) < length(paccu), then (x * y) % pn == (x * y)
		if (alength < plength)
			r.copy(accu);
		else {
			mpn_divrem(accu + plength, 0, accu, alength, paccu, plength);
			r.copy(accu, 0, plength);
		}

		// shift down
		if (n)
			mpn_rshift(r.m, r.m, C_LG, n);
	}

	friend void invert_mod(udigit128 &r, const udigit128 & x,
			       const udigit128 & p)
	{
		digit tg[C_LG], tx[C_LG], tp[C_LG];
		udigit128 h;

		h = x;
		add(r, x, p);
		int xl = r.length();

		for (int i = 0; i < C_LG; i++) {
			tg[i] = 0;
			tx[i] = r.m[i];
			tp[i] = p.m[i];
		}
		r.clear();

		mpn_gcdext(tg, r.m, tx, xl, tp, p.length());

		multiply_mod(h, r, h); // this test is necessary due to a gmp bug.
		if (!h.is_one())
			subtract(r, p, r);
	}

	friend void divide_by_2_mod(udigit128 &r, const udigit128 &a,
				    const udigit128 &p)
	{
		if (!(a.m[0] & 1))
			mpn_rshift(r.m, a.m, C_LG, 1);
		else {
			add(r, a, p);
			r.divide_by_2();
		}
	}




	//--------------------------------------------------------------
	// input / output
	//--------------------------------------------------------------

	void assign_string(const char *s)
	{
		// eat sign and find out the base
		const char *p = s;
		int base = 10, dig;
		while (*p == ' ' || *p == '\t')
			++p;
		if (*p == '+')
			++p;
		if (*p == '0') {
			base = 8, ++p;
			if (tolower(*p) == 'x')
				base = 16, ++p;
		}

		// assign the value of s to (*this)
		clear();
		for (;; ++p) {
			dig = (base == 16) ? tolower(*p) : *p;
			if (dig >= 'a' && dig <= 'f')
				dig -= 'a', dig += 10;
			else
				dig -= '0';
			if (dig >= 0 && dig < base) {
				mpn_mul_1(m, m, C_LG, base);
				mpn_add_1(m, m, C_LG, dig);
			}
			else
				break;
		}
	}

	char *as_string(int base = 10) const
	{
		static char hex_digit[] = {'0', '1', '2', '3', '4', '5', '6', '7',
					   '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
		static char buffer[BITS_PER_MP_LIMB];
		int i = BITS_PER_MP_LIMB - 1;
		buffer[i] = '\0';
		if (is_zero()) {
			buffer[--i] = '0';
			if (base == 16) {
				buffer[--i] = 'x';
				buffer[--i] = '0';
			}
		}
		else {
			udigit128 n(*this), q;
			digit d;
			while (!n.is_zero()) {
				d = divide(n, n, base);
				if (base == 16)
					buffer[--i] = hex_digit[static_cast<int>(d)];
				else
					buffer[--i] = '0' + static_cast<int>(d);
			}
		}
		if (base == 8)
			buffer[--i] = '0';
		else if (base == 16) {
			buffer[--i] = 'x';
			buffer[--i] = '0';
		}
		return buffer + i;
	}

	friend void udigit_to_string(const udigit128 &x, char *s)
	{
		strcpy(s, x.as_string());
	}

	friend std::ostream & operator << (std::ostream & out, const udigit128 & x)
	{
		out << x.as_string();
		return out;
	}

	void print() const
	{
		std::cout << "[ ";
		for (int i = 0; i < length(); i++)
			std::cout << m[i] << " ";
		std::cout << " ]" << std::endl;
	}



	//--------------------------------------------------------------------
	// and now specialised modular function which use precomputed
	// normalized modulus. Note that these variables have to set correctly.


	friend void add_mod(udigit128 &r, const udigit128 & x, const udigit128 & y)
	{
		add_mod(r, x, y, modulus);
	}

	friend void add_mod(udigit128 &r, const udigit128 & x, digit y)
	{
		add_mod(r, x, y, modulus);
	}

	friend void subtract_mod(udigit128 &r, const udigit128 & x,
				 const udigit128 & y)
	{
		subtract_mod(r, x, y, modulus);
	}

	friend void subtract_mod(udigit128 &r, const udigit128 & x, digit y)
	{
		subtract_mod(r, x, y, modulus);
	}

	friend void negate_mod(udigit128 &r, const udigit128 & a)
	{
		negate_mod(r, a, modulus);
	}

	friend void divide_by_2_mod(udigit128 &r, const udigit128 & a)
	{
		divide_by_2_mod(r, a, modulus);
	}

	//--------------------------------------------------------------
	// compute (x * y) mod pn. pn is assumed to be normalized, i.e
	// its most significant bit is set. The result is shifted down
	// by n bits. This function is used to avoid normalizing when
	// computing mod p. Instead, we compute pn := p * 2^n and
	// shift down each time appropriately (by n). This function
	// assumes length(pn) == plength and n < BITS_PER_MP_LIMB


	friend void multiply_mod(udigit128 & r, const udigit128 & x,
				 const udigit128 & y)
	{
		// compute x * y * 2^n in the C_LG2-digit-buffer accu
		// since x < p and y < p the result fits into 2 * C_LG words

		digit accu[C_LG2];
		for (int i = 0; i < C_LG2; i++)
			accu[i] = 0;
		mpn_mul_n(accu, x.m, y.m, C_LG);
		if (toshift != 0)
			mpn_lshift(accu, accu, C_LG2, toshift);

		// compute the length of x * y: this may be less than 2 * C_LG
		int alength = C_LG2;
		while (alength && accu[alength - 1] == 0)
			alength--;

		// if length(accu) < length(pn), then (x * y) % pn == (x * y)
		// have to handle this because of the mpn_divrem preconditions
		if (alength < n_mod_length)
			r.copy(accu);
		else {
			mpn_divrem(accu + n_mod_length, 0, accu, alength, n_modulus.m,
				   n_mod_length);
			r.copy(accu, 0, n_mod_length);
		}

		// shift down
		// x * y * 2^n = q * pn + r = q * p * 2^n + r ==>
		// x * y = q * p + r / 2^n
		if (toshift)
			mpn_rshift(r.m, r.m, C_LG, toshift);
	}

	friend void square_mod(udigit128 &r, const udigit128 & a)
	{
		multiply_mod(r, a, a);
	}

	friend void invert_mod(udigit128 &r, const udigit128 & x)
	{
		invert_mod(r, x, modulus);
	}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

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
//	$Id: bigint.cc,v 2.9 2004/06/15 10:19:48 lidiaadm Exp $
//
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/random_generator.h"

#include <vector>
#include <cassert>
#include <cstring>

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#if !INLINE_INTERFACE
# define inline
# include	"LiDIA/kernel/bigint_interface.h"
# undef inline
#endif



xdouble
bigint::xdbl () const
{
	long l = I._mp_size;

	l = ((l < 0) ? -l : l);
	xdouble d = 0.0;
	if (l) {
		int i = 1;

		d = static_cast<double>(I._mp_d[0]);
		xdouble base = bigint::radix();
		xdouble power = base;
		while (i < l) {
			d += I._mp_d[i] * power;
			power *= base;
			i++;
		}
		if (I._mp_size < 0) {
			d = -d;
		}
	}
	return d;
}



void
power (bigint & c, const bigint & a, const bigint & b)
{
	bigint exponent, multiplier;

	if ((a.I._mp_size == 1 || a.I._mp_size == -1) && a.I._mp_d[0] == 1) {
		// a = +-1->a^(|b|)
		if (b.is_odd()) {
			c = a;
		}
		else {
			c.assign_one();
		}
	}
	else if (b.is_negative()) {
		// b < 0, a != +-1->0
		c.assign_zero();
	}
	else if (b.is_zero()) {
		// b == 0->1
		c.assign_one();
	}
	else if (b.is_one()) {
		// b == 1->a
		c.assign(a);
	}
	else {
		exponent.assign(b);

		multiplier = a;
		c = a;

		// Note: The bitlength is at least 2.
		lidia_size_t length = exponent.bit_length()-2;
		bigint tmp;

		shift_left(tmp, 1, length);

		while (tmp.is_gt_zero()) {
			square(c, c);
			if (!((exponent&tmp).is_zero())) {
				multiply(c, c, multiplier);
			}
			tmp.divide_by_2();
		}
	}
}



void
power (bigint & c, const bigint & a, long i)
{
	if ((a.I._mp_size == 1 || a.I._mp_size == -1) && a.I._mp_d[0] == 1) {
		// a = +-1->a^(|b|)
		if (i&1) {
			c = a;
		}
		else {
			c.assign_one();
		}
	}
	else if (i < 0) {
		// i < 0, a != +-1->0
		c.assign_zero();
	}
	else if (i == 0) {
		// i == 0->1
		c.assign_one();
	}
	else {
		mpz_pow_ui(&c.I, &a.I, i);
	}
}



//
// random numbers
//

static gmp_randstate_t *
get_randstate ()
{
	static gmp_randstate_t rstate;
	static bool initialized = false;

	if (!initialized) {
		gmp_randinit(rstate, GMP_RAND_ALG_DEFAULT, 32L);
		initialized = true;
	}
	return &rstate;
}



void
bigint::seed ()
{
	bigint_rep_t	tmp;
	unsigned long	hi, lo;
	random_generator	rg;

	rg >> hi >> lo;

	mpz_init_set_ui(&tmp, hi);
	mpz_mul_2exp(&tmp, &tmp, BITS_PER_LONG);
	mpz_add_ui(&tmp, &tmp, lo);
	gmp_randseed(*get_randstate(), &tmp);
	mpz_clear(&tmp);

	is_seeded = true;
}



void
bigint::seed (const bigint & a)
{
	bigint_rep_t	tmp;

	mpz_init_set(&tmp, &a.I);
	gmp_randseed(*get_randstate(), &tmp);
	mpz_clear(&tmp);

	is_seeded = true;
}



void
bigint::randomize (const bigint & a)
{
	if (!is_seeded) {
		seed();
	}
	if (a.is_zero()) {
		lidia_error_handler("bigint", "Bound must not be equal to zero.");
	}
	else {
		bigint tmp;

		mpz_urandomb(&tmp.I, *get_randstate(), mpz_sizeinbase(&a.I, 2));

		if (a.is_lt_zero()) {
			if (tmp <= a)
				remainder(tmp, tmp, a);
		}
		else if (tmp >= a)
			remainder(tmp, tmp, a);
		swap(tmp);
	}
}



//
// input / output
//

int
string_to_bigint (const char *s, bigint & a)
{
        assert(s);

        char const* ss = s;
        if(ss[0] == '+' || ss[0] == '-') {
	        ss += 1;
        }

	if(ss[0] == '0' && ss[1] == '\0') {
	    a.assign_zero();
	    return std::strlen(s);
	}
	// now any leading '0' indicates a base != 10

        int base;
	if(ss[0] != '0') {
	        base = 10;
	}
	else if(ss[1] == 'x') {
	        base = 16;
		ss += 2;
	}
	else {
	        base = 8;
		ss += 1;
	}

	if (ss[0] != '\0' && !mpz_set_str(&a.I, const_cast<char*>(ss), base)) {
	        if(s[0] == '-') {
		        a.negate();
		}
		return std::strlen(s);
	}
	else {
		return 0;
	}
}


std::istream &
operator >> (std::istream & in, bigint & a)
{
	a.scan(in);
	return in;
}



std::ostream &
operator << (std::ostream & out, const bigint & a)
{
    if(out.good()) 
    {
	std::ios::fmtflags streamflags = out.flags();
	int base = ((streamflags & std::ios::hex) ? 16 :
		    (streamflags & std::ios::oct) ? 8 : 10);
	bool printbase = streamflags & std::ios::showbase;
	
	typedef std::vector<char> CharVector;
	typedef CharVector::size_type VectorSize;
	VectorSize l = mpz_sizeinbase(&(a.bigint_rep()), base);
	CharVector v(l + 10); // allow for minus sign, base indicator,
	                      // trailing '\0' etc. 

	if(printbase && base != 10) {
	    l = (base == 16) ? 2 : 1;
	}
	else {
	    l = 0;
	}
	mpz_get_str(&v[l], base, &(a.bigint_rep()));
	if(printbase && base != 10) {
	    if(a.is_negative()) {
		v[0] = '-';
		l = 1;
	    }
	    else {
		l = 0;
	    }
	    v[l] = '0';
	    if(base == 16) {
		v[l+1] = 'x';
	    }
	}
	a.print(out, &v[0]);

// Changed by G.A. - we must return in any case, otherwise gcc 4.1.2
// misunderstands this trying to optimize and the program crashes!
	//return out;
    } // if(out.good()) 
    return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

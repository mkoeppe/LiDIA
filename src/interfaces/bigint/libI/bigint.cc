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

#include        <iostream>
#include        <vector>
#include        <cassert>
#include        <cstring>

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#if !INLINE_INTERFACE
# define inline
# include	"LiDIA/kernel/bigint_interface.h"
# undef inline
#endif



double
bigint::dbl () const
{
	int l = I.length;
	double d = 0.0;

	if (l != 0) {
		int i = 1;
		d = I.vec[0];
		double base = bigint::radix();
		double power = base;
		while (i < l) {
			d += I.vec[i] * power;
			power *= base;
			i++;
		}
		if (I.sign == MINUS) {
			d = -d;
		}
	}
	return d;
}



xdouble
bigint::xdbl () const
{
	int l = I.length;
	xdouble d = 0.0;

	if (l) {
		int i = 1;
		d = static_cast<double>(I.vec[0]);
		xdouble base = bigint::radix();
		xdouble power = base;
		while (i < l) {
			d += I.vec[i] * power;
			power *= base;
			i++;
		}
		if (I.sign == 1) {
			d = -d;
		}
	}
	return d;
}



void
power (bigint & c, const bigint & a, const bigint & b)
{
	bigint exponent, multiplier;

	if ((a.I.length == 1) && (*(a.I.vec) == 1)) {
		// a = +-1 -> a^(|b|)
		if (b.is_odd()) {
			c = a;
		}
		else {
			c.assign_one();
		}
	}
	else if (b.is_negative()) {		// b < 0, a != +-1->0
		c.assign_zero();
	}
	else if (b.is_zero()) {			// b == 0->1
		c.assign_one();
	}
	else if (b.is_one()) {			// b == 1->a
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
			if (!((exponent&tmp).is_zero()))
				multiply(c, c, multiplier);
			tmp.divide_by_2();
		}
	}
}



void
power (bigint & c, const bigint & a, long i)
{
	if ((a.I.length == 1) && (*(a.I.vec) == 1)) {
		// a = +-1 -> a^(|i|)
		if (i&1) {
			c = a;
		}
		else {
			c.assign_one();
		}
	}
	else if (i < 0) {			// i < 0, a != +-1->0
		c.assign_zero();
	}
	else if (i == 0) {			// i == 0->1
		c.assign_one();
	}
	else {
		IasIpowD(&c.I, &a.I, i);
	}
}



void
bigint::swap (bigint & b)
{
	int tmp;
	DigitType *tmpp;

	tmp = I.sign;
	I.sign = b.I.sign;
	b.I.sign = tmp;

	tmp = I.length;
	I.length = b.I.length;
	b.I.length = tmp;

	tmp = I.maxlength;
	I.maxlength = b.I.maxlength;
	b.I.maxlength = tmp;

	tmpp = I.vec;
	I.vec = b.I.vec;
	b.I.vec = tmpp;
}



//
// random numbers
//

void
bigint::seed ()
{
	unsigned long		lo;
	random_generator	rg;

	rg >> lo;

	IseedD(lo);

	is_seeded = true;
}



void
bigint::seed (const bigint & a)
{
	IseedD(a.least_significant_digit());

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
		IasrandomI(&I, &a.I);
	}
}



//
// input / output
//

int
string_to_bigint (const char *s, bigint & a)
{

    assert(s);

    int pos = 0;
    bool negative = false;
    if(s[pos] == '+' || s[pos] == '-') {
	negative = (s[pos] == '-');
	pos += 1;
    }
    
    if(s[pos] == '0' && s[pos+1] == '\0') {
	a.assign_zero();
	return pos + 1;
    } 
    // now any leading '0' indicates a base != 10

    long base;
    if((s[pos] == '0') && (s[pos + 1] == 'x' || s[pos + 1] == 'X')) {
	base = 16;
	pos += 2;
    }
    else if(s[pos] == '0') {
	base = 8;
	pos += 1;
    }
    else {
	base = 10;
    }

    if(s[pos] == '+' || s[pos] == '-') {
	negative = (negative != (s[pos] == '-'));
	pos += 1;
    }
    
    a.assign_zero();
    bool finished = false;
    do {
	switch(s[pos++]) {
	    case '0':
		if(a.is_zero()) {
		    finished = true;
		}
		else {
		    multiply(a, a, base);
		}
		break;
	    case '1':
		multiply(a, a, base);
		a.inc();
		break;
	    case '2':
		multiply(a, a, base);
		add(a, a, 2L);
		break;
	    case '3':
		multiply(a, a, base);
		add(a, a, 3L);
		break;
	    case '4':
		multiply(a, a, base);
		add(a, a, 4L);
		break;
	    case '5':
		multiply(a, a, base);
		add(a, a, 5L);
		break;
	    case '6':
		multiply(a, a, base);
		add(a, a, 6L);
		break;
	    case '7':
		multiply(a, a, base);
		add(a, a, 7L);
		break;
	    case '8':
		if(base <= 8) {
		    finished = true;
		}
		else {
		    multiply(a, a, base);
		    add(a, a, 8L);
		}
		break;
	    case '9':
		if(base <= 8) {
		    finished = true;
		}
		else {
		    multiply(a, a, base);
		    add(a, a, 9L);
		}
		break;
	    case 'a':
	    case 'A':
		if(base <= 10) {
		    finished = true;
		}
		else {
		    shift_left(a, a, 4);
		    add(a, a, 10L);
		}
		break;
	    case 'b':
	    case 'B':
		if(base <= 10) {
		    finished = true;
		}
		else {
		    shift_left(a, a, 4);
		    add(a, a, 11L);
		}
		break;
	    case 'c':
	    case 'C':
		if(base <= 10) {
		    finished = true;
		}
		else {
		    shift_left(a, a, 4);
		    add(a, a, 12L);
		}
		break;
	    case 'd':
	    case 'D':
		if(base <= 10) {
		    finished = true;
		}
		else {
		    shift_left(a, a, 4);
		    add(a, a, 13L);
		}
		break;
	    case 'e':
	    case 'E':
		if(base <= 10) {
		    finished = true;
		}
		else {
		    shift_left(a, a, 4);
		    add(a, a, 14L);
		}
		break;
	    case 'f':
	    case 'F':
		if(base <= 10) {
		    finished = true;
		}
		else {
		    shift_left(a, a, 4);
		    add(a, a, 15L);
		}
		break;
	    default:
		finished = true;
	}
    } while(!finished);

    if(negative) {
	a.negate();
    }

    return pos;
}



std::istream &
operator >> (std::istream & in, bigint & a)
{
	a.scan (in);
	return (in);
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
	VectorSize l = Ilog(&a.I) + 2;
	l /= (base == 16) ? 4 : 3;
	l += 10;

	CharVector v;
	v.resize(l);
	v[--l] = '\0';

	if(a.is_zero()) {
	    v[--l] = '0';
	}
	else {
	    char const digit[] = "0123456789abcdef";

	    bigint tmp(a);
	    tmp.abs();
	    while(!tmp.is_zero()) {
		long rem;
		div_rem(tmp, rem, tmp, base);
		v[--l] =  digit[rem];
	    }
	}

	if(printbase) {
	    if(base == 16) {
		v[--l] = 'x';
		v[--l] = '0';
	    }
	    else if (base == 8) {
		v[--l] = '0';
	    }
	}

	if(a.is_negative()) {
	    v[--l] = '-';
	}

	// if l < 0 then the program may already hve crashed, but anyway...
	assert(l >= 0);
	a.print (out, &v[l]);
// Changed by G.A. - we must return in any case, otherwise gcc 4.1.2
// misunderstands this trying to optimize and the program crashes!
//	return out;
    } // if(out.good()) 
    return out;
}



#ifdef C_STDIO

//
// using fread/fwrite
//


void
bigint::read_from_file (FILE * fp)
{
	int soi, soD;

	soi = sizeof(int);
	soD = sizeof(DigitType);
	if (feof(fp)) {
		return;
	}
	delDigitVec(I.vec, I.maxlength);
	fread ((char *)&I.maxlength, soi, 3, fp);
	I.maxlength = I.length;
	I.vec = (DigitType *) newDigitVec (&I.maxlength);
	fread ((char *)I.vec , soD, I.length, fp);
}



void bigint::write_to_file (FILE * fp)
{
	int soi, soD;

	soi = sizeof(int);
	soD = sizeof(DigitType);
	fwrite ((char *)&I.maxlength, soi, 3, fp);
	fwrite ((char *)I.vec, soD, I.length, fp);
}



//
// using fscanf/fprintf
//

void
bigint::scan_from_file (FILE * fp)
{
	fscanI(fp, &I);
}



void
bigint::print_to_file (FILE * fp)
{
	fprintI(fp, &I);
}



#endif	// C_STDIO



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif

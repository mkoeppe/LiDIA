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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/random_generator.h"
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#if !INLINE_INTERFACE
# define inline
# include	"LiDIA/kernel/bigint_interface.h"
# undef inline
#endif



//
// assigners
//

void
bigint::assign (double d)
{
	double mantisse;
	int exponent;
	unsigned long lmantisse;
	int sign = 0, count = 0;

	if (d < 0.0) {
		sign = 1;
		d = -d;
	}

	if (d >= 0 && d < 0.5) {
		I = 0;
	}
	else {
		I = 0;
		mantisse = frexp(d, &exponent);
		while (mantisse != 0.0) {
			I <<= BETA;
			mantisse = std::ldexp(mantisse, BETA);
			lmantisse = static_cast<unsigned long>(mantisse);
			mantisse -= lmantisse;
			I += Natural(lmantisse);
			count++;
		}
		exponent -= count * BETA;
		if (exponent >= 0)
			I <<= exponent;
		else
			I >>= (-exponent);
	}
	if (sign)
		I = -I;
}



//
// conversion
//

double
dbl (const bigint & a)
{
	int l = a.I.length();
	double d = 0.0;
	if (!a.is_zero()) {
		Integer z = a.I;
		int i = 1;
		d = z.lowest();
		double base = bigint::radix();
		double power = base;
		while (i < l) {
			z >>= BETA;
			d += z.lowest() * power;
			power *= base;
			i++;
		}
		if (a.is_negative())
			d = -d;
	}
	return d;
}



xdouble
xdbl (const bigint & a)
{
	int l = a.I.length();
	xdouble d = 0.0;
	if (a.is_zero()) {
		Integer z = a.I;
		int i = 1;
		d = z.lowest();
		xdouble base = bigint::radix();
		xdouble power = base;
		while (i < l) {
			z >>= BETA;
			d += z.lowest() * power;
			power *= base;
			i++;
		}
		if (a.is_negative())
			d = -d;
	}
	return d;
}



//
//
//

void
power (bigint & c, const bigint & a, const bigint & b)
{
	if ((a.I == 1) || (a.I == -1)) {
		if (b.is_odd()) {
			c.I = a.I;
		}
		else {
			c.assign_one();
		}
	}
	else if (b.is_negative()) {
		c.assign_zero();
	}
	else if (b.is_zero()) {
		c.assign_one();
	}
	else {
		c.I = pow(a.I, b.I);
	}
}



void
power (bigint & c, const bigint & a, long i)
{
	if ((a.I == 1) || (a.I == -1)) {
		if (i&1) {
			c.I = a.I;
		}
		else {
			c.assign_one();
		}
	}
	else if (i < 0) {
		c.assign_zero();
	}
	else if (i == 0) {
		c.assign_one();
	}
	else {
		c.I = pow(a.I, i);
	}
}



bigint
xgcd (bigint & u, bigint & v, const bigint & a, const bigint & b)
{
	if (a.I == 0) {
		u.I = 0;
		v = b.sign();
		return abs(b);
	}
	if (b.I == 0) {
		v.I = 0;
		u = a.sign();
		return abs(a);
	}
	if (abs(a.I) == abs(b.I)) {
		u.I = 0;
		v = b.sign();
		return abs(b);
	}

	Integer g;
	gcd(a.I, b.I, u.I, v.I, g);

	if (g < 0) {
		g = -g;
		u.I = - u.I;
		v.I = - v.I;
	}

#if 0
	if (abs(2*u.I*g) > abs(b.I)) {
		if (u.sign() == b.sign()) {
			u.I = u.I - exquo(b.I,g);
			v.I = v.I + exquo(a.I,g);
		}
		else {
			u.I = u.I + exquo(b.I,g);
			v.I = v.I - exquo(a.I,g);
		}
	}
	else if (abs(2*v.I*g) > abs(a.I)) {
		if (v.sign() == a.sign()) {
			u.I = u.I + exquo(b.I,g);
			v.I = v.I - exquo(a.I,g);
		}
		else {
			u.I = u.I - exquo(b.I,g);
			v.I = v.I + exquo(a.I,g);
		}
	}
#endif

	return bigint(g);
}



//
// random numbers
//

void
bigint::seed ()
{
	random_generator rg;
	unsigned long	aa;

	rg >> aa;
	srand(aa);

	is_seeded = true;
}



void
bigint::seed (const bigint & a)
{
	random_generator rg;
	unsigned long aa;

	rg >> aa;
	srand(a.least_significant_digit() ^ aa);

	is_seeded = true;
}



void
bigint::randomize (const bigint & a)
{
	if (!is_seeded) {
		seed();
	}
	Natural y;
	y.rand(log2(a.I)+1);
	Integer x = y;
	if (a.is_negative()) {
		x %= -a.I;
		x = -x;
	}
	else
		x %= a.I;
	*this = x;
}



//
// input / output
//

std::istream &
operator >> (std::istream & in, bigint & a)
{
	a.scan (in);
	return (in);
}



std::ostream &
operator << (std::ostream & out, const bigint & a)
{
	out << a.I;
	return out;
}



#ifdef C_STDIO

//
// using fread/fwrite
//

void
bigint::read_from_file (FILE * fp)
{
	scan_from_file(fp);
	int c = getc(fp);
	if (c != EOF && c != '\n') {
		ungetc(c, fp);
	}
}



void
bigint::write_to_file (FILE * fp)
{
	print_to_file(fp);
	putc('\n', fp);
}



//
// using fscanf/fprintf
//

void bigint::scan_from_file(FILE * fp)
{
	int alloc = 32; // must be > 0
	char* buffer = 0;
	bigint::allocate(buffer, static_cast<int>(0), alloc);
	int index = 0;

	int c;
	while ((c = getc(fp)) != EOF) {
		if (!isspace(c)) {
			break;
		}
	}

	switch(c) {
		//case EOF: lidia_error_handler("bigint", "scan_from_file(...)::end of file");
	case '-':
		bigint::append_char(buffer, alloc, index++, c);
	case '+':
		break;
	default :
		ungetc(c, fp);
		break;
	}

	while ((c = getc(fp)) != EOF) {
		if (isdigit(c)) {
			bigint::append_char(buffer, alloc, index++, c);
		}
		else {
			ungetc(c, fp);
			break;
		}
	}

	bigint::append_char(buffer, alloc, index, '\0');
	string_to_bigint(buffer, *this);

	delete[] buffer;
}



void
bigint::print_to_file (FILE * fp)
{
	char* s = new char[static_cast<int>(bit_length()*std::log(2.0)/std::log(10.0) + 10)];
	Itoa(I, s, 10);
	fprintf(fp, "%s", s);
	delete[] s;
}

#endif	// C_STDIO



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
